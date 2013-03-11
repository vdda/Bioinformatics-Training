#!/bin/bash

# A simple tool that runs HGAP using 1.4 software.
# Must be run from the same directory as your input.xml
# Directory must be writeable.
# Must have your smrtpipe environment setup.

QUEUE='secondary'
DEBUG=0
RUN_LOCAL=0
SEYMOUR_HOME=${SEYMOUR_HOME:?"SMRT Pipe environment not detected."}
RunCA=${RunCA:-'runCA'}

srcdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
qsub=$(which qsub)

usage() {
    echo -e "USAGE: $(basename $0) [params] <input.xml>\n"          \
            "-p    Optional path to a preassembler params file\n"   \
            "-r    Optional path to a resequencing params file\n"   \
            "-s    Optional path to a celera-assembler spec file\n" \
            "-x    Override default options to smrtpipe\n"          \
            "-l    Run everything locally (e.g., no cluster)\n"     \
            "-d    Keep all the intermediate files and output\n"    \
            "      some statistics at the end\n"
    exit 1
}

cbase() {
    if gzip -t $1 >& /dev/null
    then
        eval $2=$(gunzip -c $1 | perl -ne 'chomp; $l+=length($_) if !/^>/;END{print "$l\n"}')
    else
        eval $2=$(fastalength $1 | cut -d' ' -f1 | awk '{t+=$1}END{print t}')
    fi
}

n50() {
    eval $3=$(fastalength $1 | cut -d' ' -f1 | sort -nr | awk -v corrc=$2 '{t+=$1;if (t >= corrc / 2){print $1; exit;}}')
}

stats() {
    cbase data/filtered_subreads.fasta rawbc 
    cbase filtered_longreads.fasta seedc 
    cbase data/corrected.fasta corrc 
    cbase data/consensus.fasta.gz asmbc
    n50 data/corrected.fasta $corrc rfift 
    nread=$(grep -c '>' data/corrected.fasta)
    bloss=$(echo "scale=3;1-${corrc}/${seedc}" | bc)
    mread=$(echo "$corrc / $nread" | bc)
    ncont=$(gunzip -c data/consensus.fasta.gz | grep -c '>')

    printf "%-15s $rawbc\n" "Raw Bases:"
    printf "%-15s $seedc\n" "Seed Bases:"
    printf "%-15s $corrc\n" "Corr Bases:"
    printf "%-15s $bloss\n" "Base Loss:"
    printf "%-15s $nread\n" "# reads:"
    printf "%-15s $mread\n" "Mean Read:"
    printf "%-15s $rfift\n" "Read N50:"
    printf "%-15s $ncont\n" "# Contigs:"
    printf "%-15s $asmbc\n" "Asmbl Bases:"
}

debug() {
    if [ $DEBUG -eq 1 ]
    then
        echo $1 >&2
    fi
}

timeit() {
    desc=$1
    shift
    echo "Running $desc ..." >&2 
    debug "$@"
    
    cmd="$@"
    /usr/bin/time -f"${desc} time: %E real, %U user, %S sys" $cmd
    ec=$?
    if [ $ec -gt 0 ]
    then
        echo "$desc failed" >&2
        exit $ec
    fi 
}

if [ $# -lt 1 -o "$1" == "--help" ]
then
    usage
fi

while getopts ":p:r:s:x:dl" opt
do
    case $opt in
      p)
        p_preasm="$OPTARG"
        ;;
      r)
        p_reseq="$OPTARG"
        ;;
      s)
        p_caspec="$OPTARG"
        ;;
      x)
        x_opts="$OPTARG"
        ;;
      l)
        RUN_LOCAL=1
        ;;
      d)
        DEBUG=1
        ;;
      \?)
        echo "Invalid option: -$OPTARG" >&2
        usage
        ;;
      :)
        echo "Option -$OPTARG requires an argument." >&2
        usage
        ;;
    esac
done

input=${@:$OPTIND:1}

p_preasm=${p_preasm="${srcdir}/params_preasm.xml"}
p_reseq=${p_reseq="${srcdir}/params_reseq.xml"}

x_opts_def="-D MAX_THREADS=30 -D HEARTBEAT_FREQ=-1"
[ $RUN_LOCAL -eq 0 ] && x_opts_def="${x_opts_def} --distribute"
x_opts=${x_opts=${x_opts_def}}

maybe_qsub=
if [[ $RUN_LOCAL -eq 0 && -x "$qsub" ]]
then
    maybe_qsub="${qsub} -cwd -sync y -S /bin/bash -V -q ${QUEUE} -N CA -o ./CA.err -b y -j y -pe smp 15 "
fi

ca_opts=
if [ \! -z $p_caspec ]
then
    ca_opts="-s $p_caspec"
else
    if [ -e "${srcdir}/ca_default.spec" ]
    then 
        ca_opts="-s ${srcdir}/ca_default.spec"
    fi
fi

[ $DEBUG -eq 1 ] && x_opts="--debug ${x_opts}"

debug "p_preasm = ${p_preasm}"
debug "p_reseq = ${p_reseq}"
debug "input = ${input}"
debug "ca_opts = ${ca_opts}"
debug "x_opts = ${x_opts}"
debug "qsub = ${maybe_qsub}"
debug "run_local = ${RUN_LOCAL}"

if [ -z "$input" ]
then
    echo "INVALID cmd line: Must provide an input.xml."
    exit 1
fi

timeit "PreAssembler" "smrtpipe.py ${x_opts} --params=${p_preasm} xml:${input}" > /dev/null
if [ \! -s data/corrected.fasta ] 
then
    echo "No corrected reads, consider lowering minLongReadLength in $p_preasm" >&2 
    exit 1
fi
timeit "CA prep" "fastqToCA -technology sanger -type sanger -reads data/corrected.fastq -libraryname reads" > reads.frg
timeit "CA" "${maybe_qsub} ${RunCA} reads.frg -d assembly -p assembly ${ca_opts}"
ln -s assembly/9-terminator/assembly.scf.fasta reference.fasta
timeit "Create Reference Repository" "referenceUploader -p. -f reference.fasta -c -n reference"
timeit "Resequencing" "smrtpipe.py ${x_opts} --params=${p_reseq} xml:${input}" > /dev/null
[ $DEBUG -eq 1 ] && stats | tee stats.txt
