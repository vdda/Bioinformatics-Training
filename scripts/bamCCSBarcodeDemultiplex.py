from pbcore.io import BamReader
from pbcore.io.BarcodeH5Reader import BarcodeH5Fofn
from pbcore.io.FastqIO import FastqWriter
import os

def main(parser):

    args = parser.parse_args()

    bam    = BamReader(args.ccsBAM)
    bcFofn = BarcodeH5Fofn(args.barcodeFofn)

    oFiles =  { bc:FastqWriter('{dir}/{bc}.fastq'.format(dir=args.outDir,bc=bc)) for bc in bcFofn.barcodeLabels }
    for rec in bam:
        try:
            lZmw = bcFofn.labeledZmwFromName(rec.readName)
        except KeyError:
            #catch zmws with no barcode and skip
            continue
        if       rec.readScore     >= args.minPredictedAccuracy \
             and lZmw.averageScore >= args.minAvgBarcodeScore \
             and rec.numPasses     >= args.minNumPasses:
            header = rec.readName
            if args.extendedHeader:
                header +=  ' predictedAccuracy={predAcc} numPasses={numPasses} barcodeScore={bcScore}'\
                           .format(predAcc=rec.readScore, numPasses=rec.numPasses, bcScore=lZmw.averageScore)
            qual = [ ord(q)-33  for q in rec.peer.qual ]
            writer = oFiles[bcFofn.barcodeLabels[lZmw.bestIdx]]
            writer.writeRecord(header, rec.read(aligned=False), qual)
    
    for f in oFile.values():
        f.close()

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='bamCCSBarcodeDemultiplex.py', description='Export barcoded fastqs from unaligned bam and barcode h5 fofn')
    parser.add_argument('ccsBAM', metavar='ccsBAM', type=str,
                    help='BAM file containing CCS data')
    parser.add_argument('barcodeFofn', metavar='barcodeFofn', type=str,
                    help='fofn of barcode score files in H5 format')
    parser.add_argument('-a,--minPredictedAccuracy', dest='minPredictedAccuracy', type=float, default=0.9,
                    help='minimum read score. default 0.90' )
    parser.add_argument('-b,--minAveBarcodeScore', dest='minAvgBarcodeScore', default=0, type=int,
                    help='minimum average barcode score. default 0' )
    parser.add_argument('-n,--minNumPasses', dest='minNumPasses', default=3, type=int,
                    help='minimum number of passes. default 3' )
    parser.add_argument('-o,--outDir', dest='outDir', type=str, default=os.getcwd(),
                    help='output directory. default cwd' )
    parser.add_argument('-x', dest='extendedHeader', action='store_false', default=True,
                    help='Do not add accuracy information to fastq headers')

    main(parser)

 
