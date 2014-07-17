#!/usr/bin/env python

from __future__ import division, print_function

import sys
import os.path
import argparse
import itertools
import json
import csv
import warnings

import numpy as np
import h5py

from operator import itemgetter
from pbcore.io import CmpH5Reader, FastaReader

baseDict = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "AC": "M", "AG": "R", "AT": "W",
    "CG": "S", "CT": "Y", "GT": "K",
    "ACG": "V", "ACT": "H", "AGT": "D", "CGT": "B",
    "ACGT": "N", "": ""}

gapBase = "-"

def rle(string):
    return [(char, sum(1 for _ in group)) for char, group in itertools.groupby(string)]

def rleTwo(string):
    res = []
    for c, n in rle(string):
        if n > 1:
            res.append((c, n - 1))
        res.append((c, 1))
    return res

def rleNot(string):
    return [(c, 1) for c in string]

def rld(encoding):
    return "".join(c * n for c, n in encoding)

class Aligner(object):
    bases = ["A", "C", "G", "T"]
    # matchCost, errorCost, ambigCost = 4, -1, 2
    matchCost, errorCost, ambigCost = 3, -1, 2
    deletion, insertion, match, mismatch = xrange(4)

    def __init__(self, reference, query, whichRow=-1):
        self.refRle = rleNot(reference)
        self.qryRle = rleNot(query)

        self.ref = "".join(b for b, _ in self.refRle)
        self.qry = "".join(b for b, _ in self.qryRle)

        # map whichRow into new RLE space
        if whichRow >= 0:
            rowIdx = 0
            for i, (b, c) in enumerate(self.refRle):
                rowIdx += c
                if rowIdx > whichRow:
                    whichRow = i
                    break

        self.whichRow = whichRow
        self.tbl = [[[] for col in xrange(len(self.qry))] for row in xrange(len(self.ref))]
        self.buildTable()

    def bestMoves(self, refIdx, qryIdx, matchCost=None, refBase=None):
        if matchCost is None:
            matchCost = Aligner.matchCost

        if refBase is None:
            refBase = self.ref[refIdx]

        allMoves = []

        deletionScore  = max(self.tbl[refIdx - 1][qryIdx]) if refIdx > 0 else 0
        allMoves.append((Aligner.deletion, deletionScore + Aligner.errorCost))

        insertionScore = max(self.tbl[refIdx][qryIdx - 1]) if qryIdx > 0 else 0
        allMoves.append((Aligner.insertion, insertionScore + Aligner.errorCost))

        matchScore = max(self.tbl[refIdx - 1][qryIdx - 1]) if refIdx > 0 and qryIdx > 0 else 0

        if (len(self.tbl[refIdx][qryIdx]) == len(Aligner.bases)):
            allMoves.append((Aligner.match, matchScore + Aligner.ambigCost))
        elif (refBase == self.qry[qryIdx]):
            allMoves.append((Aligner.match, matchScore + matchCost))
        else:
            # allMoves.append((Aligner.mismatch, matchScore + 2 * Aligner.errorCost))
            allMoves.append((Aligner.mismatch, matchScore + Aligner.errorCost))

        bestScore = max(s for _, s in allMoves)

        return [(m, s) for m, s in allMoves if s >= bestScore]

    def updateTable(self, refIdx, qryIdx, matchCost, refBase=None):
        bestScore = self.bestMoves(refIdx, qryIdx, matchCost, refBase)[0][1]
        self.tbl[refIdx][qryIdx].append(bestScore)

    def buildTable(self):
        for qryIdx in xrange(len(self.qry)):
            for refIdx in xrange(len(self.ref)):
                if refIdx == self.whichRow:
                    for refBase in Aligner.bases:
                        self.updateTable(refIdx, qryIdx, Aligner.ambigCost, refBase=refBase)
                else:
                    self.updateTable(refIdx, qryIdx, Aligner.matchCost)

    def bestEntry(self):
        lRefIdx = len(self.ref) - 1
        lQryIdx = len(self.qry) - 1

        bestRefIdx = max(xrange(len(self.ref)), key=lambda i: max(self.tbl[i][lQryIdx]))
        bestQryIdx = max(xrange(len(self.qry)), key=lambda i: max(self.tbl[lRefIdx][i]))

        if max(self.tbl[bestRefIdx][lQryIdx]) > max(self.tbl[lRefIdx][bestQryIdx]):
            return (bestRefIdx, lQryIdx)

        return (lRefIdx, bestQryIdx)

    def countRow(self, minAccuracy=0.9):
        if not self.ref or not self.qry:
            return []

        refIdx, qryIdx = self.bestEntry()
        obsAcc = []
        nMatch, nError = 0, 0

        obsList = []
        moveStates = []

        while True:
            newMoves = self.bestMoves(refIdx, qryIdx)

            for moveScore in newMoves:
                newMove   = moveScore[0]
                newRefIdx = refIdx - (1 if newMove != Aligner.insertion else 0)
                newQryIdx = qryIdx - (1 if newMove != Aligner.deletion else 0)
                newObsAcc = obsAcc
                newNMatch = nMatch + (1 if newMove == Aligner.match else 0)
                newNError = nError + (0 if newMove == Aligner.match else 1)

                if newMove == Aligner.match and refIdx == self.whichRow:
                    bestScore = max(self.tbl[refIdx][qryIdx])
                    newObsAcc = [v for v in itertools.chain(obsAcc, [
                        Aligner.bases[i]
                        for i in xrange(len(Aligner.bases))
                        if self.tbl[refIdx][qryIdx][i] >= bestScore])]

                if newRefIdx < 0 or newQryIdx < 0:
                    newNError += max(newRefIdx, newQryIdx) + 1
                    obsList.extend([(obs, newNMatch, newNError) for obs in newObsAcc])
                else:
                    moveStates.append((newRefIdx, newQryIdx, newObsAcc, newNMatch, newNError))

            if not moveStates:
                break

            refIdx, qryIdx, obsAcc, nMatch, nError = moveStates.pop()

        if not obsList:
            return []

        minAccuracy = max([minAccuracy] + [nM/(nM + nE) for _, nM, nE in obsList])
        bases = uniq(o for o, nM, nE in obsList if nM/(nM + nE) >= minAccuracy)

        return bases

    def localAlign(self):
        if not self.ref and not self.qry:
            return "", ""

        if not self.ref:
            qry = rld(self.qryRle)
            return gapBase * len(qry), qry

        if not self.qry:
            ref = rld(self.refRle)
            return ref, gapBase * len(ref)

        refIdx, qryIdx = self.bestEntry()

        refAligned = []
        qryAligned = []

        for i in xrange(len(self.ref) - 1, refIdx, -1):
            nBases = self.refRle[i][1]
            refAligned.append(self.ref[i] * nBases)
            qryAligned.append(gapBase * nBases)

        for i in xrange(len(self.qry) - 1, qryIdx, -1):
            nBases = self.qryRle[i][1]
            refAligned.append(gapBase * nBases)
            qryAligned.append(self.qry[i] * nBases)

        while (refIdx >= 0 and qryIdx >= 0):
            move = max(self.bestMoves(refIdx, qryIdx), key=itemgetter(1))[0]

            if move == Aligner.deletion:
                nBases = self.refRle[refIdx][1]
                refAligned.append(self.ref[refIdx] * nBases)
                qryAligned.append(gapBase * nBases)
                refIdx -= 1
            elif move == Aligner.insertion:
                nBases = self.qryRle[qryIdx][1]
                refAligned.append(gapBase * nBases)
                qryAligned.append(self.qry[qryIdx] * nBases)
                qryIdx -= 1
            else:
                nRefBases = self.refRle[refIdx][1]
                nQryBases = self.qryRle[qryIdx][1]
                nBases = max(nRefBases, nQryBases)
                refAligned.append(self.ref[refIdx] * nRefBases + (gapBase * (nBases - nRefBases)))
                qryAligned.append(self.qry[qryIdx] * nQryBases + (gapBase * (nBases - nQryBases)))
                refIdx -= 1
                qryIdx -= 1

        while refIdx >= 0:
            nBases = self.refRle[refIdx][1]
            refAligned.append(self.ref[refIdx] * nBases)
            qryAligned.append(gapBase * nBases)
            refIdx -= 1

        while qryIdx >= 0:
            nBases = self.qryRle[qryIdx][1]
            refAligned.append(gapBase * nBases)
            qryAligned.append(self.qry[qryIdx] * nBases)
            qryIdx -= 1

        return "".join(refAligned[::-1]), "".join(qryAligned[::-1])

    @staticmethod
    def test1():
        ref = "ACACACTA"
        qry = "AGCACACA"

        refAln, qryAln = Aligner(ref, qry).localAlign()

        print(refAln)
        print(qryAln)

        assert refAln == "A-CACACTA"
        assert qryAln == "AGCACAC-A"

    @staticmethod
    def test2():
        ref = "CACACTA"
        qry = "AGCACAC"

        refAln, qryAln = Aligner(ref, qry).localAlign()

        print(refAln)
        print(qryAln)

        assert refAln == "--CACACTA"
        assert qryAln == "AGCACAC--"

    @staticmethod
    def test3():
        ref = "ACACACA"
        qry = "ACATACA"

        refAln, qryAln = Aligner(ref, qry).localAlign()

        print(refAln)
        print(qryAln)

        assert refAln == "ACACACA"
        assert qryAln == "ACATACA"

    @staticmethod
    def test4():
        ref = "ACGCACA"
        qry = "ACGACA"

        refAln, qryAln = Aligner(ref, qry).localAlign()

        print(refAln)
        print(qryAln)

        assert refAln == "ACGCACA"
        assert qryAln == "ACG-ACA"

    @staticmethod
    def test5():
        ref = "ACACACA"
        qry = "ACATACA"

        aln = Aligner(ref, qry, 3)
        refAln, qryAln = aln.localAlign()
        bases = aln.countRow()

        print(refAln)
        print(qryAln)
        print(bases)

        assert refAln == "ACACACA"
        assert qryAln == "ACATACA"
        assert bases == ["T"]

    @staticmethod
    def test6():
        ref = "CAGGAATGGAGTGGCCCAAAAG"
        qry = "CAGGAATGGATGGCCCAAAAG"

        aln = Aligner(ref, qry, 10)
        refAln, qryAln = aln.localAlign()
        bases = aln.countRow()

        print(refAln)
        print(qryAln)
        print(bases)

        assert not bases

    @staticmethod
    def test7():
        ref = "GGAGTGG"
        qry = "GGATGG"

        aln = Aligner(ref, qry, 3)
        refAln, qryAln = aln.localAlign()
        bases = aln.countRow()

        print(refAln)
        print(qryAln)
        print(bases)

        assert not bases

    @staticmethod
    def test8():
        ref = "CCGTTTCTCCTGGCTCAGTTT"
        qry = "CCGTTTCTCATGGCTTCAGTTT"

        aln = Aligner(ref, qry, 10)
        refAln, qryAln = aln.localAlign()
        bases = aln.countRow()

        print(refAln)
        print(qryAln)
        print(bases)

        assert bases and bases == ["T"]

    @staticmethod
    def test9():
        ref = "CCGTTTCTCCTGGCTCAGTTT"
        qry = ""

        aln = Aligner(ref, qry, 10)
        refAln, qryAln = aln.localAlign()
        bases = aln.countRow()

        print(refAln)
        print(qryAln)
        print(bases)

        assert refAln == ref
        assert qryAln == ("-" * len(ref))
        assert not bases

    @staticmethod
    def testAll():
        Aligner.test1()
        Aligner.test2()
        Aligner.test3()
        Aligner.test4()
        Aligner.test5()
        Aligner.test6()
        Aligner.test7()
        Aligner.test8()
        Aligner.test9()
        print("PASSED")

def parseArguments(args):
    if "--testAligner" in args:
        Aligner.testAll()
        sys.exit(0)

    parser = argparse.ArgumentParser(description="Call which reads support which variants")
    parser.add_argument("--minAccuracy", type=float, default=0.9, help="local realignment minimum accuracy (default=0.9)")
    parser.add_argument("--radius", type=int, default=10, help="local realignment radius (default=10)")
    parser.add_argument("--cmp", required=True, help="input cmp.h5 file")
    parser.add_argument("--variants", required=True, help="input minor_variants.csv file")
    parser.add_argument("--reference", required=True, help="input reference FASTA file")
    parser.add_argument("--output", help="output directory (default=.)", default=".")

    return parser.parse_args(args)

def loadCmpH5(filename):
    filename = os.path.abspath(os.path.expanduser(filename))
    handle = h5py.File(filename, "r")
    return CmpH5Reader(handle)

def uniq(iterable):
    return sorted(set(iterable))

def contigVariantPositions(variantsCsv):
    # pull out the contig, position, and reference base
    contigRefPositions = zip(variantsCsv["CONTIG"], variantsCsv["POSITION"], variantsCsv["REF"], variantsCsv["ALT"])
    # group by contig
    byContig = itertools.groupby(contigRefPositions, key=itemgetter(0))
    # build a dictionionary of {contig: [(position, refBase, altBases)]} where position is 0-indexed
    contigVariants = dict(
        (contig, sorted(set(
            (p - 1, "".join(uniq(cpra[2] for cpra in posData)), "".join(uniq(cpra[3] for cpra in posData)))
            for p, posData in (
                (a, list(b))
                for a, b in itertools.groupby(contigData, key=itemgetter(1)))), key=itemgetter(0)))
        for contig, contigData in byContig)
    return contigVariants

def ccsName(read):
    return "{0:s}/{1:d}/ccs".format(read.movieInfo.Name, read.HoleNumber)

def dispAlign(ref, qry):
    alnChars = []

    for r, q in zip(ref, qry):
        if r == gapBase or q == gapBase:
            alnChars.append(" ")
        elif r == q:
            alnChars.append("|")
        else:
            alnChars.append("*")

    # qryChars = [q.upper() if a == "|" else q.lower() for a, q in zip(alnChars, qry)]

    return "\n".join((ref, "".join(alnChars), qry))


def callReadVariants(reference, read, variantSites, minAccuracy, radius, nRead):
    readName = ccsName(read) if read.cmpH5.readType == "CCS" else read.readName
    resultRow = [readName]

    # print("processing read: {0:d}\r".format(nRead), end="", file=sys.stderr)

    readString = read.read(orientation="genomic")
    refPositions = read.referencePositions(orientation="genomic")
    alnLen = len(readString)

    def baseInWindow(start, end):
        return lambda i: (
            refPositions[i] >= windowStart and
            refPositions[i] < windowEnd and
            readString[i] != gapBase)

    assert alnLen == len(refPositions)

    for site, refBase, altBases in variantSites:
        windowStart = max(site - radius, 0)
        windowEnd = min(site + radius + 1, len(reference))

        if not read.spansReferenceRange(windowStart, windowEnd):
            resultRow.append("")
            continue

        isValidBase = baseInWindow(windowStart, windowEnd)

        refWindow = reference[windowStart:windowEnd]
        readWindow = "".join(readString[i] for i in xrange(alnLen) if isValidBase(i))

        aln = Aligner(refWindow, readWindow, site - windowStart)
        bases = "".join(aln.countRow(minAccuracy))
        ambigBase = baseDict.get(bases, "X")

        if len(bases) == 1 and bases != refBase and bases in altBases:
            print(">{0:s}/{1:s}{2:d}{3:s}".format(readName, refBase, site + 1, bases))
            print(dispAlign(*aln.localAlign()), end="\n\n")

        resultRow.append(ambigBase)

    return resultRow

def callVariants(cmpH5, contigVariants, referenceFasta, minAccuracy, radius):
    readVariants = {}

    for contigName, variantSites in contigVariants.iteritems():
        contigName = contigName.strip("\"")
        try:
            refSeq = None
            with FastaReader(referenceFasta) as fr:
                for contig in fr:
                    if contig.name == contigName:
                        refSeq = contig.sequence
                        break
            if refSeq is None:
                warnings.warn("Cannot find contig: \"{0:s}\" in the FASTA: \"{1:s}".format(contigName, referenceFasta))
                continue
            refInfo = cmpH5.referenceInfo(contigName)
            refId = refInfo.ID
            refLen = refInfo.Length
            readVariants[contigName] = [["READ"] + ["{0:s}{1:d}{2:s}".format(ref, pos + 1, baseDict.get(alt, "X")) for pos, ref, alt in variantSites]]
            readVariants[contigName].extend(
                callReadVariants(refSeq, cmpH5[readRow], variantSites, minAccuracy, radius, nRead)
                for nRead, readRow
                in enumerate(cmpH5.readsInRange(refId, 0, refLen, justIndices=True), start=1))
        except KeyError:
            warnings.warn("Cannot find contig: \"{0:s}\" in the cmpH5: \"{1:s}\"".format(contig, cmpH5.filename))

    return readVariants

def main(args=sys.argv[1:]):
    options = parseArguments(args)
    cmpH5 = loadCmpH5(options.cmp)
    variantsCsv = np.recfromtxt(options.variants, delimiter=",", names=True)
    contigVariants = contigVariantPositions(variantsCsv)

    with open(os.path.join(options.output, "read_variants.json"), "w") as handle:
        result = callVariants(cmpH5, contigVariants, options.reference, options.minAccuracy, options.radius)
        json.dump(result, handle, indent=4, separators=(",", ": "))

    input = open("read_variants.json")
    data = json.load(input)
    input.close()

    with open('read_variants.csv', 'w') as csvfile:
        output = csv.writer(csvfile, delimiter=',')

        for row in data:
            for i in range(0,len(data[row])):
                a = data[row][i]
                a.insert(0, row)
                output.writerow(a)

    return 0

if __name__ == "__main__":
    sys.exit(main())
