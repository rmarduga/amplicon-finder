#!/usr/bin/python
"""
Created on Fri Feb 24 18:06:55 2017

@author: rmarduga
"""

import pysam
import os
import numpy as np
from pomegranate import DiscreteDistribution, State, HiddenMarkovModel
import multiprocessing
import sys
import argparse
os.system('taskset -p 0xffffffff %d' % os.getpid())


def getChromosomeListFromBam(path):
    rows = pysam.view("-H", path)
    r = []
    xrows = []
    for c in rows:
        if c == '\n':
            xrows.append(''.join(r))
            r = []
        else:
            r.append(c)
    for x in xrows:
        cols = x.split()
        if cols[0][:3] != '@SQ':
            continue
        yield((cols[1][3:], int(cols[2][3:])))



def calculateChromosomeCoverageInBamFile(pathToBam, chromosome, size, blockSize):
    samfile = pysam.AlignmentFile(pathToBam, "rb" )
    cov = []
    coord = []
    try:
        start = -1
        interval_cov = 0
        start = 0
        for pileupcolumn in samfile.pileup(chromosome, 0, size):
            if (pileupcolumn.pos // blockSize) * blockSize  != start:
                cov.append(interval_cov / blockSize)
                coord.append(start)
                start = (pileupcolumn.pos // blockSize) * blockSize
                interval_cov = 0
            interval_cov += pileupcolumn.n
    except ValueError as err:
        print( err )
        cov   = []
        coord = []
    samfile.close()
    return (cov, coord)
    


def predictAmpliconesInChromosome(cov, coord, coverageThreshold):
    covArr = np.array([cov]);
    covAvg = np.average(covArr)
    covStd = np.std(covArr)

    amp = [1 if x > covAvg + coverageThreshold * covStd else 0 for x in cov]

    bkgdDist = DiscreteDistribution({0: 0.85, 1: 0.15})
    ampDist  = DiscreteDistribution({0: 0.15, 1: 0.85})

    s1 = State( bkgdDist, name='background' )
    s2 = State( ampDist, name='amplicon' )

    alpha = 1e-20
    beta = 1e-20

    hmm = HiddenMarkovModel()
    hmm.add_states(s1, s2)
    hmm.add_transition( hmm.start, s1, 1-alpha )
    hmm.add_transition( hmm.start, s2, alpha )
    hmm.add_transition( s1, s1, 1 - alpha)
    hmm.add_transition( s1, s2, alpha)
    hmm.add_transition( s2, s1, beta)
    hmm.add_transition( s2, s2, 1 - beta)

    hmm.bake()
    hmm.freeze_distributions()

    hmm_predictions = hmm.predict( amp, algorithm="viterbi" )
    start = 0
    amps = []
    for i in range(1, len(amp)):
        if hmm_predictions[i] == 0 and ( hmm_predictions[i-1] == 1 or hmm_predictions[i-1] == 2 ):
            start = i
        elif ( hmm_predictions[i] == 1 or hmm_predictions[i] == 2 ) and hmm_predictions[i-1] == 0:
            amps.append((coord[start], coord[i]))
    return amps

def printAmpliconInBedFormat(ch, amps):
    print('\n'.join("{}\t{}\t{}".format(ch,  amps[i][0], amps[i][1]) for i in range(0,len(amps))))


def getAmpliconsInChromosome(hg19ChromosomeSizesEntry):
    chromosome, size = hg19ChromosomeSizesEntry
    cov, coord = calculateChromosomeCoverageInBamFile(args.bam, chromosome, size, args.window)
    amps = predictAmpliconesInChromosome(cov, coord, args.cov_thsh)
    return (chromosome, amps)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find the amplicons in a BAM file.')
    parser.add_argument('bam', metavar='<bam>', type=str, help='Path to the tumor bam file')
    parser.add_argument('-t','--threads', dest='thr_num', type=int, default=2, help='number of threads (default is 1)')
    parser.add_argument('-k','--coverage-threshold', dest='cov_thsh', type=int, default=4, help='Coverage threshold in standard diviation units. (default is 3)')
    parser.add_argument('-i','--interval', dest='window', type=int, default=100,
                        help='Length of the interval of DNA that is being averaged. \
                        It could be seen as the resolution of the method. (default is 100)')
    args = parser.parse_args()
    chromosomeList = [ch for ch in getChromosomeListFromBam(args.bam)]
    partinionSize = (len(chromosomeList) + 1) / args.thr_num
    n = (len(chromosomeList) + 1) / args.thr_num
    pool = multiprocessing.pool.ThreadPool(processes=args.thr_num)
    for ch, amps in pool.map(getAmpliconsInChromosome, chromosomeList):
        printAmpliconInBedFormat(ch, amps)

