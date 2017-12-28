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
import argparse

#def setAffinity():
#    proc_num_f = os.popen( 'cat /proc/cpuinfo | grep ''^processor'' | wc -l')
#    proc_num = int(proc_num_f.read())
#    proc_num_f.close()
#    affinity_flag = pow(2, proc_num) - 1
#    #print (hex(affinity_flag))
#    os.system(('taskset -p {} {}'.format(hex(affinity_flag), os.getpid())))
#    affinity_f = os.popen('taskset -p {} {}'.format(hex(affinity_flag), os.getpid()))
#    affinity_f.close()
#
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
    
def calculateChromosomeRelativeCoverage(pathToTumorBam, pathToGermlineBam, chromosome, size, blockSize):
    block_size = 100
    germSamfile = pysam.AlignmentFile(pathToGermlineBam, "rb" )
    tumorSamfile = pysam.AlignmentFile(pathToTumorBam, "rb" )
    relative_cov = []
    coord = []
    try:
        start = -1
        germ_cov = [0 for i in range(0, (size // block_size) + 1 )]
        tumor_cov = [0 for i in range(0, (size // block_size) + 1)]
        start = 0
        for pileupcolumn in germSamfile.pileup(chromosome, 0, size):
            #if pileupcolumn.pos <= size:
            germ_cov[pileupcolumn.pos // block_size] += pileupcolumn.n
        for pileupcolumn in tumorSamfile.pileup(chromosome, 0, size):
            #if pileupcolumn.pos <= size:
            tumor_cov[pileupcolumn.pos // block_size] += pileupcolumn.n
        for i in range(0, size // block_size):
            relative_interval_cov = 1.0
            if germ_cov[i] != 0 and tumor_cov != 0:
                relative_interval_cov = float(tumor_cov[i]) / float(germ_cov[i])
            start = i * block_size
            relative_cov.append( (relative_interval_cov) )
            coord.append(start)
    except ValueError as err:
        cov   = []
        coord = []
        print err
    germSamfile.close()
    tumorSamfile.close()
    return (relative_cov, coord)


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
    if args.germline_bam:
        cov, coord = calculateChromosomeRelativeCoverage(args.tumor_bam, args.germline_bam, chromosome, size, args.window)
    else:
        cov, coord = calculateChromosomeCoverageInBamFile(args.tumor_bam, chromosome, size, args.window)
    amps = predictAmpliconesInChromosome(cov, coord, args.cov_thsh)
    return (chromosome, amps)
    
if __name__ == '__main__':
#    setAffinity()
    parser = argparse.ArgumentParser(description='Find the amplicons in a BAM file.')
    parser.add_argument('tumor_bam', metavar='<tumor_bam>', type=str, help='Path to the tumor bam file')
    parser.add_argument('germline_bam', metavar='<germline_bam>', type=str, nargs='?', help='Path to the germline bam file. If provided, relative coveage between tumor and germline samples will be used. (Optional)')
    parser.add_argument('-t','--threads', dest='thr_num', type=int, default=2, help='number of threads (default is 1)')
    parser.add_argument('-k','--coverage-threshold', dest='cov_thsh', type=int, default=4, help='Coverage threshold in standard diviation units. (default is 3)')
    parser.add_argument('-i','--interval', dest='window', type=int, default=100,
                        help='Length of the interval of DNA that is being averaged. \
                        It could be seen as the resolution of the method. (default is 100)')

    args = parser.parse_args()
    chromosomeList = [ch for ch in getChromosomeListFromBam(args.tumor_bam)]
    partinionSize = (len(chromosomeList) + 1) / args.thr_num
    n = (len(chromosomeList) + 1) / args.thr_num
    pool = multiprocessing.pool.ThreadPool(processes=args.thr_num)
    for ch, amps in pool.map(getAmpliconsInChromosome, chromosomeList):
        printAmpliconInBedFormat(ch, amps)

