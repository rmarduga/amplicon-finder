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



def calculateChromosomeCoverageInBamFile(pathToBam, chromosome, chromosomeSize, window):
    samfile = pysam.AlignmentFile(pathToBam, "rb" )
    cov = []
    coord = []
    try:
        start = -1
        interval_cov = 0
        start = 0
        for pileupcolumn in samfile.pileup(chromosome, 0, chromosomeSize):
            if (pileupcolumn.pos // window) * window  != start:
                cov.append(interval_cov / window)
                coord.append(start)
                start_prev = start
                start = (pileupcolumn.pos // window) * window
                for i in range(start_prev + window, start, window):
                    cov.append(0)
                    coord.append(i)
                interval_cov = 0
            interval_cov += pileupcolumn.n
    except ValueError as err:
        print( err )
        cov   = []
        coord = []
    samfile.close()
    return (cov, coord)
    
def calculateChromosomeRelativeCoverage(pathToTumorBam, pathToGermlineBam, chromosome, chromosomeSize, window):
    germSamfile = pysam.AlignmentFile(pathToGermlineBam, "rb" )
    tumorSamfile = pysam.AlignmentFile(pathToTumorBam, "rb" )
    relative_cov = []
    coord = []
    try:
        start = -1
        germ_cov = [0 for i in range(0, (chromosomeSize // window) + 1 )]
        tumor_cov = [0 for i in range(0, (chromosomeSize // window) + 1)]
        start = 0
        for pileupcolumn in germSamfile.pileup(chromosome, 0, chromosomeSize):
            germ_cov[pileupcolumn.pos // window] += pileupcolumn.n
        for pileupcolumn in tumorSamfile.pileup(chromosome, 0, chromosomeSize):
            tumor_cov[pileupcolumn.pos // window] += pileupcolumn.n
        for i in range(0, chromosomeSize // window):
            relative_interval_cov = 1.0
            if germ_cov[i] != 0 and tumor_cov != 0:
                relative_interval_cov = float(tumor_cov[i]) / float(germ_cov[i])
            start = i * window
            relative_cov.append( (relative_interval_cov) )
            coord.append(start)
    except ValueError as err:
        cov   = []
        coord = []
        print err
    germSamfile.close()
    tumorSamfile.close()
    return (relative_cov, coord)


def buildHmm(minAmpliconLength, maxGap, window):
    b_bkgd_1 = 0.1
    a_interstate = b_bkgd_1 ** (2 * minAmpliconLength / window)
    b_amp_0  = ( a_interstate ) ** (0.5 * window / maxGap)
    b_amp_1  = 1 - b_amp_0
    b_bkgd_0 = 1 - b_bkgd_1
    bkgdDist = DiscreteDistribution({0: b_bkgd_0, 1: b_bkgd_1})
    ampDist  = DiscreteDistribution({0: b_amp_0,  1: b_amp_1})
    s_bkgd = State( bkgdDist, name='background' )
    s_amp = State( ampDist, name='amplicon' )
    hmm = HiddenMarkovModel()
    hmm.add_states(s_bkgd, s_amp)
    hmm.add_transition( hmm.start, s_bkgd, 1 - a_interstate )
    hmm.add_transition( hmm.start, s_amp,      a_interstate )
    hmm.add_transition( s_bkgd, s_bkgd,    1 - a_interstate)
    hmm.add_transition( s_bkgd, s_amp,         a_interstate)
    hmm.add_transition( s_amp,  s_bkgd,        a_interstate)
    hmm.add_transition( s_amp,  s_amp,     1 - a_interstate)
    hmm.bake()
    return hmm

def buildObservedSequence(cov, coord, coverageThreshold):
    covArr = np.array([cov]);
    covAvg = np.average(covArr)
    covStd = np.std(covArr)
    observedSeq = [1 if x > covAvg + coverageThreshold * covStd else 0 for x in cov]
    return observedSeq

def predictAmpliconesInChromosome(cov, coord, coverageThreshold, minAmpliconLength, maxGap, window):
    hmm = buildHmm(minAmpliconLength, maxGap, window)
    observedSeq = buildObservedSequence(cov, coord, coverageThreshold)
    hmm_predictions = hmm.predict( observedSeq, algorithm="viterbi" )
    start = 0
    predictedAmplicons = []
    for i in range(1, len(observedSeq)):
        if hmm_predictions[i] == 0 and ( hmm_predictions[i-1] == 1 or hmm_predictions[i-1] == 2 ):
            start = i
        elif ( hmm_predictions[i] == 1 or hmm_predictions[i] == 2 ) and hmm_predictions[i-1] == 0:
            predictedAmplicons.append((coord[start], coord[i]))
    return predictedAmplicons

def printAmpliconInBedFormat(ch, predictedAmplicons):
    if len(predictedAmplicons) > 0:
        print('\n'.join("{}\t{}\t{}".format(ch,  predictedAmplicons[i][0], predictedAmplicons[i][1]) for i in range(0,len(predictedAmplicons))))


def getAmpliconsInChromosome(hg19ChromosomeSizesEntry):
    chromosome, size = hg19ChromosomeSizesEntry
    if args.germline_bam:
        cov, coord = calculateChromosomeRelativeCoverage(args.tumor_bam, args.germline_bam, chromosome, size, args.window)
    else:
        cov, coord = calculateChromosomeCoverageInBamFile(args.tumor_bam, chromosome, size, args.window)
    predictedAmplicons = predictAmpliconesInChromosome(cov, coord, args.cov_thsh, args.min_amplicon_length, args.max_gap, args.window)
    return (chromosome, predictedAmplicons)
    
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
    parser.add_argument('-g','--max-gap', dest='max_gap', type=int, default=1000,
                        help='Maximum gap in basepairs between amplicons \
                        that are considered as a single amplicon. (default is 1000)')
    parser.add_argument('-l','--min-amplicon-length', dest='min_amplicon_length', type=int, default=1000,
                        help='Minimum amplicon length in basepairs to be captured. (default is 1000)')

    args = parser.parse_args()
    chromosomeList = [ch for ch in getChromosomeListFromBam(args.tumor_bam)]
    partinionSize = (len(chromosomeList) + 1) / args.thr_num
    n = (len(chromosomeList) + 1) / args.thr_num
    pool = multiprocessing.pool.ThreadPool(processes=args.thr_num)
    for chromosome, predictedAmplicons in pool.map(getAmpliconsInChromosome, chromosomeList):
        printAmpliconInBedFormat(chromosome, predictedAmplicons)

