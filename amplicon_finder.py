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
import csv


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



def calculateChromosomeCoverageInBamFile(tumorCoverageData, chromosome, chromosomeSize, windowSize):
    start = 0
    cov   = [0 for i in range(0, (chromosomeSize // windowSize) + 1 )]
    coord = [100 * i for i in range(0, (chromosomeSize // windowSize) + 1 )]
    for i in range(0, chromosomeSize):
        cov[ i // windowSize] += tumorCoverageData[i]
    return (cov, coord)
    
def calculateChromosomeRelativeCoverage(tumorCoverageData, germlineCoverageData, chromosome, chromosomeSize, windowSize):
    relative_cov = []
    coord = []
    start = 0
    germ_cov = [0 for i in range(0, (chromosomeSize // windowSize) + 1 )]
    tumor_cov = [0 for i in range(0, (chromosomeSize // windowSize) + 1)]
    for i in range(0, chromosomeSize):
        germ_cov[ i // windowSize] += germlineCoverageData[i]
        tumor_cov[i // windowSize] += tumorCoverageData[i]
    for i in range(0, chromosomeSize // windowSize):
        relative_interval_cov = 1.0
        if germ_cov[i] != 0 and tumor_cov[i] != 0:
            relative_interval_cov = float(tumor_cov[i]) / float(germ_cov[i])
        start = i * windowSize
        relative_cov.append( (relative_interval_cov) )
        coord.append(start)
    return (relative_cov, coord)


def buildHmm(minAmpliconLength, maxGap, windowSize):
    b_bkgd_1 = 0.1
    a_interstate = b_bkgd_1 ** (2 * minAmpliconLength / windowSize)
    b_amp_0  = ( a_interstate ) ** (0.5 * windowSize / maxGap)
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

def predictAmpliconesInChromosome(cov, coord, coverageThreshold, minAmpliconLength, maxGap, windowSize):
    hmm = buildHmm(minAmpliconLength, maxGap, windowSize)
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

def LoadChromosomeCoverageFromBamFile(pathToBam, chromosome, chromosomeSize):
    samfile = pysam.AlignmentFile(pathToBam, "rb" )
    cov = []
    try:
        cov = [0 for i in range(0, chromosomeSize + 1 )]
        for pileupcolumn in samfile.pileup(chromosome, 0, chromosomeSize):
            cov[pileupcolumn.pos] += int(pileupcolumn.n)
    except ValueError as err:
        cov   = []
        print(err)
    samfile.close()
    return cov

def constructPathToCacheFile(pathToBam, chromosome, pathToCacheDir):
    cacheFileName = os.path.basename(pathToBam) + '.' +  chromosome + '.txt'
    pathToCacheFile = os.path.join( pathToCacheDir, cacheFileName)
    return pathToCacheFile

def LoadChromosomeCoverageFromCache(pathToCacheFile):
    cov = []
    if os.path.exists(pathToCacheFile):
        with open(pathToCacheFile, 'r') as cacheFile:
            cov = cacheFile.readlines()
            cov = [ int(x) for x in cov ] 
    return cov


def StoreChromosomeCoverageToCache(pathToCacheFile, coverageData):
    pathToCacheDir = os.path.dirname(pathToCacheFile)
    if not os.path.exists(pathToCacheDir):
        os.makedirs(pathToCacheDir)
    with open(pathToCacheFile, 'w') as cacheFile:
        for item in coverageData:
            cacheFile.write('{}\n'.format(item))

def LoadChromosomeCoverage(pathToBam, chromosome, chromosomeSize, pathToCacheDir):
    pathToCacheFile = constructPathToCacheFile(pathToBam, chromosome, pathToCacheDir)
    if not os.path.exists(pathToCacheFile):
        coverageData = LoadChromosomeCoverageFromBamFile(pathToBam, chromosome, chromosomeSize)
        StoreChromosomeCoverageToCache(pathToCacheFile, coverageData)
        return coverageData
    return LoadChromosomeCoverageFromCache(pathToCacheFile)
        

def getAmpliconsInChromosome(hg19ChromosomeSizesEntry):
    chromosome, chromosomeSize = hg19ChromosomeSizesEntry
    tumorCoverageData = LoadChromosomeCoverage(args.tumor_bam, chromosome, chromosomeSize, args.path_to_cache_dir)
    if args.germline_bam:
        germlineCoverageData = LoadChromosomeCoverage(args.germline_bam, chromosome, chromosomeSize, args.path_to_cache_dir)
        cov, coord = calculateChromosomeRelativeCoverage(tumorCoverageData, germlineCoverageData, chromosome, chromosomeSize, args.window)
    else:
        cov, coord = calculateChromosomeCoverageInBamFile(tumorCoverageData, chromosome, chromosomeSize, args.window)
    predictedAmplicons = predictAmpliconesInChromosome(cov, coord, args.cov_thsh, args.min_amplicon_length, args.max_gap, args.window)
    return (chromosome, predictedAmplicons)
    
def checkInputArguments(args):
    if not os.path.exists(args.tumor_bam):
        print("Specified tumor BAM file does not exist")
        return False

    if args.germline_bam and not os.path.exists(args.germline_bam):
        print("Specified germline BAM file does not exist")
        return False

    return True

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
    parser.add_argument('-c','--cache-dir', dest='path_to_cache_dir', type=str, default='./.amplicon_finder_cache',
                        help='Change path to cache dir (default is \'./.amplicon_finder\')')
    

    args = parser.parse_args()

    chromosomeList = [ch for ch in getChromosomeListFromBam(args.tumor_bam)]
    if not os.path.exists( args.path_to_cache_dir ):
        for chromosome, chromosomeSize in chromosomeList:
            LoadChromosomeCoverage(args.tumor_bam, chromosome, chromosomeSize, args.path_to_cache_dir)
            if args.germline_bam:
                LoadChromosomeCoverage(args.germline_bam, chromosome, chromosomeSize, args.path_to_cache_dir)


    partinionSize = (len(chromosomeList) + 1) / args.thr_num
    n = (len(chromosomeList) + 1) / args.thr_num
    pool = multiprocessing.pool.ThreadPool(processes=args.thr_num)
    for chromosome, predictedAmplicons in pool.map(getAmpliconsInChromosome, chromosomeList):
        printAmpliconInBedFormat(chromosome, predictedAmplicons)

