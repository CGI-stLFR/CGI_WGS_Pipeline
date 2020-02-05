#!/usr/bin/env python
#coding: utf8

###### Import Modules
import sys, os

# import json
# import pprint
# import cPickle as pickle

###### Document Decription
'''  '''

###### Version and Date
prog_version = '0.1.0'
prog_date = '2015-09-24'

###### Usage
usage = '''

     Version %s  by Vincent Li  %s

     Usage: %s <gcBlock> >STDOUT
''' % (prog_version, prog_date, os.path.basename(sys.argv[0]))

######## Global Variable


#######################################################################
############################  BEGIN Class  ############################
#######################################################################
class gcBin(object):
    """ Read GC Block file, calculate normalized GC for ploting """
    def __init__(self, gcFile):
        self.gcFile = gcFile
        self.data = []
        self.plotArr = []
        self.winSize = 0
        self.sumCov = 0
        self.sumCount = 0

        self.readFile(gcFile)
        self.normalizeCov()
    
    def readFile(self, f):
        ########### Open and read input file
        try:
            fhIn = open(f, 'r')
        except Exception as e:
            raise e
        
        for line in fhIn:
            info = line.split()
            intArr = map(int, info)
            self.data.append(intArr)
            self.winSize = intArr[0]
            self.sumCount += intArr[1]
            self.sumCov += intArr[2]
        fhIn.close
        return self.data

    def calGcCov(self):
        thresh = 0.6 * self.avgSumCov
        basesBelowThreshold = 0
        for v in self.data:
            if v[1] == 0:
                continue
            if float(v[2]) / v[1] < thresh:
                # print [v[0], float(v[2]) / v[1], thresh]
                basesBelowThreshold += v[1]
        return 100.0 * basesBelowThreshold / self.sumCount

    def normalizeCov(self):
        currSum = 0.0
        self.avgSumCov = 1.0 * self.sumCov / self.sumCount
        for e in self.data:
            currSum += e[1]
            xAxis = 100 * currSum/self.sumCount
            if e[1] == 0:
                yAxis = 0
            else:
                yAxis = (e[2] * 1.0 / e[1]) / self.avgSumCov ## (gc_coverage[k]/gc_count[k]) / (sum_gc_coverage/sum_gc_count)
            self.plotArr.append([xAxis, yAxis])
        return self.plotArr

    def printPlot(self):
        for e in self.plotArr:
            print "%.4f\t%.4f" % tuple(e)

##########################################################################
############################  BEGIN Function  ############################
##########################################################################


######################################################################
############################  BEGIN Main  ############################
######################################################################
#################################
##
##   Main function of program.
##
#################################
def main():
    
    ######################### Phrase parameters #########################
    import argparse
    ArgParser = argparse.ArgumentParser(usage = usage, version = prog_version)

    (para, args) = ArgParser.parse_known_args()

    if len(args) != 1:
        ArgParser.print_help()
        print >>sys.stderr, "\nERROR: The parameters number is not correct!"
        sys.exit(1)
    else:
        (gcBlockFile,) = args

    ############################# Main Body #############################
    gcObj = gcBin(gcBlockFile)
    # print gcObj.calGcCov()
    gcObj.printPlot()



#################################
##
##   Start the main program.
##
#################################
if __name__ == '__main__':
    main()

################## God's in his heaven, All's right with the world. ##################