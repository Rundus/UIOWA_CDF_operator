#!/usr/bin/env python

__author__ = "Scott Bounds"
__date__ = "2015-09-29"
__version__ = "1.0.1"

import sys
import os
import os.path
import struct as S
import datetime as dt
import numpy
#os.putenv('CDF_LIB', '/usr/local/cdf/lib')
from spacepy import pycdf
import glob
#import time
project="CAPERII_52005"
#g_sSync = "\xfa\xf3\x34\x00"  #CAPER The real frame sync but the data is little endian so...
g_sSync = "\x34\xf3\xfa"  #CAPER  removed the \x00 because compacted data had 2  supportingDocs bits of information  compressed into the sync word
#g_sSync = "\x02\x04\xf3\xfa"  #CAPER2?
#  "\xfe\x6b\x28\x40"    #TRICE
#g_sSync = "\x40\x28\x6b\xfe"  #TRICE
version='v01'

basedate=dt.datetime(2000, 1, 1)


##############################################################################
def prnHelp(fOut):
    fOut.write("""tad_file_convert Usage:
	tad_file_convert DIR
    DIR: Directory of the data to convert

Description:
	Prints binary information from CAPER rocket telemetry files\n""")


##############################################################################
def timestamp(year, nFrmHdr):
    'Converts frame header into a datetime timestamp (year, frmhdr, type)'
    'where type:"Text" returns a text string, else it returns a datetime tuple'
    doy = (nFrmHdr[0] >> 24 & 15) * 100 + (nFrmHdr[0] >> 20 & 15) * 10 + (nFrmHdr[0] >> 16 & 15)
    hours = (nFrmHdr[0] >> 12 & 15) * 10 + (nFrmHdr[0] >> 8 & 15)
    mins = (nFrmHdr[0] >> 4 & 15) * 10 + (nFrmHdr[0] & 15)
    secs = (nFrmHdr[1] >> 28 & 15) * 10 + (nFrmHdr[1] >> 24 & 15)
    a = nFrmHdr[1] >> 20 & 15
    b = nFrmHdr[1] >> 16 & 15
    c = nFrmHdr[1] >> 12 & 15
    d = nFrmHdr[1] >> 8 & 15
    e = nFrmHdr[1] >> 4 & 15
    f = nFrmHdr[1] & 15

    #if a > 9:
	#print "Value " + str(a) + " in 100s msec"
       # a = 7
    #if b > 9:
	#print "Value " + str(b) + " in 10s msec"
       # b = 7
    #if c > 9:
        #print "Value " + str(c) + " in 1s msec"
        #c = 7

    #if d > 9:
	#print "Value " + str(d) + " in 100s usec"
       #d = 7
    #if e > 9:
	#print "Value " + str(e) + " in 10s usec"
       #e = 7
    if f > 9:
	#print "Value " + str(f) + " in 1s usec"
        f = 7

    msec = a * 100 + b * 10 + c
    usec = d * 100 + e * 10 + f

#    if nFrmHdr[1] >> 20 & 15 > 9 or nFrmHdr[1] >> 16 & 15 >9 or nFrmHdr[1] >> 12 & 15 > 9:
#	print msec, usec

    if hours > 23:
        doy += 1
        hours = 0
    try:
	d = dt.datetime.strptime(year + ' ' + str(doy), '%Y %j')
    	t = dt.time(hours, mins, secs, msec * 1000 + usec, tzinfo=None)
    	tt = dt.datetime.combine(d, t)

    except ValueError:
	tt = dt.datetime(2000, 1, 1)
	print year, doy, hours, mins, secs, msec, usec



    return tt
    #return pycdf.lib.datetime_to_tt2000(tt)


###############################################################################
def unpackData(sRec, nNomDataLen):
    'Unpacks the reacord into integer values for all data in a minor frame'

    nData = numpy.array(S.unpack("59I", sRec[16:]))
    nSync = long(numpy.array(S.unpack("I", sRec[12:16])))

    # Unpack XXX data numbers

    newData = [0]*nNomDataLen
    newData[0::2] = nData>>16
    newData[1::2] = nData & 65535
    newData[nNomDataLen-1:] = [int(nSync>>22), int((nSync>>12) & 1023), int((nSync>>2 ) & 1023)]

    newData = [int(i) for i in newData]	
    return newData

def unpack10bitData(sRec, nNomDataLen):
    'Unpacks the record into integer values for all data in a major frame'
    '16 bit words into 10 bit words'

    nData = numpy.array(S.unpack("38I", sRec[12:]))
    nSync = long(numpy.array(S.unpack("I", sRec[12:16])))
    #print ":".join("{:02x}".format(c) for c in nData)

    newData = [0]*76
    newData[0::2] = nData>>16
    newData[1::2] = nData &  65535
    newData = numpy.concatenate( (newData, [0,0,0,0]), axis=None)
    #print len(newData)
    # Unpack XXX data numbers
    #print ":".join("{:04x}".format(c) for c in newData[0:5])
    #print ":".join("{:02x}".format(nSync))
    newnewData = [0]*128
    #print ":".join("{:02x}".format(c) for c in numpy.right_shift(newData[0::5],6))
    newnewData[0::8] = numpy.right_shift(newData[0::5],6)
    newnewData[1::8] = numpy.left_shift(numpy.bitwise_and(newData[0::5] , 63),4) + numpy.right_shift(newData[1::5], 12)
    newnewData[2::8] = numpy.bitwise_and(numpy.right_shift(newData[1::5], 2), 1023)
    newnewData[3::8] = numpy.left_shift(numpy.bitwise_and(newData[1::5], 3),  8) + numpy.right_shift(newData[2::5], 8)
    newnewData[4::8] = numpy.left_shift(numpy.bitwise_and(newData[2::5], 255), 2) + numpy.right_shift(newData[3::5], 14)
    newnewData[5::8] = numpy.bitwise_and(numpy.right_shift(newData[3::5], 4), 1023)
    newnewData[6::8] = numpy.left_shift(numpy.bitwise_and(newData[3::5], 15), 6) + numpy.right_shift(newData[4::5], 10)
    newnewData[7::8] = numpy.bitwise_and(newData[4::5], 1023)
    #newnewData[117:] = [int(nSync>>22), int((nSync>>12) & 1023), int((nSync>>2 ) & 1023)]

    newnewData = [int(i) for i in newnewData] 

    newnewnewData = [0]*120
    newnewnewData[0:117]=newnewData[3:120]
    newnewnewData[117:] = [int(nSync>>22), int((nSync>>12) & 1023), int((nSync>>2 ) & 1023)]

   


    #print ":".join("{:04x}".format(c) for c in newData)
    #print ":".join("{:03x}".format(c) for c in newnewnewData)

    #print "{:03x}".format(newnewnewData[24])
    #raw_input("Press to continue")

    return newnewnewData

# Main Routine

def main(lArgs):
    plog = sys.stderr.write
    pout = sys.stdout.write

    g_nWords = 120
    #g_nNomRecLen = g_nWords * 2 + 12    # ---12 records/words for time tag
    g_nNomRecLen =  164  #changed for CAPERII 120 10 bit words packed into 75 16 bit words---  12 records/words for time tag
    nNomDataLen = g_nWords - 2  # less 2 sync words for 16 bit data and 3 words for 10 bit data
    nNomDataLen = 120  #changed for CAPER II   

    lFiles = lArgs[1:]
    if len(lFiles) != 1:
        prnHelp(sys.stderr)
        return 3

    for sFile in lFiles:
        if sFile.startswith('-h') or sFile.startswith('--help'):
            prnHelp(sys.stderr)
            return 0

    if not os.path.isdir(lArgs[1]):
        print "Directory does not exits."
        return 4

    os.chdir(lArgs[1])
    #masterfile='../output/' + lArgs[2]
    masterfiles=files=glob.glob("./output/*rawtm*.cdf") 	
    if not os.path.isfile(masterfiles[0]):
        print "No master rawtm cdf file in the output base directory."
        return 4


    lFiles = glob.glob('./input/*.tad')
    for sfile in lFiles:
        fIn = file(sfile, 'rb')
        pycdf.lib.set_backward(False)

        sFileHdr = fIn.read(10 + 12 + 22 + 260 + 12 + 4 + 4 + 4)
        pout("File:       %s\n" % sfile)
        pout("Signature:  %s\n" % sFileHdr[:10].strip('\x00'))
        pout("Version:    %s\n" % sFileHdr[10:22].strip('\x00'))
        pout("Time:       %s\n" % sFileHdr[22:44].strip('\x00'))
        pout("Config:     %s\n" % sFileHdr[44:304].strip('\x00'))
        pout("Source:     %s\n" % sFileHdr[304:316].strip('\x00'))
        pout("\n")

        timestmp = sFileHdr[22:44].strip('\x00')
        year = timestmp.split('/')[2][0:4]

        #dataFormat = ", %4d" * (nNomDataLen + 2) + "\n"

        sRec = fIn.read(g_nNomRecLen)
        ts=timestamp(year, S.unpack("3I", sRec[:12]))
        fIn.close()
        version=1
        print ts.isoformat()
        filedatetime=ts.isoformat().translate(None, '-:.')
        dirout='/home/srb/' #'./output/rawtm/'
        #fileout='trice003_k0_rawtm_'+pycdf.lib.tt2000_to_datetime(ts).isoformat()+'_v'+str(version).zfill(2)+'.cdf'
        fileout=str(project)+'_k0_rawtm_'+filedatetime[:15]+'_v'+str(version).zfill(2)+'.cdf'
	done = 0
        if not os.path.isfile(dirout+fileout):
            pout ("Creating File: %s\n" % fileout)
            iRec = 0L
            fIn = file(sfile, 'rb')
            sFileHdr = fIn.read(10+12+22+260+12+4+4+4)
            fcdf = pycdf.CDF(dirout+fileout, masterfiles[0])
            while True:

                sRec = fIn.read(g_nNomRecLen)
                if not sRec:
                    pout("ERROR: No record found\n")
                    fIn.close()
                    break
		(nMFcount, ) = S.unpack("<H", sRec[10:12])
                ts=timestamp(year, S.unpack("3I", sRec[:12]))
                if sRec[13:16] != g_sSync or ts == basedate:   #changed to start with 13 rather than 12 for CAPER II

                    pout("ERROR: Frame sync lost or bad time code at record %s\n" % str(iRec))
                    pout("\nTIME CODE: " + ts.isoformat()+"\n")
                    while sRec[13:16] != g_sSync or ts == basedate:  #changed to start with 13 rather than 12 for CAPERII
                        #keep reading until a syncword is found
                        sRec = fIn.read(g_nNomRecLen)
                        #print ''.join(['%s' %b for b in sRec])
                        if not sRec:
                            #fIn.close()
                            print "ERROR: No record found\n"
                            done = 1
                            break
			(nMFcount, ) = S.unpack("<H", sRec[10:12])
			ts=timestamp(year, S.unpack("3I", sRec[:12]))

                    if done == 1:
                        break
                    #iRec = 0L
                    #(nMFcount, ) = S.unpack("<H", sRec[10:12])
                    #ts=timestamp(year, S.unpack("3I", sRec[:12]))
                    #fileout='trice003_k0_rawtm_'+pycdf.lib.tt2000_to_datetime(ts).isoformat()+'_v'+str(version).zfill(2)+'.cdf'
		    #fileout=str(project)+'_k0_rawtm_'+ts.isoformat()+'_v'+str(version).zfill(2)+'.cdf'
                    #fcdf = pycdf.CDF(dirout+fileout, masterfiles[0])
		    #iRec += 1L
                #else:
		#(nMFcount, ) = S.unpack("<H", sRec[10:12])
                #ts = timestamp(year, S.unpack("3I", sRec[:12]))
		#pout(ts.isoformat()+"\n")
		fcdf['Epoch'][iRec:] = [ts]
                fcdf['minorframe'][iRec:] = [unpack10bitData(sRec, nNomDataLen)] #changed for CAPERII to 10bit version

                fcdf['sfid'][iRec:] = [nMFcount]
                if iRec % 100000 == 0:
                    pout("Number of Records written: %10d\n" % iRec)
		iRec += 1L
            fIn.close()
            pout("Number of Records written: %10d\n" % iRec)
            fcdf.close()
            print "CDF File Closed"

    return 0

##############################################################################
# kicker stub
if __name__ == "__main__":
    sys.exit(main(sys.argv))

