
# !/usr/bin/env python
__author__ = "Scott Bounds"
__date__ = "2014-09-01"
__version__ = "1.0.0"

import sys
import os
import os.path
import struct as S
import datetime as dt
import numpy

# os.environ["CDF_LIB"] = "~/CDF/lib"
# os.putenv('CDF_LIB', '/usr/local/cdf/lib')
import pycdf



import glob

# g_sSync = "\xfa\xf3\x34\x00"  #The real frame sync but the data is little endian so...
g_sSync = "\x00\x34\xf3\xfa"
version = 'v00'


##############################################################################
def prnHelp(fOut):
    fOut.write("""tad_file_convert Usage:

	tad_file_convert DIR FILE
    DIR: Directory of the data to convert	
    FILE CDF master file 

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
    msec = (nFrmHdr[1] >> 20 & 15) * 100 + (nFrmHdr[1] >> 16 & 15) * 10 + (nFrmHdr[1] >> 12 & 15)
    usec = (nFrmHdr[1] >> 8 & 15) * 100 + (nFrmHdr[1] >> 4 & 15) * 10 + (nFrmHdr[1] & 15)

    if hours > 23:
        doy += 1
        hours = hours - 24

    d = dt.datetime.strptime(year + ' ' + str(doy), '%Y %j')
    t = dt.time(hours, mins, secs, msec * 1000 + usec, tzinfo=None)
    tt = dt.datetime.combine(d, t)

    return pycdf.lib.datetime_to_tt2000(tt)


###############################################################################
def unpackData(sRec, nNomDataLen):
    'Unpacks the record into integer values for all data in a minor frame'
    nData = numpy.array(S.unpack("59I", sRec[16:]))
    nSync = numpy.array(S.unpack("I", sRec[12:16]))
    newData = [0] * nNomDataLen
    newData[0::2] = nData >> 16
    newData[1::2] = nData & 1023
    newData[nNomDataLen - 1:] = [nSync >> 22, (nSync >> 12) & 1023, (nSync >> 2) & 1023]
    return newData


# Main Routine

def main(lArgs):
    plog = sys.stderr.write
    pout = sys.stdout.write

    g_nWords = 120
    g_nNomRecLen = g_nWords * 2 + 12 # 252
    nNomDataLen = g_nWords - 2       # 250

    lFiles = lArgs[1:]
    if len(lFiles) != 2:
        prnHelp(sys.stderr)
        return 3

    for sFile in lFiles:
        if sFile.startswith('-h') or sFile.startswith('--help'):
            prnHelp(sys.stderr)
            return 0

    if not os.path.isdir(lArgs[1]):
        print
        "Directory does not exits."
        return 4

    os.chdir(lArgs[1])
    masterfile = '../output/' + lArgs[2]
    if not os.path.isfile(masterfile):
        print
        "No master file in this directory."
        return 4

    lFiles = glob.glob('*.tad')
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

        dataFormat = ", %4d" * (nNomDataLen + 2) + "\n"

        sRec = fIn.read(g_nNomRecLen)
        ts = timestamp(year, S.unpack("3I", sRec[:12]))

        version = 0

        dirout = '../output/rawtm/'
        fileout = 'caper_k0_rawtm_' + pycdf.lib.tt2000_to_datetime(ts).isoformat() + '_v' + str(version).zfill(
            2) + '.cdf'


        if not os.path.isfile(dirout + fileout):
            fIn.close()
            iRec = 0
            fIn = file(sfile, 'rb')
            sFileHdr = fIn.read(10 + 12 + 22 + 260 + 12 + 4 + 4 + 4)
            fcdf = pycdf.CDF(dirout + fileout, masterfile)
            while True:
                sRec = fIn.read(g_nNomRecLen)
                if not sRec:
                    fIn.close()
                    break
                if sRec[12:16] != g_sSync:
                    plog("ERROR: Frame sync lost at record %s\n" % str(iRec))
                    fcdf.close()
                    done = 0
                    while sRec[12:16] != g_sSync:
                        sRec = fIn.read(g_nNomRecLen)
                        if not sRec:
                            fIn.close()
                            done = 1
                            break
                    if done == 1:
                        done = 0
                        break
                    iRec = 0
                    ts = timestamp(year, S.unpack("3I", sRec[:12]))
                    fileout = 'caper_k0_rawtm_' + pycdf.lib.tt2000_to_datetime(ts).isoformat() + '_v' + str(
                        version).zfill(2) + '.cdf'
                    fcdf = pycdf.CDF(dirout + fileout, masterfile)

                (nMFcount,) = S.unpack("<H", sRec[10:12])

                fcdf['Epoch'][iRec:] = [timestamp(year, S.unpack("3I", sRec[:12]))]
                fcdf['minorframe'][iRec:] = [unpackData(sRec, nNomDataLen)]
                fcdf['sfid'][iRec:] = [nMFcount]
                if iRec % 10000 == 0:
                    pout("Number of Records written: %10d\n" % iRec)

                iRec += 1

            fIn.close()
            pout("Number of Records written: %10d\n" % iRec)
            fcdf.close()

    return 0


##############################################################################
# kicker stub
if __name__ == "__main__":
    sys.exit(main(sys.argv))