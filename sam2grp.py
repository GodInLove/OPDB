#! /usr/bin/python
import sys
import re

def main():
    fragmentSize = 0
    if len(sys.argv) == 3:
        fragmentSize = int(sys.argv[2])
    posDict, idDict = samReader(fragmentSize)
    startDict, covDict = readCounter(posDict, idDict)
    grpWriter(startDict, covDict)

def samReader(fragmentSize):
    """parses the sam file for read coordinates and IDs
    """
    print 'parsing sam file...'
    samFile = open(sys.argv[1],"r")
    parsed, count = 0, 0
    posDict, idDict = {}, {}
    old = []
    for line in samFile:
        column = line.split('\t')
        if len(column) < 11:
            if column[0] == '@SQ':
                posDict[column[1].split(':')[1]] = {0:{}, 1:{},'maxPos': 0}
            continue
        count += 1
        if (count % 1000000) == 0:
            parsed += 1
            print str(parsed) + ' million lines parsed'
        if column[6] in posDict and not column[6] == column[2]:
            continue
        decoded    = hexDecoder(int(column[1])) # decoded form of bitwise flag
        strand     = decoded[-5]
        ID         = column[0]
        cigar      = column[5]
        refName    = column[2]
        start, end = getPositions(cigar, int(column[3]))
        if decoded[-2] == 1: # paired end
            if old == []:
                old = [start, end, strand, ID]
            else:
                pairedStart     = min([start, end, old[0], old[1]])
                pairedEnd       = max([start, end, old[0], old[1]])
                if pairedEnd - pairedStart <= fragmentSize:
                    posDict[refName], idDict = addIDs(pairedStart, pairedEnd, old[2], ID, posDict[refName], idDict)
                    if end > posDict[refName]['maxPos']: # get the maximum genome position with reads
                       posDict[refName]['maxPos'] = end
                old = []
        elif decoded[-1] == 0: # unpaired read (if paired end, only first read mapped)
            posDict[refName], idDict = addIDs(start, end, strand, ID, posDict[refName], idDict)
            if end > posDict[refName]['maxPos']: # get the maximum genome position with reads
                posDict[refName]['maxPos'] = end
    samFile.close()
    print 'sam-file successfully parsed'
    return posDict, idDict

def readCounter(posDict, idDict):
    """counts the number of reads for all covered positions
    """
    covDict   = {}
    startDict = {}
    for refName, refDict in posDict.items():
        if refDict['maxPos'] == 0:
            continue
        print 'calculating coverage for ' + refName + '...'
        covDict[refName]   = {0:{}, 1:{}, 'maxPos':refDict['maxPos']} # used for final count of coverage
        startDict[refName] = {0:{}, 1:{}, 'maxPos':refDict['maxPos']}
        for strand, starts in refDict.items():
            if strand != 'maxPos':
                for start, IDs in starts.items():
                    startCount = 0
                    for ID, end in IDs.items():
                        count       = float(1)/idDict[ID]
                        startCount += count
                        for position in range(min([start, end]),max([start, end])+1):
                            if position in covDict[refName][strand]:
                                covDict[refName][strand][position] += count
                            else:
                                covDict[refName][strand][position]  = count
                    startDict[refName][strand][start] = startCount
        print 'coverage for ' + refName +' successfully calculated.'
    return startDict, covDict

def grpWriter(startDict, covDict):
    """writes starts and coverage into grp-files
    """
    print 'writing grp-files...'
    for refName, refDict in startDict.items():
        filename  = re.findall('(.*)\.',sys.argv[1])[0]
        fwdFile = open(filename + '_' + refName + "_fwd.grp", "w")
        revFile = open(filename + '_' + refName + "_rev.grp", "w")
        for position in range(1,refDict['maxPos']+1):
            fwd = ["0","0"]
            rev = ["0","0"]
            if position in startDict[refName][0]:
                fwd[0] = str(float(startDict[refName][0][position]))
            if position in covDict[refName][0]:
                fwd[1] = str(float(covDict[refName][0][position]))
            if position in startDict[refName][1]:
                rev[0] = '-' + str(float(startDict[refName][1][position]))
            if position in covDict[refName][1]:
                rev[1] = '-' + str(float(covDict[refName][1][position]))
            fwdFile.write('\t'.join(fwd) + '\n')
            revFile.write('\t'.join(rev) + '\n')
        fwdFile.close()
        revFile.close()
    print 'grp-files successfully written'

def addIDs(start, end, strand, ID, posDict, idDict):
    """counts ID-frequencies and saves starting-positions
    """
    if strand == 1:
        start, end = end, start
    if ID in idDict:
        idDict[ID] += 1
    else:
        idDict[ID]  = 1
    if start in posDict[strand]:
        posDict[strand][start][ID] = end
    else:
        posDict[strand][start] = {ID:end}
    return posDict, idDict

def hexDecoder(hexNumber):
    """decodes bitwise FLAG, decoded starts with 0x400 -> 0x1
    """
    hexArray  = [1024,512,256,128,64,32,16,8,4,2,1]
    decoded   = [0,0,0,0,0,0,0,0,0,0,0]
    for i in range(len(hexArray)):
        if hexNumber - hexArray[i] >= 0:
            decoded[i] = 1
            hexNumber -= hexArray[i]
    return decoded 

def getPositions(cigar, pos1):
    """counts the alignment length and return start an stop
    """
    m = sum([int(x) for x in re.findall("(\d+)M",cigar)])
    d = sum([int(x) for x in re.findall("(\d+)D",cigar)])
    i = sum([int(x) for x in re.findall("(\d+)I",cigar)])
    s = sum([int(x) for x in re.findall("(\d+)S",cigar)])
    pos2 = pos1 + m + d
    return pos1, pos2
    
if __name__ == '__main__':
    main()
