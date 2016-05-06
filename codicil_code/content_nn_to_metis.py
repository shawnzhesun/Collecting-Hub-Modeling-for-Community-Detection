#!/usr/bin/env python

import sys

args = sys.argv
usageStr = "{0} knn_file".format(args[0])
if len(args) < 2 :
    print >> sys.stderr, usageStr
    sys.exit(1)

print >> sys.stderr, "loading {0}".format(args[1])
adjListDict = {}
lineCount = 0
for line in open(args[1], 'r') :
    lineCount += 1
    if lineCount % 100000 == 0 :
        print >> sys.stderr, lineCount
    splits = map(int, line.split())
    vid = splits[0]
    if vid not in adjListDict :
        adjListDict[vid] = []
    for uid in splits[1:] :
        adjListDict[vid].append(uid)
        if uid not in adjListDict :
            adjListDict[uid] = []
        adjListDict[uid].append(vid)

maxMid = max(adjListDict.iterkeys())
nodeCount = maxMid
edgeCount = 0
print >> sys.stderr, "sorting"
lineCount = 0
for mid in xrange(1, maxMid+1) :
    lineCount += 1
    if lineCount % 100000 == 0 :
        print >> sys.stderr, lineCount
    if mid in adjListDict :
        adjListDict[mid] = sorted(list(set(adjListDict[mid])))
        edgeCount += len(adjListDict[mid])
edgeCount /= 2
print >> sys.stderr, "sorting done"

print >> sys.stderr, "nodeCount", nodeCount, "edgeCount", edgeCount
print >> sys.stderr, "writing to file"
outFile = open(args[1]+".conv", 'w')
print >> outFile, nodeCount, edgeCount
for mid in xrange(1, maxMid+1) :
    if mid not in adjListDict :
        print >> outFile
    else :
        print >> outFile, " ".join(map(str, adjListDict[mid]))
outFile.close()

