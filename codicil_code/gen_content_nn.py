#!/usr/bin/env python

import sys
from math import log
from math import sqrt
from operator import itemgetter

### Let's play with argparse ###
import argparse
argpar = argparse.ArgumentParser(description='Convert content graph.')
argpar.add_argument('cfile', type = argparse.FileType('r'), help = 'content list file name')
argpar.add_argument('tfile', type = argparse.FileType('r'), help = 'topology graph file name')
argpar.add_argument('ncount', type = int, help = 'node count')
argpar.add_argument('ccount', type = int, help = 'content vocabulary count')
argpar.add_argument('--hop', type = int, choices = (0, 1, 2), default = 0, help = 'search radius (hop)')
argpar.add_argument('-m', type = int, default = sys.maxint, help = 'tag per node')
argpar.add_argument('-k', type = int, default = 50, help = 'near neighbor count')
args = argpar.parse_args()

nodeCount = args.ncount
content_voc_size = args.ccount
h = args.hop
m = args.m
k = args.k

print >> sys.stderr, "h", h, "m", m, "k", k

print >> sys.stderr, "Reading topo graph"
inFile = args.tfile
inFile.readline()
adj = [set()] * (nodeCount+1)
lineCount = 0
for line in inFile :
    lineCount += 1
    adj[lineCount] = set(map(int, line.split()))
inFile.close()
print >> sys.stderr, "Reading topo graph done"

print >> sys.stderr, "Computing TF-IDF"
inFile = args.cfile
inFile.readline()
tfDict = [{}] * (nodeCount+1)
idfDict = [0] * (content_voc_size+1)
lineCount = 0
for line in inFile :
    lineCount += 1
    if lineCount % 10000 == 0 :
        print >> sys.stderr, lineCount
    nodeId = lineCount
    tf = {}
    for cid in map(int, line.split()) :
        if cid in tf :
            tf[cid] += 1
        else :
            tf[cid] = 1
            idfDict[cid] += 1
    tfDict[nodeId] = dict((cid, sqrt(tf)) for cid, tf in tf.iteritems())
inFile.close()
nodeCount = lineCount

idfDict = dict(((cid, 1+log(nodeCount * 1.0 / (idfDict[cid]+1))) for cid in xrange(1, content_voc_size+1)))
lenDict = [0.0] * (nodeCount+1)
vecDict = [{}] * (nodeCount+1)
for nodeId in xrange(1, nodeCount+1) :
    tfVec = tfDict[nodeId]
    if nodeId % 10000 == 0 :
        print >> sys.stderr, nodeId
    tfidfVec = dict(((cid, tf * idfDict[cid]) for cid, tf in tfVec.iteritems()))
    if len(tfVec) <= m :
        vecDict[nodeId] = tfidfVec
    else :
        vecDict[nodeId] = dict(sorted(tfidfVec.iteritems(), key=itemgetter(1), reverse = True)[:m])
    lenDict[nodeId] = sqrt(sum(tfidf**2 for tfidf in vecDict[nodeId].itervalues()))
del idfDict, tfDict
print >> sys.stderr, "Computing TF-IDF done"

# Instead of checking if a node is in the 2-hop neighbor of not
# We should just scan the list and add elements if not computed before
print >> sys.stderr, "Searching for nearest neighbors"
outFile = open("{0}.h{3}.m{1}.k{2}".format(args.cfile.name, m, k, h), 'w')
for nodeId in xrange(1, nodeCount+1) :
    if nodeId % 100 == 0 :
        print >> sys.stderr, nodeId
    nodeTagDict = vecDict[nodeId]
    cosDict = {}

    sl1 = None
    sl2 = None

    if h == 0 :
        sl1 = xrange(1, nodeCount+1)
    elif h == 1 :
        sl1 = adj[nodeId]
    elif h == 2 :
        sl1 = adj[nodeId]
        sl2 = (adj[nbrId] for nbrId in adj[nodeId])

    for neighbor in sl1 :
        if nodeId == neighbor or lenDict[neighbor] == 0.0 :
            continue
        neighborTagDict = vecDict[neighbor]
        dotProduct = 0.0
        for tag in nodeTagDict :
            if tag in neighborTagDict :
                dotProduct += nodeTagDict[tag] * neighborTagDict[tag]
        dotProduct /= lenDict[neighbor]
        if dotProduct > 0 :
            pass
        else :
            dotProduct = -1
        cosDict[neighbor] = dotProduct

    if sl2 != None :
        for gen in sl2 :
            for neighbor in gen :
                if nodeId == neighbor or lenDict[neighbor] == 0.0 :
                    continue
                if neighbor not in cosDict :
                    neighborTagDict = vecDict[neighbor]
                    dotProduct = 0.0
                    for tag in nodeTagDict :
                        if tag in neighborTagDict :
                            dotProduct += nodeTagDict[tag] * neighborTagDict[tag]
                    dotProduct /= lenDict[neighbor]
                    if dotProduct > 0 :
                        pass
                    else :
                        dotProduct = -1
                    cosDict[neighbor] = dotProduct

    nnList = (v for v, cos in sorted(((v, cos) for v, cos in cosDict.iteritems() if cos != -1), key=itemgetter(1), reverse = True)[:k])
    print >> outFile, nodeId, " ".join(map(str, nnList))
outFile.close()
print >> sys.stderr, "Searching for nearest neighbors done"

