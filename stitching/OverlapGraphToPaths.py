#!/usr/bin/env python

import networkx as nx
import sys
import argparse
import Overlap
#import IPython

ap = argparse.ArgumentParser(description="Find a path through overlap graphs")
ap.add_argument("ovp", help="Overlap file")
ap.add_argument("graph",help="Graph file")
ap.add_argument("paths", help="Output paths")
ap.add_argument("--tgraph", help="Graph with transitive edge removal", default=None)
args = ap.parse_args()


def OrderDestByOverlap(n,g,o, extend=False):
    adj = []
    for dest in g[n]:
        ovp = o[g.edge[n][dest]['index']]
        if extend == True and ovp.Extends(wiggle=2000) == False:
            continue
        dist = ovp.DistanceToStart()
        if dist != -1:
            adj.append((dist, dest))
    return sorted(adj)

def RemoveTips(g,l):


    nTips = 0
    toRemove = []
    for n in g.nodes():
        if g.in_degree(n) == 0 or g.out_degree(n) == 0:
            # Found a source or sink
            src = n
            srcLen = 0
            srcOut = n
            while g.in_degree(src) == 1 and g.out_degree(src) <= 1 and srcLen <= l:
                srcOut = src
                src = g.in_edges(src)[0]
                srcLen += 1

            dest = n
            destLen = 0
            destOut = n
            while g.out_degree(dest) == 1 and g.in_degree(dest) <= 1 and destLen + srcLen <= l:
                destOut = dest
                dest    = g[dest].keys()[0]
                destLen += 1

            if srcLen + destLen <= l:
                # Found a tip, remove the path from the src node to the
                # end of the dest
                if srcOut != src:
                    toRemove.append(srcOut)
                    while srcOut != n:
                        toRemove.append(srcOut)
                        src = srcOut
                        srcOut = g.in_edges(src)[0]

                dest = n
                while g.out_degree(dest) == 1 and g.in_degree(dest) <= 1:
                    toRemove.append(n)
                    destOut = g[dest].keys()[0]
                    dest = destOut

    g.remove_nodes_from(toRemove)
    print("removed {} tips".format(len(toRemove)))
    return nTips

def SoleOut(g,n):
    return g[n].keys()[0]

def PathToBranch(g, e):
    path = [e[0], e[1]]
    d = e[1]
    while g.in_degree(d) == 1 and g.out_degree(d) == 1:
        s = d
        d = SoleOut(g,d)
        path.append(d)
    return path


def GetSourceNodes(g):
    #
    # Return a list of source nodes in topological order.
    #
    nodes = []
    for n in nx.topological_sort(g):
        if g.in_degree(n) == 0:
            nodes.append(n)
    return nodes

BRANCH_OUT=0
BRANCH_IN=1
BRANCH_INOUT=2
def GetBranchingNodes(g, branch=BRANCH_OUT):
    forkNodes = []
    for n in g.nodes():
        if branch == BRANCH_OUT or branch == BRANCH_INOUT:
            if g.out_degree(n) == 2:
                forkNodes.append(n)
        if branch == BRANCH_IN or branch == BRANCH_INOUT:
            if g.in_degree(n) == 2:
                forkNodes.append(n)
    return forkNodes



def RemoveSimpleBulges(g, maxBulgeLength):
    forkNodes = GetBranchingNodes(g, BRANCH_OUT)
    for source in forkNodes:
        dests = g[source].keys()
        if len(dests) != 2:
            print("ERROR, should be on a branch")
        paths = [PathToBranch(g, (source, dest)) for dest in dests]
        if paths[0][-1] != paths[1][-1]:
            print("Path from " + source + " is not a bulge, " + str(len(paths[0])) + "\t" + str(len(paths[1])))
            continue
        print("Paths: ")
        print(str(paths[0]))
        print(str(paths[1]))

def DFSOverlaps(g,n,o,c,d, prevOverlap=None):
    adj = OrderDestByOverlap(n, g, o, extend=True)
    if len(adj) == 0:
        return (c, n)
    depths = [0]*len(adj)
    dest   = [None]*len(adj)

    for i in range(0,len(adj)):
        curOverlap = o[g.edge[n][adj[i][1]]['index']]
        #
        # If the middle regions of the alignments do not have proper overlap (the prev mid must
        # end before the cur mid).
        #
        import pdb
        pdb.set_trace()
        if prevOverlap is not None and prevOverlap.bMidOvp[1] >= curOverlap.aMidOvp[1]:
            sys.stderr.write("Skipping " + str(curOverlap) + "\n")
            continue
        points = Overlap.GetOverlapPoints(prevOverlap, curOverlap)
        if points[0] < points[1]:
            if d > 1:
                (depths[i], tempMaxDest) = DFSOverlaps(g,adj[i][1],o,c+1,d-1, curOverlap)
                dest[i] = adj[i][1]
            else:
                depths[i] = c+1
    maxDepth = max(depths)
    maxI =   depths.index(maxDepth)
    maxDest = dest[maxI]
    return (maxDepth, maxDest)




def GreedyPath(g,n,o):
    source = n
    path = [n]
    prevEdge = None
    curEdge  = None
    prevOverlap = None
    foundOverlap = False
    while len(g[n].keys()) > 0 and g.node[n]['visited'] == False:
        g.node[n]['visited'] = True
        adj = OrderDestByOverlap(n, g, o, extend=True)
        if len(adj) == 0:
            break


        foundOverlap = False
        if prevOverlap is not None:
            # Search 4 nodes ahead for maximum path.
            (depth, dest) = DFSOverlaps(g,n,o,0,4, prevOverlap)
            if dest is not None:
                curOverlap = o[g.edge[n][dest]['index']]
                foundOverlap = True
                n = dest
        else:
            curOverlap = o[g.edge[n][adj[0][1]]['index']]
            n = adj[0][1]
            foundOverlap = True
        if foundOverlap == True:
            path.append(n)
            foundOverlap = False
        else:
            print("ending search " + str(len(adj)))
        prevOverlap = curOverlap
    if foundOverlap == True:
        path.append(n)
    # Mark last node as visited
    g.node[n]['visited'] = True

    print("Path length: " + str(len(path)) + " Ended at: " + path[-1])
    return path

import pdb
def DFSLongestPath(g, n, lp, lpd):
    # at dest node

    if len(g[n]) == 0:
        lp[n] = 1
        lpd[n] = None
        return
    
    for d in g[n]:
        if lp[d] == -1:
            DFSLongestPath(g, d, lp, lpd)
        if lp[n] < lp[d] + 1:
            lp[n] = lp[d] + 1
            lpd[n] = d


def StoreLongestPath(g,n,lp,lpd):
    path = [n]
    dest = lpd[n]
    while dest is not None:
        n = lpd[n]
        path.append(n)
        dest = lpd[n]

    return path
                
def LongestPath(g):
    ts = nx.topological_sort(g)
    lp = { i: -1 for i in g.nodes() }
    lpd = { i: None for i in g.nodes() }

    for n in ts:
        if lp[n] == -1:
            DFSLongestPath(g, n, lp, lpd)

    longestPathLength = 0
    longestPathStart = -1
    for n in g.nodes():
        if lp[n] > longestPathLength:
            longestPathStart = n
            longestPathLength = lp[n]
    path = StoreLongestPath(g, longestPathStart, lp, lpd)
#    print("Stored longest path " + str(longestPathLength) + " stored " + str(len(path)))
#    print(path)
    return path
        
    

def GreedyPaths(g, o):
    nx.set_node_attributes(g, 'visited', {n: False for n in g.nodes()})
    srcNodes = GetSourceNodes(g)
    paths = [GreedyPath(g,n,o) for n in srcNodes]
    for n in g.nodes():
        if 'visited' in g.node[n]:
            del g.node[n]['visited']

    return paths

def LinearPaths(g):
    srcNodes = GetSourceNodes(g)
    nx.set_node_attributes(g, 'visited', {n: False for n in g.nodes()})
    optPaths = []
    for src in srcNodes:
        branchStack = [(src, 1)]
        branchPrev = [src]
        terminalNodes = []
        terminalNodeLengths = []
        path = []
        optPath = [src]
        tmp = src

        while len(branchStack) > 0:
            (node, index) = branchStack.pop()
            print("Branching - node: {}, index: {}, cur path: {}".format(node, index, len(path)))
            path = path[0:index]
            path.append(node)
            g.node[node]['visited'] = True
            # if this is a branching path, add dest to the Queue
            # process nonbranching paths.

            while g.out_degree(node) == 1:
                g.node[node]['visited'] = True
                node = SoleOut(g, node)
                path.append(node)
                if g.node[node]['visited'] == True:
                    break



            #
            # First case, the search ended at a branching node
            #
            if g.node[node]['visited'] == False and g.out_degree(node) > 1:
                for dest in g[node].keys():
                    branchStack.append((dest, len(path)))
                continue

            # Reached either the end of a path, or a branch.  Handle both separately
            if g.out_degree(node) == 0:
                if len(path) > len(optPath):
                    optPath = path
        print("source: " + src + " has path of length " + str(len(optPath)))
        optPaths.append(optPath)
    return optPaths




def RemoveTransitiveEdges(g, o):
    #
    # perform simple a->b, b->c, a-/->c transitive edge removal.
    # g - digraph
    # o - overlaps, indexed by g[node]['index']
    #
    transitive = []
    nThird = 0
    nFourth = 0
    for n in nx.topological_sort(g):
        inPlay = { dest for dest in g[n].keys() }
        adj = OrderDestByOverlap(n, g, o)
        for (dist, neighbor) in adj:
            for second in g[neighbor].keys():
                if second in inPlay:
                    transitive.append((n,second))
                for third in g[second].keys():
                    if third in inPlay:
                        transitive.append((n,third))
                        nThird+=1
                    for fourth in g[third].keys():
                        if fourth in inPlay:
                            transitive.append((n,fourth))
                            nFourth+=1
    g.remove_edges_from(transitive)

g = nx.read_gml(args.graph)
overlapFile = open(args.ovp)
pathFile = open(args.paths,'w')
overlaps = []
i=1
for line in overlapFile:
    # hack to get around blank lines
    if len(line) > 1:
        overlaps.append(Overlap.Overlap(line))
    i+=1



#RemoveTransitiveEdges(g, overlaps)
#RemoveTips(g,1)

components = nx.weakly_connected_components(g)
paths = []
for comp in components:

    subgraph = g.subgraph(comp)
#    greedyPaths=GreedyPaths(subgraph, overlaps)
    longestPath=LongestPath(subgraph)
    paths.append(longestPath)

for p in paths:
    pathFile.write("\t".join(p) + "\n")

if args.tgraph is not None:
    nx.write_gml(g, args.tgraph)
