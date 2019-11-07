import copy
import random
import re
from Bio import SeqIO


class Edge:

    def __init__(self, NodePoint):
        self.Node = NodePoint
        self.Node_r = self.reverse_complement(NodePoint)
        self.visited = False
        self.count = 1

    def __eq__(self, other):
        return self.Node == other.Node or self.Node == other.Node_r

    def add_one(self):
        self.count += 1

    @staticmethod
    def reverse_complement(dna):
        dna = re.sub('A', 'B', dna)
        dna = re.sub('T', 'A', dna)
        dna = re.sub('B', 'T', dna)
        dna = re.sub('C', 'B', dna)
        dna = re.sub('G', 'C', dna)
        dna = re.sub('B', 'G', dna)
        return dna[::-1]


def Suffix(pattern):
    return pattern[:(len(pattern) - 1)]


def Prefix(pattern):
    return pattern[1:]


def adjacencyList(kmerlist):
    adjacencies = {}
    for i in range(len(kmerlist)):
        for x in range(len(kmerlist)):
            if (i != x) & (Suffix(kmerlist[x]) == Prefix(kmerlist[i])):
                adjacencies[kmerlist[i]] = kmerlist[x]
    return adjacencies


# edited to use Edge Class
def adjacencyListBrujn(kmerlist):
    adjacencies = {}
    for i in range(len(kmerlist)):
        if Suffix(kmerlist[i]) in adjacencies.keys():
            adjacencies[Suffix(kmerlist[i])].append(Edge(Prefix(kmerlist[i])))
        else:
            adjacencies[Suffix(kmerlist[i])] = [Edge(Prefix(kmerlist[i]))]
    return adjacencies


# kmer string version
def buildGraphWithFalseEdge(adlist):
    edgeCount = 0
    edgelist = []
    end = ""
    start = ""

    for edge in adlist.values():
        for val in edge:
            edgeCount += 1
            edgelist.append(val)

    # find start. Needs to be a node with 1 more outdegree than indegree
    for key in adlist.keys():
        indegreeCount = 0
        for lol in edgelist:
            if lol.Node == key:
                indegreeCount += 1
        if len(adlist[key]) > indegreeCount:
            start = key
            break

    # find end. Is node with one more indegree than outdegree
    # we are going to try a shortcut and just find the end of the path
    for edge in edgelist:
        if edge.Node not in adlist.keys():
            end = edge.Node
            break

    # add false edge
    adlist[end] = [Edge(start)]
    edgeCount += 1

    return adlist, edgeCount, start, end


def findPath(nodeList, edgeCount, start):
    pathLength = 1

    curNodekey = start

    curEdges = nodeList[curNodekey]
    # curEdges[0].visited = True

    # keeping track of which visited nodes have unvisited edges
    hasMoreEdges = None
    overallPath = []

    # this cycle iteration
    curPath = [curNodekey]

    # keeping track of when to stop by checking number of edges
    while edgeCount > len(curPath):
        # find unvisited edge from current node
        found = False
        for ind in range(len(curEdges)):
            if not curEdges[ind].visited:
                curNodekey = curEdges[ind].Node
                curEdges[ind].visited = True

                # in eulerian paths, this is to check if we are at the end node
                # if curNodekey not in nodeList.keys():
                #     end = curNodekey
                #     curPath.append(curNodekey)
                #     break

                curEdges = nodeList[curNodekey]
                curPath.append(curNodekey)
                found = True
                break
        if not found:  # we are stuck
            # we need to "reset" the ant but keep the progress we have made
            # first find node with unexplored edges from along the path
            # iterate through node keys and find node with untouched path
            found = False
            for ind in range(len(curPath)):
                curNodekey = curPath[ind]
                curEdges = nodeList[curNodekey]
                # iterate through this nodes edges and see if you can find an unexplored node
                for ind2 in range(len(curEdges)):
                    if not curEdges[ind2].visited:
                        # we found one
                        found = True
                        # restructure overall path to fit information we have found so far
                        hlep = 0
                        overallPath = []

                        # restructure algorithm.
                        while hlep < len(curPath):
                            overallPath.append(curPath[ind])
                            hlep += 1
                            ind += 1
                            if ind >= len(curPath):
                                ind = 1
                        curPath = overallPath
                        break
                if found:
                    break
            # if we iterate through all the nodes in the path and don't find an index that will work we screwed
            if not found:
                print("No cycle found")
                return None

    return curPath


def formatPath(cycleList, end):
    cycleList = [str(i) for i in cycleList]
    foundEnd = False
    endpath = []
    startpath = []
    for el in cycleList:
        if el == end:
            foundEnd = True
            endpath.append(el)
        elif foundEnd:
            startpath.append(el)
        else:
            endpath.append(el)

    startpath.extend(endpath)
    return startpath


def stringFromPath(pathList):
    dna = pathList[0]
    for i in range(1, len(pathList)):
        dna += pathList[i][(len(pathList[i]) - 1)]

    return dna

if __name__ == '__main__':
    # with open("dnastring.txt") as dnaFile:
        # kmers = dnaFile.read().rstrip()
        # kmers = kmers.split("\n")
        # kmers = kmers[1:]

    fasta_sequences = SeqIO.parse(open('synthetic.example.noerror.small.fasta'), 'fasta')

    kmers = []
    for record in fasta_sequences:
        kmers += [str(record.seq)]

    # kmers = {}
    # for dna in dnas:
    #     for i in range(len(dna) - k):
    #         new_node = Node(dna[i:i + k])
    #         kmers

    adj = adjacencyListBrujn(kmers)

    nodelist, edgeCount, start, end = buildGraphWithFalseEdge(adj)
    answer = findPath(nodelist, edgeCount, start)
    answer = formatPath(answer, end)
    answer = stringFromPath(answer)
    print(answer)
