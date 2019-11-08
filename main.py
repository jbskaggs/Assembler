import copy
import random
import re
import string
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


def hamming_distance(dna0, dna1):
    count = 0
    for i in range(len(dna0)):
        if dna0[i] != dna1[i]:
            count += 1
    return count


def generate_the_d_neighborhood(dna, d):
    if d == 0:
        return {dna}
    if len(dna) == 1:
        return ['A', 'C', 'G', 'T']
    neighborhood = []
    suffix_neighbors = generate_the_d_neighborhood(dna[1:], d)
    for text in suffix_neighbors:
        if hamming_distance(dna[1:], text) < d:
            for nucleotide in ['A', 'C', 'G', 'T']:
                neighborhood += [nucleotide + text]
        else:
            neighborhood += [dna[0] + text]
    return neighborhood


def dfs(nodes, node, lastNode, path):
    for i in range(len(node.ys)):
        if node.visited[i] == 0:
            node.visited[i] = 1
            dfs_result = dfs(nodes, nodes[node.ys[i]], lastNode, path + [node.ys[i]])
            if dfs_result[0] == lastNode:  # and len(dfs_result[1]) == num_vals * 3 + 1:
                return dfs_result
            node.visited[i] = 0
    return node.x, path


def find_an_eulerian_cycle(nodes):
    # start_node = int(random() * len(nodes))
    start_node = 0
    _, cycle = dfs(nodes, nodes[start_node], start_node, [start_node])
    i = 0
    while i < len(cycle):
        p = cycle[i]
        if len(nodes[p].ys) > 1:
            for j in range(len(nodes[p].ys)):
                k = nodes[p].ys[j]
                if nodes[p].visited[j] == 0:
                    _, c = dfs(nodes, nodes[p], nodes[p].x, [nodes[p]])
                    cycle = cycle[:i+1] + c[1:] + cycle[i+1:]
                    break
        i += 1
    return cycle


def find_an_eulerian_path(nodes):
    start_node = 0
    last_node = 0
    for node in nodes:
        if len(node.ys) > len(node.inv_ys):
            start_node = node.x
            break
    for node in nodes:
        if len(node.ys) < len(node.inv_ys):
            last_node = node.x
            break

    _, e_path = dfs(nodes, nodes[start_node], last_node, [start_node])
    i = 0
    while i < len(e_path):
        p = e_path[i]
        if len(nodes[p].ys) > 1:
            for j in range(len(nodes[p].ys)):
                k = nodes[p].ys[j]
                if nodes[p].visited[j] == 0:
                    _, c = dfs(nodes, nodes[p], nodes[p].x, [nodes[p]])
                    e_path = e_path[:i+1] + c[1:] + e_path[i+1:]
                    break
        i += 1
    return e_path



def HammingDistance(string1, string2):
    hammingCount = 0
    for c1, c2 in zip(string1, string2):
        if c1 != c2:
            hammingCount +=1
    return hammingCount


def Neighbors(pattern, d):
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}
    neighborhood = set()
    suffixneighbors = Neighbors(pattern[1:], d)
    for st in suffixneighbors:
        if HammingDistance(pattern[1:], st) < d:
            neighborhood.add("A{0}".format(st))
            neighborhood.add("C{0}".format(st))
            neighborhood.add("G{0}".format(st))
            neighborhood.add("T{0}".format(st))
        else:
            neighborhood.add(pattern[0]+st)
    return neighborhood


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

                # in real sequencing graphs, we need to account for error
                if curNodekey in nodeList.keys():
                    curEdges = nodeList[curNodekey]
                    found = True
                else:
                    #account for error. We decided to do this by generating the D neighborhood of the key
                    neighbors = Neighbors(curNodekey, 1)
                    for neigh in neighbors:
                        if neigh in nodeList.keys():
                            curEdges = nodeList[neigh]
                            found = True
                            break
                    found = False

                if found:
                    curPath.append(curNodekey)
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

    k = 21

    fasta_sequences = SeqIO.parse(open('real.error.small.fasta'), 'fasta')

    dnas = []
    for record in fasta_sequences:
        dnas += [str(record.seq)]

    kmers = []
    for dna in dnas:
        for i in range(len(dna) - k):
            kmers += [dna[i:i + k]]

    adj = adjacencyListBrujn(kmers)

    nodelist, edgeCount, start, end = buildGraphWithFalseEdge(adj)

    for node in nodelist.keys():
        for i in range(len(nodelist[node])):
            if nodelist[node][i].Node not in nodelist.keys():
                for j in range(1, k):
                    for neighbor in generate_the_d_neighborhood(nodelist[node][i].Node, 1):
                        if neighbor in nodelist.keys():
                            nodelist[node][i].Node = neighbor
                            break
                    if nodelist[node][i].Node in nodelist.keys():
                        break


    answer = findPath(nodelist, edgeCount, start)
    answer = formatPath(answer, end)
    answer = stringFromPath(answer)
    print(answer)
