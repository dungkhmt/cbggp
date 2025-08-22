from math import log2
import random
import sys
import time
import sys

checkPrint = False

def readline():
    return sys.stdin.readline()

class VertexForGreedyAlgorithm:
    def __init__(self, id: int):
        self.id = id
        self.degree = 0
        self.frequency = 0
        self.inSolution = False
        self.tabuTenure = 0

class GraphForGreedyAlgorithm:
    def __init__(self, numberVertices: int, vertices: list, edges: list):
        self.numberVertices = numberVertices
        self.numberEdges = len(edges)
        self.vertices = vertices
        self.edges = edges
        self.adjacency = [[] for _ in range(numberVertices)]
        self.connected = [[False for _ in range(numberVertices)] for _ in range(numberVertices)]
        for edge in self.edges:
            u, v = edge
            self.adjacency[u].append(v)
            self.adjacency[v].append(u)
            self.connected[u][v] = True
            self.connected[v][u] = True
        self.density = self.numberEdges / (numberVertices * (numberVertices - 1) / 2) if numberVertices > 1 else 0

    def findClique(self):
        solution = []
        k = 1
        while k <= self.numberVertices:
            newSolution = self.amts(k, k * 4, k * 10)
            if len(newSolution) == k:
                solution = newSolution
                k += 1
            else:
                break
        return solution

    # adaptive multistart tabu search
    def amts(self, k: int,  L: int, iterMax: int):
        currentSolution = self.initialize(k)
        self.updateSolution(currentSolution)
        iter = 0
        while iter < iterMax:
            newSolution, iter = self.tabuSearch(currentSolution, k, L, iter)
            self.updateSolution(newSolution)
            costOfNewSolution = self.getCostOfSolution(newSolution)
            if costOfNewSolution == k * (k - 1) // 2:
                return newSolution
            currentSolution = self.frequencyBasedInitialize(k)
            self.updateSolution(currentSolution)
        return []  # Return an empty solution if no solution is found within the iterations
    
    def initialize(self, k: int):
        for u in range(self.numberVertices):
            self.vertices[u].degree = 0
            self.vertices[u].frequency = 0
            self.vertices[u].inSolution = False
            self.vertices[u].tabuTenure = 0
        solution = []
        for i in range(k):
            maxDegree = 0
            for v in range(self.numberVertices):
                if not self.vertices[v].inSolution and self.vertices[v].degree > maxDegree:
                    maxDegree = self.vertices[v].degree
            listOfVertices = [v for v in range(self.numberVertices) if not self.vertices[v].inSolution and self.vertices[v].degree == maxDegree]
            u = listOfVertices[random.randint(0, len(listOfVertices) - 1)]
            solution.append(u)
            self.vertices[u].inSolution = True
            for v in self.adjacency[u]:
                self.vertices[v].degree += 1
        return solution
    
    def tabuSearch(self, solution: list, k: int, L: int, iter: int):
        costOfSolution = self.getCostOfSolution(solution)
        bestSolution = solution[:]
        costOfBestSolution = self.getCostOfSolution(bestSolution)
        i = 0
        while i < L:
            iter += 1
            u, v, delta = self.chooseSwapMove(k, iter, costOfSolution)
            if (u < 0 or v < 0):
                break
            self.updateSolutionBySwap(solution, u, v)
            costOfSolution += delta
            self.updateTabuList(u, v, k, iter, costOfSolution)
            if costOfSolution == k * (k - 1) // 2:
                return solution, iter
            if costOfSolution > costOfBestSolution:
                bestSolution = solution[:]
                costOfBestSolution = costOfSolution
                i = 0
            else:
                i += 1
        return bestSolution, iter
    
    # update tabu list
    def updateTabuList(self, u: int, v: int,  k: int, iter: int, costOfSolution: int):
        l =  min(k * (k - 1) // 2 - costOfSolution, 10)
        C = max(int(k / 40), 6)
        self.vertices[u].tabuTenure = iter + l + random.randint(0, C)
        self.vertices[v].tabuTenure = iter + int(l * 0.6 + random.randint(0, int(C * 0.6)))

    # update the solution by swapping vertices
    def updateSolutionBySwap(self, solution: list, u: int, v: int):
        solution.remove(u)
        self.vertices[u].inSolution = False
        for w in self.adjacency[u]:
            self.vertices[w].degree -= 1
        self.vertices[u].frequency += 1
        solution.append(v)
        self.vertices[v].inSolution = True
        for w in self.adjacency[v]:
            self.vertices[w].degree += 1
        self.vertices[v].frequency += 1

    def chooseSwapMove(self, k: int, iter: int, costOfSolution: int):
        A = []
        B = []
        minInSolution = k
        maxOutSolution = 0
        for u in range(self.numberVertices):
            if self.vertices[u].tabuTenure < iter:
                if self.vertices[u].inSolution:
                    minInSolution = min(minInSolution, self.vertices[u].degree)
                else:
                    maxOutSolution = max(maxOutSolution, self.vertices[u].degree)
        for u in range(self.numberVertices):
            if self.vertices[u].tabuTenure < iter:
                if self.vertices[u].inSolution and self.vertices[u].degree == minInSolution:
                    A.append(u)
                elif not self.vertices[u].inSolution and self.vertices[u].degree == maxOutSolution:
                    B.append(u)
        if not A or not B:
            return -1, -1, 0
        T = self.bestSwapMoves(A, B)
        l = k * (k - 1) // 2 - costOfSolution
        p = min((l + 2) / self.numberEdges, 0.1)
        value = random.random()
        if (maxOutSolution <= minInSolution or (not T and maxOutSolution == minInSolution + 1)) and value <= p:
            return self.pmsr(k, iter)
        if T:
            u, v = T[random.randint(0, len(T) - 1)]
            delta = maxOutSolution - minInSolution
            return u, v, delta
        else:
            u = A[random.randint(0, len(A) - 1)]
            v = B[random.randint(0, len(B) - 1)]
            delta = maxOutSolution - minInSolution - 1
            return u, v, delta
    
    def bestSwapMoves(self, A: list, B: list):
        if not A or not B:
            return []
        T = []
        for u in A:
            for v in B:
                if not self.connected[u][v]:
                    T.append((u, v))
        return T

    # probabilistic move selection rule
    def pmsr(self, k: int, iter: int):
        listOfVertices = [v for v in range(self.numberVertices)]
        random.shuffle(listOfVertices)
        u = -1
        v = -1
        for w in listOfVertices:
            if self.vertices[w].tabuTenure >= iter:
                continue
            if u < 0 and self.vertices[w].inSolution:
                u = w
            if v < 0 and not self.vertices[w].inSolution and self.vertices[w].degree < int(k * self.density):
                v = w
            if u >= 0 and v >= 0:
                break
        delta = self.vertices[v].degree - self.vertices[u].degree - (1 if self.connected[u][v] else 0)
        return u, v, delta

    def updateSolution(self, solution: list):
        for id in range(self.numberVertices):
            self.vertices[id].degree = 0
            self.vertices[id].frequency = 0
            self.vertices[id].inSolution = False
            self.vertices[id].tabuTenure = 0
        for u in solution:
            self.vertices[u].inSolution = True
            for v in self.adjacency[u]:
                self.vertices[v].degree += 1

    def getCostOfSolution(self, solution: list) -> int:
        cost = 0
        for u in solution:
            cost += self.vertices[u].degree
        return cost / 2
    
    def frequencyBasedInitialize(self, k: int):
        self.refreshSolution(k)
        solution = []
        minFrequency = self.numberVertices
        for v in self.vertices:
            if v.frequency < minFrequency:
                minFrequency = v.frequency
        listVerticesMinFrequency = [i for i in range(self.numberVertices) if self.vertices[i].frequency == minFrequency]
        u = listVerticesMinFrequency[random.randint(0, len(listVerticesMinFrequency) - 1)]
        solution.append(u)
        self.vertices[u].inSolution = True
        for v in self.adjacency[u]:
            self.vertices[v].degree += 1
        while len(solution) < k:
            minDegree = self.numberVertices
            minFrequencyOfMinDegree = self.numberVertices
            for v in range(self.numberVertices):
                if not self.vertices[v].inSolution:
                    if self.vertices[v].degree < minDegree:
                        minDegree = self.vertices[v].degree
                        minFrequencyOfMinDegree = self.vertices[v].frequency
                    elif self.vertices[v].degree == minDegree and self.vertices[v].frequency < minFrequencyOfMinDegree:
                        minFrequencyOfMinDegree = self.vertices[v].frequency
            listOfVertices = [i for i in range(self.numberVertices) if not self.vertices[i].inSolution and self.vertices[i].degree == minDegree and self.vertices[i].frequency == minFrequencyOfMinDegree]
            u = listOfVertices[random.randint(0, len(listOfVertices) - 1)]
            solution.append(u)
            self.vertices[u].inSolution = True
            for v in self.adjacency[u]:
                self.vertices[v].degree += 1
        return solution
    
    def refreshSolution(self, k: int):
        for id in range(self.numberVertices):
            self.vertices[id].degree = 0
            self.vertices[id].frequency %= k
            self.vertices[id].inSolution = False

class Vertex:
    def __init__(self, id: int):
        self.id = id
        self.positionInListSorted = id

class Graph:
    def __init__(self, numberVertices: int, vertices: list, adjacency: list):
        self.numberVertices = numberVertices
        self.vertices = vertices
        self.adjacency = adjacency
        self.connected = [[False] * numberVertices for _ in range(numberVertices)]
        # for i in range(1, len(self.vertices)):
        #     if self.vertices[i].positionInListSorted < self.vertices[i - 1].positionInListSorted:
        #         print("Error: Input vertices are not sorted correctly.")
    
    def checkSolution(self, solution: list):
        if len(solution) == 0:
            return False
        for i in range(len(solution)):
            for j in range(i + 1, len(solution)):
                if not self.adjacency[solution[i]] & (1 << solution[j]):
                    return False
        return True

    # lấy tập hợp các đỉnh kề với nút u trong subgraph hiện tại
    def listVerticesOfNeighbors(self, adjacencyOfU: int):
        listOfVertices = []
        while adjacencyOfU > 0:
            lowestBit = adjacencyOfU & -adjacencyOfU
            v = lowestBit.bit_length() - 1
            listOfVertices.append(v)
            adjacencyOfU ^= lowestBit
        return listOfVertices

def degSort(graph: Graph, edges: list):
    degree = [0] * graph.numberVertices
    for edge in edges:
        u, v = edge
        degree[u] += 1
        degree[v] += 1
    graph.vertices.sort(key=lambda v: degree[v.id], reverse=True)
    positions = [0] * graph.numberVertices
    for i in range(graph.numberVertices):
        graph.vertices[i].positionInListSorted = i
        positions[graph.vertices[i].id] = i
    graph.adjacency = [0 for _ in range(graph.numberVertices)]
    graph.connected = [[False] * graph.numberVertices for _ in range(graph.numberVertices)]
    for i in range(len(edges)):
        edges[i] = (positions[edges[i][0]], positions[edges[i][1]])
        u, v = edges[i]
        graph.adjacency[u] |= (1 << v)
        graph.adjacency[v] |= (1 << u)
        graph.connected[u][v] = True
        graph.connected[v][u] = True
    # for i in range(graph.numberVertices):
    #     print(f"Vertex {self.vertices[i].id}:", end=' ')
    #     for j in self.adjacency[i]:
    #         print(f"{j} ", end='')
    #     print()
    # print();

def buildInitialUpperBound(graph: Graph):
    upperBound = [1] * graph.numberVertices
    for u in range(graph.numberVertices):
        adjacencyOfU = graph.adjacency[u]
        while adjacencyOfU > 0:
            lowestBit = adjacencyOfU & -adjacencyOfU
            v = lowestBit.bit_length() - 1
            if v >= u:
                break
            upperBound[u] = max(upperBound[u], upperBound[v] + 1)
            adjacencyOfU ^= lowestBit
    return upperBound

class GraphBitString:
    def __init__(self, numberVertices: int, connected: list, vertices: list):
        self.numberVertices = numberVertices
        self.connected = connected
        self.vertices = vertices
        self.vertices.sort()
        self.positionInVertices = {}
        for indexOfU in range(len(self.vertices)):
            self.positionInVertices[self.vertices[indexOfU]] = indexOfU
        self.bitStringOfGraph = (1 << len(self.vertices)) - 1
        self.adjacency = [0] * len(self.vertices)
        for indexOfU in range(len(self.vertices)):
            u = self.vertices[indexOfU]
            for indexOfV in range(indexOfU + 1, len(self.vertices)):
                v = self.vertices[indexOfV]
                if not self.connected[u][v]:
                    self.adjacency[indexOfU] |= (1 << indexOfV)

    def iseq(self, numberColors: int, vertices: list):
        # vertices.sort(key=lambda u: self.numberNeighbors[u])
        colorOfVertices = [-1] * self.numberVertices
        bitStringOfVertices = 0
        for u in vertices:
            bitStringOfVertices |= (1 << self.positionInVertices[u])
        for color in range(numberColors):
            if bitStringOfVertices == 0:
                return False, []
            copyBitStringOfVertices = bitStringOfVertices
            while copyBitStringOfVertices > 0:
                lowestBit = copyBitStringOfVertices & -copyBitStringOfVertices
                indexOfU = lowestBit.bit_length() - 1
                colorOfVertices[self.vertices[indexOfU]] = color
                bitStringOfVertices ^= (1 << indexOfU)
                copyBitStringOfVertices &= self.adjacency[indexOfU]
        if bitStringOfVertices == 0:
            return False, []
        return True, colorOfVertices
    
    def filt(self, numberColors: int, vertex: int, lastVertexOfColorHasVertex: list, verticesOfNewSubGraph: list, isInNewPrunedSet: list):
        bitStringOfVertices = 0
        for u in verticesOfNewSubGraph:
            if isInNewPrunedSet[u] and self.connected[vertex][u]:
                bitStringOfVertices |= (1 << self.positionInVertices[u])
        for color in range(numberColors):
            if bitStringOfVertices == 0:
                return False
            copyBitStringOfVertices = bitStringOfVertices
            lowestBit = copyBitStringOfVertices & -copyBitStringOfVertices
            indexOfU = lowestBit.bit_length() - 1
            lastVertex = lastVertexOfColorHasVertex[self.vertices[indexOfU]]
            bitStringOfVertices ^= (1 << indexOfU)
            copyBitStringOfVertices &= self.adjacency[indexOfU]
            while copyBitStringOfVertices > 0:
                lowestBit = copyBitStringOfVertices & -copyBitStringOfVertices
                indexOfU = lowestBit.bit_length() - 1
                u = self.vertices[indexOfU]
                if u < lastVertex:
                    copyBitStringOfVertices &= self.adjacency[indexOfU]
                else:
                    isInNewPrunedSet[u] = False
                    copyBitStringOfVertices ^= (1 << indexOfU)
                bitStringOfVertices ^= (1 << indexOfU)
        return True

class SoftClause:
    def __init__(self, bitStringOfVertices: int):
        self.bitStringOfVertices = bitStringOfVertices
        self.state = 0
        self.numberOfUndefinedLiterals = 0
        self.isInProperSet = False

class GraphSAT:
    def __init__(self, bitStringSubGraph: GraphBitString, vertices: list):
        self.maxNumber = bitStringSubGraph.numberVertices
        self.numberVertices = len(vertices)
        self.vertices = vertices
        bitStringOfVertices = 0
        for u in vertices:
            bitStringOfVertices |= (1 << bitStringSubGraph.positionInVertices[u])
        # self.isVertexInGraph = [False] * numberVertices
        # for vertex in vertices:
        #     self.isVertexInGraph[vertex] = True
        self.positionInVertices = {}
        self.indexOfSoftClause = {}
        for indexOfU in range(len(vertices)):
            self.positionInVertices[vertices[indexOfU]] = indexOfU
            self.indexOfSoftClause[indexOfU] = -1
        self.adjacency = [0 for _ in range(self.numberVertices)]
        for indexOfU in range(len(vertices)):
            u = vertices[indexOfU]
            adjacencyOfU = bitStringSubGraph.adjacency[bitStringSubGraph.positionInVertices[u]] & bitStringOfVertices
            while adjacencyOfU > 0:
                lowestBit = adjacencyOfU & -adjacencyOfU
                indexOfV = self.positionInVertices[bitStringSubGraph.vertices[lowestBit.bit_length() - 1]]
                self.adjacency[indexOfU] |= (1 << indexOfV)
                self.adjacency[indexOfV] |= (1 << indexOfU)
                adjacencyOfU ^= lowestBit
        self.softClauses = []

    def addSoftClause(self, bitStringOfVertices: int):
        if bitStringOfVertices == 0:
            return
        clauseIndex = len(self.softClauses)
        self.softClauses.append(SoftClause(bitStringOfVertices))
        while bitStringOfVertices > 0:
            lowestBit = bitStringOfVertices & -bitStringOfVertices
            indexOfU = lowestBit.bit_length() - 1
            self.indexOfSoftClause[indexOfU] = clauseIndex
            bitStringOfVertices ^= lowestBit

    def unitPropagation(self, vertex: int):
        # if not self.isVertexInGraph[vertex]:
        #     print(f"Error: Vertex {vertex} is not in the graph.")
        #     return False
        self.stateOfVertex = [-1 for _ in range(self.numberVertices)]
        for softClause in self.softClauses:
            softClause.state = 0
            softClause.numberOfUndefinedLiterals = softClause.bitStringOfVertices.bit_count()
            softClause.isInProperSet = False
        indexOfVertex = self.positionInVertices[vertex]
        self.stateOfVertex[indexOfVertex] = 1
        stack = [(indexOfVertex, 1)]
        self.softClauses[self.indexOfSoftClause[indexOfVertex]].isInProperSet = True
        while stack:
            indexOfU, state = stack.pop()
            # if checkPrint and vertex == 29:
            #     print(f"Processing vertex {u} with state {state}.")
            if state == 1:
                adjacencyOfU = self.adjacency[indexOfU]
                while adjacencyOfU > 0:
                    lowestBit = adjacencyOfU & -adjacencyOfU
                    indexOfV = lowestBit.bit_length() - 1
                    if self.stateOfVertex[indexOfV] == -1:
                        self.updateStateOfVertex(indexOfV, 0, self.indexOfSoftClause[indexOfV])
                        stack.append((indexOfV, 0))
                    adjacencyOfU ^= lowestBit
            if self.indexOfSoftClause[indexOfU] >= 0:
                softClause = self.softClauses[self.indexOfSoftClause[indexOfU]]
                bitStringOfSoftClause = softClause.bitStringOfVertices
                if state == 1 and softClause.numberOfUndefinedLiterals > 0:
                    while bitStringOfSoftClause > 0:
                        lowestBit = bitStringOfSoftClause & -bitStringOfSoftClause
                        indexOfV = lowestBit.bit_length() - 1
                        if self.stateOfVertex[indexOfV] == -1:
                            self.updateStateOfVertex(indexOfV, 0, self.indexOfSoftClause[indexOfV])
                            stack.append((indexOfV, 0))
                        bitStringOfSoftClause ^= lowestBit
                if softClause.state == 0:
                    if softClause.numberOfUndefinedLiterals == 0:
                        # if checkPrint and vertex == 29:
                        #     print(f"Error: Soft clause {softClause.vertices} is not satisfied.")
                        softClause.isInProperSet = True
                        return False
                    if softClause.numberOfUndefinedLiterals == 1:
                        while bitStringOfSoftClause > 0:
                            lowestBit = bitStringOfSoftClause & -bitStringOfSoftClause
                            indexOfV = lowestBit.bit_length() - 1
                            if self.stateOfVertex[indexOfV] == -1:
                                self.updateStateOfVertex(indexOfV, 1, self.indexOfSoftClause[indexOfV])
                                stack.append((indexOfV, 1))
                                break
                            bitStringOfSoftClause ^= lowestBit
        return True

    def updateStateOfVertex(self, indexOfVertex: int, state: int, indexOfsoftClause: int):
        self.stateOfVertex[indexOfVertex] = state
        if indexOfsoftClause >= 0:
            self.softClauses[indexOfsoftClause].state += state
            self.softClauses[indexOfsoftClause].numberOfUndefinedLiterals -= 1
            if state == 1:
                self.softClauses[indexOfsoftClause].isInProperSet = True

    def removeVertex(self, u: int):
        indexOfU = self.positionInVertices[u]
        adjacencyOfU = self.adjacency[indexOfU]
        while adjacencyOfU > 0:
            lowestBit = adjacencyOfU & -adjacencyOfU
            indexOfV = lowestBit.bit_length() - 1
            self.adjacency[indexOfV] ^= (1 << indexOfU)
            adjacencyOfU ^= lowestBit
        clauseIndex = self.indexOfSoftClause[indexOfU]
        if clauseIndex >= 0:
            self.indexOfSoftClause[indexOfU] = -1
            self.softClauses[clauseIndex].bitStringOfVertices ^= (1 << indexOfU)

    def removeLastSoftClause(self):
        if not self.softClauses:
            return
        lastClause = self.softClauses.pop()
        bitStringOfVertices = lastClause.bitStringOfVertices
        while bitStringOfVertices > 0:
            lowestBit = bitStringOfVertices & -bitStringOfVertices
            indexOfU = lowestBit.bit_length() - 1
            self.indexOfSoftClause[indexOfU] = -1
            bitStringOfVertices ^= lowestBit
    
    def transformGraph(self):
        for softClauseIndex in range(len(self.softClauses)):
            if self.softClauses[softClauseIndex].isInProperSet:
                self.addNewVertexToSoftClause(softClauseIndex)

    def addNewVertexToSoftClause(self, softClauseIndex: int):
        u = self.maxNumber
        self.positionInVertices[u] = self.numberVertices
        self.adjacency.append(0)
        self.indexOfSoftClause[self.numberVertices] = softClauseIndex
        softClause = self.softClauses[softClauseIndex]
        # for v in softClause.vertices:
        #     self.addHardClause(u, v)
        # softClause.vertices.append(u)
        bitStringOfVertices = softClause.bitStringOfVertices
        while bitStringOfVertices > 0:
            lowestBit = bitStringOfVertices & -bitStringOfVertices
            indexOfV = lowestBit.bit_length() - 1
            self.adjacency[indexOfV] |= (1 << self.numberVertices)
            self.adjacency[self.numberVertices] |= (1 << indexOfV)
            bitStringOfVertices ^= lowestBit
        self.softClauses[softClauseIndex].bitStringOfVertices |= (1 << self.numberVertices)
        self.numberVertices += 1
        self.maxNumber += 1

def filtcol(bitStringSubGraph: GraphBitString, numberColors: int, lastVertexOfColorHasVertex: list, verticesOfNewSubGraph: list, isInNewPrunedSet: list, isInNewBranchingSet: list):
    for vertex in verticesOfNewSubGraph:
        if not isInNewBranchingSet[vertex]:
            continue
        if not bitStringSubGraph.filt(numberColors - 1, vertex, lastVertexOfColorHasVertex, verticesOfNewSubGraph, isInNewPrunedSet):
            return False
    return True

def filtsat(bitStringSubGraph: GraphBitString, numberVertices: int, numberColors: int, colorOfVertices: list, verticesOfNewSubGraph: list, isInNewPrunedSet: list, isInNewBranchingSet: list):
    subGraphSAT = GraphSAT(bitStringSubGraph, verticesOfNewSubGraph)
    independentSet = [0 for _ in range(numberColors)]
    for u in verticesOfNewSubGraph:
        if isInNewBranchingSet[u] and isInNewPrunedSet[u]:
            print(f"Warning filtsat: Vertex {u} is in both branching and pruned sets.")
        if isInNewPrunedSet[u]:
            # if colorOfVertices[u] < 0:
            #     # print(f"Error filtsat: Vertex {u} has no color assigned.")
            #     continue
            independentSet[colorOfVertices[u]] |= (1 << subGraphSAT.positionInVertices[u])
        elif isInNewBranchingSet[u]:
            independentSet[numberColors - 1] |= (1 << subGraphSAT.positionInVertices[u])
    for color in range(numberColors):
        # if not independentSet[color]:
        #     print(f"Warning filtsat: No independent set found for color {color}.")
        #     continue
        subGraphSAT.addSoftClause(independentSet[color])
    for u in verticesOfNewSubGraph:
        if not isInNewPrunedSet[u] and not isInNewBranchingSet[u]:
            continue
        if not subGraphSAT.unitPropagation(u):
            indexOfU = subGraphSAT.positionInVertices[u]
            softClause = subGraphSAT.softClauses[subGraphSAT.indexOfSoftClause[indexOfU]]
            if softClause.bitStringOfVertices.bit_count() == 1:
                return False
            subGraphSAT.removeVertex(u)
            if isInNewPrunedSet[u]:
                isInNewPrunedSet[u] = False
            else:
                isInNewBranchingSet[u] = False
    return True

def satcol(bitStringSubGraph: GraphBitString, numberColors: int, colorOfVertices: list, verticesOfNewSubGraph: list, isInNewPrunedSet: list, isInNewBranchingSet: list):
    global checkPrint
    subGraphSAT = GraphSAT(bitStringSubGraph, verticesOfNewSubGraph)
    independentSet = [0 for _ in range(numberColors)]
    bitStringOfListOfVertices = 0
    for u in verticesOfNewSubGraph:
        if isInNewPrunedSet[u]:
            independentSet[colorOfVertices[u]] |= (1 << subGraphSAT.positionInVertices[u])
        else:
            bitStringOfListOfVertices |= (1 << bitStringSubGraph.positionInVertices[u])
    for color in range(numberColors - 1):
        subGraphSAT.addSoftClause(independentSet[color])
    # if not listOfVertices:
    #     print("Warning satcol: No vertices left to process.")
    # listOfVertices.sort()
    # if checkPrint:
    #     print(f"numberColors: {numberColors}")
    #     print(f"verticesOfNewSubGraph: {verticesOfNewSubGraph}")
    #     print(f"colorOfVertices:")
    #     for color in range(numberColors):
    #         print(f"Color {color}: {independentSet[color]}", end=' ')
    #         print()
    #     print(f"\nlistOfVertices: {listOfVertices}")
    while bitStringOfListOfVertices > 0:
        independentSet = []
        bitStringOfIndependentSet = 0
        copyBitStringOfListOfVertices = bitStringOfListOfVertices
        while copyBitStringOfListOfVertices > 0:
            lowestBit = copyBitStringOfListOfVertices & -copyBitStringOfListOfVertices
            indexOfU = lowestBit.bit_length() - 1
            u = bitStringSubGraph.vertices[indexOfU]
            independentSet.append(u)
            bitStringOfIndependentSet |= (1 << subGraphSAT.positionInVertices[u])
            bitStringOfListOfVertices ^= (1 << indexOfU)
            copyBitStringOfListOfVertices &= bitStringSubGraph.adjacency[indexOfU]
        subGraphSAT.addSoftClause(bitStringOfIndependentSet)
        ok = False
        for u in independentSet:
            if subGraphSAT.unitPropagation(u):
                ok = True
                break
        # if checkPrint:
        #     print(f"Processing independent set: {independentSet}, ok: {ok}")
        if ok:
            subGraphSAT.removeLastSoftClause()
        else:
            for u in independentSet:
                isInNewPrunedSet[u] = True
                isInNewBranchingSet[u] = False
            subGraphSAT.transformGraph()
            # if checkPrint:
            #     for softClause in subGraphSAT.softClauses:
            #         print(f"Soft clause: {softClause.vertices}, state: {softClause.state}, isInProperSet: {softClause.isInProperSet}")

lowerBound = 0
currentSolution = []
bestSolution = []

def findMaxClique(graph: Graph, currentPrunedSet: list, currentBranchingSet: list, currentUpperBound: list):
    global lowerBound, bestSolution, currentSolution, checkPrint
    if checkPrint:
        return
    # if len(currentSolution) == 1 and currentSolution[0] == 191:
    #     print(f"findMaxClique: Current solution size: {len(currentSolution)}, Lower bound: {lowerBound} :")
    #     print(currentSolution)
    #     print(f"Current pruned set: {currentPrunedSet}")
    #     print(f"Current branching set: {currentBranchingSet}")
    #     checkPrint = True
    # else:
    #     checkPrint = False
    numberVerticesOfGraph = graph.numberVertices
    newUpperBound = currentUpperBound[:]
    isInCurrentPrunedSet = [False] * numberVerticesOfGraph
    bitStringOfCurrentSubGraph = 0
    for u in currentPrunedSet:
        isInCurrentPrunedSet[u] = True
        bitStringOfCurrentSubGraph |= (1 << u)
    for u in currentBranchingSet:
        bitStringOfCurrentSubGraph |= (1 << u)
    bitStringSubGraph = GraphBitString(numberVerticesOfGraph, graph.connected, currentPrunedSet + currentBranchingSet)
    for u in currentBranchingSet:
        copyAdjacencyOfU = graph.listVerticesOfNeighbors(graph.adjacency[u] & bitStringOfCurrentSubGraph)
        verticesOfNewSubGraph = []
        for v in copyAdjacencyOfU:
            if isInCurrentPrunedSet[v] or v < u:
                verticesOfNewSubGraph.append(v)
        newUpperBound[u] = 1
        for v in verticesOfNewSubGraph:
            newUpperBound[u] = max(newUpperBound[u], newUpperBound[v] + 1)
        if newUpperBound[u] + len(currentSolution) <= lowerBound:
            isInCurrentPrunedSet[u] = True
        else:
            # if len(currentSolution) == 3 and currentSolution[0] == 166 and currentSolution[1] == 111 and currentSolution[2] == 26 and u == 7:
            #     print(f"Processing vertex {u} in findMaxClique.")
            #     checkPrint = True
            # else:
            #     checkPrint = False
            # for i in range(1, len(verticesOfNewSubGraph)):
            #     if verticesOfNewSubGraph[i] < verticesOfNewSubGraph[i - 1]:
            #         print("Error findMaxClique: Vertices of new subgraph are not sorted correctly.")
            if not verticesOfNewSubGraph:
                if len(currentSolution) >= len(bestSolution):
                    bestSolution = currentSolution[:]
                    bestSolution.append(u)
                    lowerBound = max(lowerBound, len(bestSolution))
                newUpperBound[u] = min(newUpperBound[u], lowerBound - len(currentSolution))
                continue
            numberColors = lowerBound - len(currentSolution)
            # if numberColors < 0:
            #     print("Error findMaxClique: Number of colors is negative.")
            isInNewPrunedSet = [False] * numberVerticesOfGraph
            isInNewBranchingSet = [True] * numberVerticesOfGraph
            if numberColors > 1:
                check, colorOfVertices = bitStringSubGraph.iseq(numberColors - 1, verticesOfNewSubGraph)
                # if check != copyCheck or copyListOfVertices != listOfVertices or colorOfVertices != copyColorOfVertices:
                #     checkPrint = True
                #     bitStringSubGraph = GraphBitString(numberVerticesOfGraph, graph.connected, currentPrunedSet + currentBranchingSet)
                #     check1, colorOfVertices1, listOfVertices1 = bitStringSubGraph.iseq(numberColors - 1, copy1ListOfVertices)
                #     check2, colorOfVertices2 = bitStringSubGraph.iseq_old(numberColors - 1, copy2ListOfVertices)
                #     return
                # else:
                #     checkPrint = False
                if not check:
                    newUpperBound[u] = min(newUpperBound[u], lowerBound - len(currentSolution))
                    continue
                isIndependentSet = True
                bitStringOfBranchingSet = 0
                for v in verticesOfNewSubGraph:
                    if colorOfVertices[v] >= 0:
                        isInNewPrunedSet[v] = True
                        isInNewBranchingSet[v] = False
                    else:
                        indexOfV = bitStringSubGraph.positionInVertices[v]
                        if bitStringSubGraph.adjacency[indexOfV] & bitStringOfBranchingSet != bitStringOfBranchingSet:
                            isIndependentSet = False
                        bitStringOfBranchingSet |= (1 << indexOfV)
                if isIndependentSet:
                    lastVertexOfColor = [-1] * numberColors
                    for v in verticesOfNewSubGraph:
                        if isInNewPrunedSet[v]:
                            if colorOfVertices[v] < 0:
                                print(f"Error findMaxClique: Vertex {v} has no color assigned.")
                            lastVertexOfColor[colorOfVertices[v]] = v
                    lastVertexOfColorHasVertex = [-1] * numberVerticesOfGraph
                    for v in verticesOfNewSubGraph:
                        if isInNewPrunedSet[v]:
                            lastVertexOfColorHasVertex[v] = lastVertexOfColor[colorOfVertices[v]]
                    if not filtcol(bitStringSubGraph, numberColors, lastVertexOfColorHasVertex, verticesOfNewSubGraph, isInNewPrunedSet, isInNewBranchingSet) or \
                       not filtsat(bitStringSubGraph, numberVerticesOfGraph, numberColors, colorOfVertices, verticesOfNewSubGraph, isInNewPrunedSet ,isInNewBranchingSet):
                        newUpperBound[u] = min(newUpperBound[u], lowerBound - len(currentSolution))
                        continue
                    # for v in newPrunedSet:
                    #     lastVertexOfColor[colorOfVertices[v]] = v
                    # lastVertexOfColorHasVertex = [-1] * bitStringSubGraph.numberVertices
                    # for v in newPrunedSet:
                    #     lastVertexOfColorHasVertex[v] = lastVertexOfColor[colorOfVertices[v]]
                    # if not filtcol(bitStringSubGraph, numberColors, lastVertexOfColorHasVertex, newPrunedSet, newBranchingSet) or \
                    #     not filtsat(bitStringSubGraph, numberVerticesOfGraph, numberColors, colorOfVertices, newPrunedSet, newBranchingSet):
                    #     newUpperBound[u] = min(newUpperBound[u], lowerBound - len(currentSolution))
                    #     continue
                    # if len(newPrunedSet) + len(newBranchingSet) < cnt:
                    #     print(f"Yeah, filtered out some vertices: {cnt} -> {len(newPrunedSet) + len(newBranchingSet)}")
                else:
                    satcol(bitStringSubGraph, numberColors, colorOfVertices, verticesOfNewSubGraph, isInNewPrunedSet, isInNewBranchingSet)
            # if len(currentSolution) == 1 and currentSolution[0] == 166 and u == 111:
            #     print(f"New pruned set: {len(newPrunedSet)}")
            #     print(f"Vertices of new subgraph: {len(verticesOfNewSubGraph)}")
            #     print(f"Number Colors: {numberColors}")
            newPrunedSet = [v for v in verticesOfNewSubGraph if isInNewPrunedSet[v]]
            newBranchingSet = [v for v in verticesOfNewSubGraph if isInNewBranchingSet[v]]
            if newBranchingSet:
                currentSolution.append(u)
                findMaxClique(graph, newPrunedSet, newBranchingSet, newUpperBound)
                currentSolution.pop()
            elif len(currentSolution) >= len(bestSolution):
                bestSolution = currentSolution[:]
                bestSolution.append(u)
                lowerBound = max(lowerBound, len(bestSolution))
        newUpperBound[u] = min(newUpperBound[u], lowerBound - len(currentSolution))

a = []

# def inputFile(filename):
#     global a
#     with open(filename, 'r') as f:
#         [numberVertices, numberEdges] = [int(x) for x in f.readline().split()]
#         vertices = [Vertex(i) for i in range(numberVertices)]
#         edges = []
#         adjacency = [0 for _ in range(numberVertices)]
#         for _ in range(numberEdges):
#             line = f.readline().split()
#             u, v = int(line[1]), int(line[2])
#             u -= 1
#             v -= 1
#             edges.append((u, v))
#             adjacency[u] |= (1 << v)
#             adjacency[v] |= (1 << u)
#         a = [int(x) for x in f.readline().split()]
#         return numberVertices, vertices, edges, adjacency

# numberVertices, vertices, edges, adjacency = inputFile('src/.inp')

[numberVertices, numberEdges] = [int(x) for x in readline().split()]
vertices = [Vertex(i) for i in range(numberVertices)]
edges = []
adjacency = [0 for _ in range(numberVertices)]
for _ in range(numberEdges):
    line = readline().split()
    u, v = int(line[0]), int(line[1])
    u -= 1
    v -= 1
    edges.append((u, v))
    adjacency[u] |= (1 << v)
    adjacency[v] |= (1 << u)
# a = [int(x) for x in readline().split()]

startTime = time.time()
graph = Graph(numberVertices, vertices, adjacency)
degSort(graph, edges)
verticesOfGraphForGreedy = [VertexForGreedyAlgorithm(vertex.positionInListSorted) for vertex in graph.vertices]
graphForGreedy = GraphForGreedyAlgorithm(numberVertices, verticesOfGraphForGreedy, edges)
initialSolution = graphForGreedy.findClique()
# initialSolution = [0]
lowerBound = len(initialSolution) 
initialUpperBound = buildInitialUpperBound(graph)
bestSolution = initialSolution[:]
# timeLimit = 60 * 60
for i in range(lowerBound, graph.numberVertices):
    # print(f"Processing subgraph starting with vertex {i}...")
    verticesOfSubGraph = []
    adjacencyOfI = graph.adjacency[i]
    while adjacencyOfI > 0:
        lowestBit = adjacencyOfI & -adjacencyOfI
        j = lowestBit.bit_length() - 1
        if j >= i:
            break
        verticesOfSubGraph.append(j)
        adjacencyOfI ^= lowestBit
    if (len(verticesOfSubGraph) < lowerBound):
        continue
    prunedSet = [verticesOfSubGraph[j] for j in range(lowerBound)]
    branchingSet = verticesOfSubGraph[lowerBound:]
    currentSolution.append(i)
    findMaxClique(graph, prunedSet, branchingSet, initialUpperBound)
    currentTime = time.time()
    # if currentTime - startTime > timeLimit:
    #     break
    currentSolution.pop()
    initialUpperBound[i] = lowerBound
endTime = time.time()
totalTime = endTime - startTime

print(len(bestSolution))
for vertex in bestSolution:
    print(graph.vertices[vertex].id + 1, end=' ')
print(f"Total Time: {totalTime} seconds")

# with open('src/.out', 'w') as f:
#     for i in range(numberVertices):
#         print(graph.vertices[i].id + 1, end=' ', file=f)
#     print(file=f)
#     # print(f"{numberVertices} {len(edges)}")
#     # for edge in edges:
#     #     print(edge[0] + 1, edge[1] + 1, file=f)
#     # for vertex in graph.vertices:
#     #     if vertex.id + 1 in a:
#     #         print(f"Vertex {vertex.id + 1} (position in sorted list: {vertex.positionInListSorted})", file=f)
#     print(f"initialSolution: {len(initialSolution)}", file=f)
#     for vertex in initialSolution:
#         print(graph.vertices[vertex].id + 1, end=' ', file=f)
#     print(f"\n", file=f)
#     print(len(bestSolution), file=f)
#     # print solution
#     for vertex in bestSolution:
#         # print(graph.vertices[vertex].id + 1, end=' ', file=f)
#         print(vertex, end=' ', file=f)
#     print(f"\nCheck Solution: {graph.checkSolution(bestSolution)}", file=f)
#     print(f"Total Time: {totalTime} seconds", file=f)