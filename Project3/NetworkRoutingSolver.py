#!/usr/bin/python3
import Proj3GUI
from CS312Graph import *
import time
import math


class NetworkRoutingSolver:
    def __init__(self):
        pass

    def initializeNetwork(self, network):
        assert (type(network) == CS312Graph)
        self.network = network

    def getShortestPath(self, destIndex):
        self.dest = destIndex
        # TODO: RETURN THE SHORTEST PATH FOR destIndex
        #       INSTEAD OF THE DUMMY SET OF EDGES BELOW
        #       IT'S JUST AN EXAMPLE OF THE FORMAT YOU'LL 
        #       NEED TO USE

        path_edges = []
        total_length = 0

        sourceNode = self.network.nodes[self.source]
        sourceID = sourceNode.node_id
        destNode = self.network.nodes[destIndex]
        currentID = destNode.node_id
        nodes = self.network.nodes

        if self.prev[currentID] == -math.inf:
            return {'cost': total_length, 'path': []}

        while currentID != sourceID:
            prevID = self.prev[currentID]
            if prevID is None:
                return {'cost': total_length, 'path': []}
            nodeDist = self.dist[currentID]
            path_edges.append((nodes[currentID].loc, nodes[prevID].loc, '{:.0f}'.format(nodeDist)))
            total_length += nodeDist
            currentID = prevID

        return {'cost': total_length, 'path': path_edges}

    def computeShortestPaths(self, srcIndex, destIndex, use_heap=False):
        self.source = srcIndex
        t1 = time.time()
        # TODO: RUN DIJKSTRA'S TO DETERMINE SHORTEST PATHS.
        #       ALSO, STORE THE RESULTS FOR THE SUBSEQUENT
        #       CALL TO getShortestPath(dest_index)
        nodes = self.network.nodes

        if not use_heap:
            dist, prev = self.computePathsArray(srcIndex, nodes, destIndex)
        else:
            dist, prev = self.computePathsHeap(srcIndex, nodes, destIndex)

        self.dist = dist
        self.prev = prev
        t2 = time.time()
        return (t2 - t1)

    # Time complexity: O((E + V)log(|V|)) because we call deleteMinHeap() |V| times worst case and we also call
    # decreaseKeyHeap() |E| times worst case. Both of these functions are O(log|V|)

    # Space complexity: O(|V|) would be O(4|V|) but simplifies
    def computePathsHeap(self, srcIndex, nodes, destIndex):
        # create dist and previous arrays
        dist = [math.inf] * len(nodes)
        prev = [-math.inf] * len(nodes)
        dist[srcIndex] = 0
        singleDist = [math.inf] * len(nodes)
        self.heapCount = len(nodes) - 1
        self.heap, self.pointer = self.makeHeap(nodes, srcIndex)
        # start while loop- run until there are no nodes left
        while self.heapCount > 0:
            index, weight = self.deleteMinHeap()
            # if index == destIndex:
            #     break
                # explore current nodes neighbors
            for edge in nodes[index].neighbors:
                currIndex = edge.dest.node_id
                currWeight = weight + edge.length
                # if current distance in dist[] array is larger- we have found a shorter path
                if dist[currIndex] > currWeight:
                    dist[currIndex] = currWeight
                    singleDist[currIndex] = edge.length
                    prev[currIndex] = index
                    # call decrease key
                    self.decreaseKeyHeap(currIndex, currWeight)
        return singleDist, prev

    # Time complexity: O(|V|)- loop over array of size |V| once to create heap
    # Space complexity: O(|V|)- creating heep and pointer arrays which are of size |V|
    def makeHeap(self, nodes, srcIndex):
        pointer = [-1] * len(nodes)
        heap = [] * len(nodes)
        # set pointer and heap for the source node
        pointer[srcIndex] = 0
        heap.append((srcIndex, 0))
        counter = 1
        for node in nodes:
            # append all nodes except for the source node
            if node.node_id != srcIndex:
                heap.append((node.node_id, math.inf))
                pointer[node.node_id] = counter
                counter += 1
        return heap, pointer

    # Time complexity: O(log|V|) because we are doing a few order one operations then calling siftDown() which is O(log|V|)
    # Space complexity: O(1)
    def deleteMinHeap(self):
        # get index and weight of root node
        index, weight = self.heap[0]
        self.pointer[index] = - 1
        # set root node to value of last node in array
        self.heap[0] = self.heap[self.heapCount]
        newIndex, newWeight = self.heap[0]
        self.pointer[newIndex] = 0
        self.heapCount -= 1
        # call siftDown()
        self.siftDown()
        return index, weight


    # Time complexity: O(log|V|) because we are doing two order one operations and calling bubbleUp which has time
        # complexity of O(log|V|)
    # Space complexity: O(1) not allocating much memory
    def decreaseKeyHeap(self, index, dist):
        # get heap index from pointers to update the node value
        updateIndex = self.pointer[index]
        self.heap[updateIndex] = (index, dist)
        self.bubbleUp(updateIndex)

    # Time complexity: O(log|V|) just traversing a tree which is logarithmic at each layer of tree
    # Space complexity: O(1)
    def siftDown(self):
        # takes root node and see if it is larger than it's children
        isBigger = False
        count = 0
        currNode = self.heap[0]
        currVal = currNode[1]
        isRightNone = False
        while not isBigger:
            # get left and right children
            leftChild, rightChild = self.getChild(count)
            if rightChild is None:
                isRightNone = True
            if leftChild is not None:
                # if current value and it's children are infinity we are done
                if currVal == math.inf and not isRightNone:
                    if leftChild[1] == math.inf and rightChild[1] == math.inf:
                        break
                # check left child
                if not isRightNone and currVal > leftChild[1] and leftChild[1] < rightChild[1]:
                    parent = self.heap[count]
                    self.heap[count] = leftChild
                    self.heap[self.pointer[leftChild[0]]] = parent
                    self.pointer[parent[0]] = self.pointer[leftChild[0]]
                    self.pointer[leftChild[0]] = count
                    count = self.pointer[parent[0]]

                # check right child
                elif not isRightNone and currVal > rightChild[1]:
                    parent = self.heap[count]
                    self.heap[count] = rightChild
                    self.heap[self.pointer[rightChild[0]]] = parent
                    self.pointer[parent[0]] = self.pointer[rightChild[0]]
                    self.pointer[rightChild[0]] = count
                    count = self.pointer[parent[0]]
                else:
                    isBigger = True
            else:
                isBigger = True
        return

    # Time Complexity: O(1)- doing a few order order one operations
    # Space Complexity: O(1)- not allocating very much memory
    def getChild(self, index):
        # formula for index of left and right children
        leftIndex = (index * 2) + 1
        rightIndex = (index * 2) + 2
        leftChild = None
        rightChild = None
        # make sure we stay in bounds of our heap
        if leftIndex < self.heapCount:
            leftChild = self.heap[leftIndex]
        if rightIndex < self.heapCount:
            rightChild = self.heap[rightIndex]
        return leftChild, rightChild

    # Time complexity: O(log|V|) because we are bubbling up a tree which is log time at each level
    # Space complexity: O(1) not really allocating any memory
    def bubbleUp(self, index):
        isGreater = False
        childIndex = index
        while not isGreater:
            # if we are at root break
            if childIndex == 0:
                break
            parentIndex = self.getParent(childIndex)
            child = self.heap[childIndex]
            parent = self.heap[parentIndex]
            # swap if child is less than parent
            if child[1] < parent[1]:
                self.heap[parentIndex] = child
                self.heap[childIndex] = parent
                # update the pointer array
                self.pointer[child[0]] = parentIndex
                self.pointer[parent[0]] = childIndex
                childIndex = parentIndex
            else:
                isGreater = True
        return

    # Time complexity: O(1)
    # Space complexity: O(1)
    def getParent(self, index):
        return math.floor((index - 1) / 2)

    # Time complexity: O(|V^2|)
    # Space complexity: O(|V|) simplifies from O(4|V|)
    def computePathsArray(self, srcIndex, nodes, destIndex):
        # Initialize array to keep track of distances
        dist = [math.inf] * len(nodes)
        singleDist = [math.inf] * len(nodes)
        singleDist[srcIndex] = 0
        dist[srcIndex] = 0
        self.pq_size = len(nodes)
        # initialize priority queue to be empty
        pq = []
        prev = [-math.inf] * len(nodes)
        pq.append((srcIndex, 0))
        # run while loop until size of our priority queue is zero
        while self.pq_size > 0:
            # call deletemin
            index, weight, pq_index = self.deletemin(pq)
            # if index == destIndex:
            #     break
            # set distance of current node in our priority queue to -1 so we don't use it again
            pq[pq_index] = (index, -1)
            # visit neighbors of current node
            for edge in nodes[index].neighbors:
                curr_index = edge.dest.node_id
                curr_weight = edge.length + weight
                # if node hasn't been visited or we found shorter path update distance array
                if (dist[curr_index] == math.inf) or (dist[curr_index] > curr_weight):
                    singleDist[curr_index] = edge.length
                    # this is essentially the decrease key for the array implementation
                    pq.append((curr_index, curr_weight))
                    dist[curr_index] = curr_weight
                    prev[curr_index] = index
        return singleDist, prev

    # FOR NUMBER TWO Decrease key is O(1) because we are appending it to the pq
    # Space complexity is also O(1)

    # Time complexity: O(|V|) because we have to search through array
    # Space complexity: O(1) because only looking through array
    def deletemin(self, pq):
        pq_index = 0
        pq_min_index = 0
        min_weight = 0
        min_index = 0
        # search in priority queue for the next smallest value to remove from the queue
        for (x, y) in pq:
            if min_weight == 0 and y != -1 or y < min_weight and y != -1:
                min_weight = y
                min_index = x
                pq_min_index = pq_index
            pq_index = pq_index + 1
        self.pq_size -= 1
        return min_index, min_weight, pq_min_index
