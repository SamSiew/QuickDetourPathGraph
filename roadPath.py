"""
Name: Siew Ming Shern
StudentID: 28098552
"""
class Edge:
    def __init__(self,u=None,v=None,w=None):
        self.u = u
        self.v = v
        self.w = w
        self.isToll = False

class Vertex:
    def __init__(self,value=None):
        self.data = value
        self.originPath = None
        self.discovered = False
        self.visited = False
        self.time = None
        self.edge = []
        self.isRedLight = False
        self.isService = False

    def reset(self):
        """
        reset value of instance variable

        Time complexity: 0(1)
        Space complexity: 0(1)
        Error handle: Nothing
        :return: Nothing
        Parameter: Nothing
        Precondition: Nothing
        """
        self.originPath = None
        self.discovered = False
        self.visited = False
        self.time= None

    def update(self,time,origin):
        """
        update vertex with new time and
        origin

        Time complexity: 0(1), update vertex needs O(1) time to set value
        Space complexity: 0(1), update vertex needs O(1) space to set value
        Error handle: Nothing
        :return: Nothing
        :param time: new time of this vertex from another vertex
        :param origin: vertexID
        Precondition: origin must be a valid vertexID and
        time must also be positive value
        """
        self.time= time
        self.originPath = origin
        self.discovered = True

    def copy(self,value):
        """
        Functionality of the function
        Time complexity: O(E), occurs when E time is needed to store edges from original vertex
        Space complexity: O(E), occurs when E space is needed to store edges from original vertex
        where E is the total number of roads/edges in the road network.
        Error handle:
        :param value: total number of vertex in graph
        :return: object reference to a copy of this vertex
        Precondition: value must be positive value
        and valid total number of vertex in graph
        """
        newcopy = Vertex(self.data + value)
        for edges in self.edge:
            newcopy.edge.append(Edge(edges.u + value, edges.v + value,edges.w))
        return newcopy

class Graph:
    def __init__(self):
        self.currentGraph = []
        self.detourGraph = []

    def __len__(self):
        return len(self.currentGraph)

    def __getitem__(self, item):
        return self.currentGraph[item]

    def buildGraph(self, filename_roads):
        """
        Functionality of the function

        Time complexity: O(V + E), occurs when V times is needed to create array with v size and
        E time is needed to store list of edges in filename and it is also needed to store edges in graph.
        Thus, this functions should takes O(V + E).

        Space complexity: O(V + E), occurs when V space is needed to create array with v size and
        all vertices will holds atmost E edges in total which will takes E space

        where V is the total number of points/ nodes/ vertices in the road network
        and E is the total number of roads/edges in the road network.

        Error handle: when error involving read and open file occurs, return and do nothing

        :param filename_roads: filename of content which contains list of edges
        :return: Nothing
        """
        try:
            file = open(filename_roads, 'r')

            edgeList = []

            for line in file:
                edgeList.append(line.strip().split())

            file.close()
        except:
            return
        #let maximum vertex number be 0
        maxV = 0
        #find the maximum vertexID in edge list
        for index in range(len(edgeList)):
            #instantiate the list item to edge object
            edgeList[index] = Edge(int(edgeList[index][0]), int(edgeList[index][1]), float(edgeList[index][2]))
            #if current max is less than current edge, u :
            #let u be max vertex
            if maxV < edgeList[index].u:
               maxV = edgeList[index].u
            # if current max is less than current edge, v :
            # let v be max vertex
            if maxV < edgeList[index].v:
                maxV = edgeList[index].v

        #create a new graph with atmost maximum vertexID
        self.currentGraph = [None] * (maxV + 1)
        # new graph will need to reset the detour graph
        self.detourGraph = []
        #initialise all vertex with index start 0...maximum vertexID
        for index in range(len(self.currentGraph)):
            self.currentGraph[index] = Vertex(index)
        #add the all edge to currentGraph correspond to edge.u
        for edges in edgeList:
            self.currentGraph[edges.u].edge.append(edges)

    def quickestPath(self,source, target):
        """
        Time complexity: O(E Log V), occurs for each Edges, rise method of heap
        is needed to add item into heap and  perform to ensure new item fits
        in the heap will takes O(Log v) which result O(E)*O(Log V)

        Space complexity: O(V), occurs when V space is needed to store V vertices in path
        and additional V space is needed to store atmost all vertex in graph
        into min heap which result O(v).

        where V is the total number of points/ nodes/ vertices in the road network
        and E is the total number of roads/edges in the road network.

        Error handle: nothing

        :param source: vertexID of start path
        :param target: vertexID of destination path
        :return: list[] contains list of path from source to target and time of that path

        Precondition: sources and target are vertex that in graph and graph is a non empty, simple graph.

        """
        #make sure source is integer
        source = int(source)
        # make sure target is integer
        target = int(target)
        #if new graph is empty or source out of range
        #or target out of range then return [[],-1]
        if len(self.currentGraph) == 0 or source < 0 or source > len(self) - 1 or  target < 0 or target > len(self) - 1 :
            return [[], -1]
        #reset all vertex in currentGraph
        for vertex in self.currentGraph:
            vertex.reset()
        #set source vertices time as 0
        self[source].time= 0
        # set origin path of source as target, such that when
        # backtracking take place, it shall end here.
        self[source].originPath = target
        #let discover be minheap holds at most all vertices
        discover = MinHeap(len(self))
        #add source to discover
        discover.add(0, source)
        #runs while discover has vertex
        while discover.hasVertexLeft():
            #get current minimum in discover
            item = discover.extract()
            #get vertex of discovered
            u = self[item[0]]
            #set u vertex has discovered
            u.visited = True
            #for each edges, choose to either add all vertex of current vertex or
            #update the distance from source to some target
            for edges in u.edge:
                v = self[edges.v]
                #when adjacent vertex is not visited yet
                if v.visited == False:
                    # when adjacent vertex is not discovered,
                    # add to discover
                    if v.discovered == False:
                        discover.add(item[1] + edges.w, edges.v)
                        v.update(item[1] + edges.w, item[0])
                    else:
                        # if current time taken is better than the vertex time taken
                        # and it is already discovered, update the time in discover
                        if v.time > item[1] + edges.w:
                            discover.update(edges.v, item[1] + edges.w)
                            v.update(item[1] + edges.w, item[0])
        #when origin path for target is none, surely, there
        #isnt a path from source, then return [[], -1]
        if self[target].originPath == None:
            return [[], -1]
        #get origin vertex of the target
        location = self[target].originPath
        #use for storing all path in reverse order
        path = [target]
        #while location is not at target yet
        while location != target:
            #add location to path
            #change location to origin vertex of a vertex
            path.append(location)
            location = self[location].originPath
        #let first pointer be first k item pointer
        fPointer = 0
        # let last pointer be last k item pointer
        lPointer = len(path) - 1
        #while firsr pointer has not reach last pointer
        #swap item on first pointer with item on last pointer at every iteration
        while fPointer <= lPointer:
            path[fPointer], path[lPointer] = path[lPointer], path[fPointer]
            fPointer += 1
            lPointer -= 1
        #to be return: path, time taken of path
        retVal = (path, self[target].time)
        return retVal

    def augmentGraph(self,filename_camera, filename_toll):
        """
        Functionality of the function

        Time complexity: O(V + E), occurs when V times is needed to create array with V size which are
        list of vertices considered unsafe,and list of edges with E size in graph considered unsafe. Furthermore,
        V time atmost to marks the unsafe vertex in graph and additional E time atmost to search for edges
        and  marks unsafe edges which results O(V + E)

        Space complexity: O(V + E), occurs when V space is needed to create array with V size
        which are considered red light and E space is needed for all edges which are toll

        where V is the total number of points/ nodes/ vertices in the road network
        and E is the total number of roads/edges in the road network.

        Error handle: when error involving read and open file occurs, return and do nothing

        :param filename_camera: filename contains list of vertices considered unsafe
        :param filename_toll: filename contains list of edges which are considered having toll and are unsafe
        :return: nothing
        """
        try:
            file = open(filename_camera, 'r')

            redLightList = []

            for line in file:
                redLightList.append(int(line.strip()))

            file.close()

            file = open(filename_toll, 'r')

            tollList = []

            for line in file:
                item = line.strip().split()
                item[0] = int(item[0])
                item[1] = int(item[1])
                tollList.append(item)

            file.close()
        except:
            return
        #reset all marks as redlight and toll on vertex to false
        for vertex in self.currentGraph:
            vertex.isRedLight = False
            for edges in vertex.edge:
                edges.isToll = False
        #marks the vertex as redlight for each redlight vertexID
        for redLight in redLightList:
            self.currentGraph[redLight].isRedLight = True
        #marks the edges as toll for each toll edges.
        for vertex in tollList:
            for edges in self.currentGraph[vertex[0]].edge:
                if edges.v == vertex[1]:
                    edges.isToll = True
                    break

    def quickestSafePath(self, source, target):
        """
        Time complexity: O(E Log V), occurs when there is no redlight and no toll which causes
        for each Edges which will be considered safe, rise method of heap is needed to perform
        to ensure new item fits in the heap which result O(E)*O(Log V)

        Space complexity: O(V), occurs when there is no redlight and no toll
        V space is needed to store V vertices in path and additional V space is needed to store
        atmost all vertex in graph into min heap which result O(v).

        where V is the total number of points/ nodes/ vertices in the road network
        and E is the total number of roads/edges in the road network.

        Error handle: nothing

        :param source: vertexID of start path
        :param target: vertexID of destination path
        :return: list[] contains list of path from source to target and time of that path

        Precondition: sources and target are vertex that in graph and graph is a non empty, simple graph.
        """
        # make sure source is integer
        source = int(source)
        # make sure target is integer
        target = int(target)
        # if new graph is empty or source out of range or target out of range then return [[],-1]
        if len(self.currentGraph) == 0 or source < 0 or source > len(self) - 1 or target < 0 or target > len(self) - 1:
            return [[], -1]
        # reset all vertex in currentGraph
        for vertex in self.currentGraph:
            vertex.reset()
        # set source vertices time as 0
        self[source].time = 0
        #if source is not redlight, continue
        #else return [[], -1]
        if self[source].isRedLight == False:
            # set origin path of source as target, such that when
            # backtracking take place, it shall end here.
            self[source].originPath = target
        else:
            return [[], -1]
        # let discover be minheap holds at most all vertices
        discover = MinHeap(len(self))
        # add source to discover
        discover.add(0, source)
        # runs while discover has vertex
        while discover.hasVertexLeft():
            # get current minimum in discover
            item = discover.extract()
            # get vertex of discovered
            u = self[item[0]]
            # set u vertex has discovered
            u.visited = True
            #when u is not a redLight, perform dijakstra algorithm
            if u.isRedLight == False:
                # for each edges, choose to either add all vertex of current vertex or
                # update the distance from source to some target
                for edges in u.edge:
                    v = self[edges.v]
                    #make sure v is not visited, v is not redlight and edges is not a toll
                    if v.visited == False and v.isRedLight == False and edges.isToll == False:
                        # when adjacent vertex is not visited yet
                        if v.discovered == False:
                            # when adjacent vertex is not discovered,
                            # add to discover
                            discover.add(item[1] + edges.w, edges.v)
                            v.update(item[1] + edges.w,item[0])
                        else:
                            # if current time taken is better than the vertex time taken
                            # and it is already discovered, update the time in discover
                            if v.time> item[1] + edges.w:
                                discover.update(edges.v, item[1] + edges.w)
                                v.update(item[1] + edges.w, item[0])

        # when origin path for target is none, surely, there
        # isnt a path from source, then return [[], -1]
        if self[target].originPath == None:
            return [[], -1]
        # get origin vertex of the target
        location = self[target].originPath
        # use for storing all path in reverse order
        path = [target]
        # while location is not at target yet
        while location != target:
            # add location to path
            # change location to origin vertex of a vertex
            path.append(location)
            location = self[location].originPath
        # let first pointer be first k item pointer
        fPointer = 0
        # let last pointer be last k item pointer
        lPointer = len(path) - 1
        # while firsr pointer has not reach last pointer
        # swap item on first pointer with item on last pointer at every iteration
        while fPointer <= lPointer:
            path[fPointer], path[lPointer] = path[lPointer], path[fPointer]
            fPointer += 1
            lPointer -= 1
        # to be return: path, time taken of path
        retVal = (path, self[target].time)
        return retVal

    def addService(self, filename_service):
        """
        Functionality of the function

        Time complexity: O(V + E), occurs when list of service on vertices
        are as much as V size is needed to be created, a new graph is created based
        on original graph created in basicGraph.txt with number of vertices and
        number of edges is doubled and because the original graph will have
        two copy of its graph on new graph with one of a copy graph have vertices which are marked
        as an service node when it is a service. At this rate, the new graph,c(g , g')
        have a graph, g which have a service node in vertices marks with another graph, g'.
        Both are exactly the same with index differ only.
        Finally, an additional atmost E edges is needed to be created to link a vertices in g'
        and vertices in g which is a services and is a vertices mirror to a vertices g'. This leads
        to O(2V + 3E) which is also equivalent to O(V + E).

        Space complexity: O(V + E), occurs when V space is needed to create array with 2V size
        and all vertices will holds atmost 3E edges in total which will takes O(V + E)

        where V is the total number of points/ nodes/ vertices in the road network
        and E is the total number of roads/edges in the road network.

        Error handle: when error involving read and open file occurs, return and do nothing

        :param filename_service: file contain vertices which are service
        :return: Nothing
        """
        try:
            file = open(filename_service, 'r')
            serviceList = []
            for line in file:
                serviceList.append(int(line.strip()))
            file.close()
        except:
            return
        #reset detour graph
        self.detourGraph = []
        #get numberofroad from length of actual graph
        numOfRoad = len(self.currentGraph)
        #add all adjencylist form current graph into detourGraph
        for item in self.currentGraph:
            self.detourGraph.append(item)
        #add a copy of current Graph into detour graph again but each vertex ID is increment to length of current graph
        for index in range(len(self.currentGraph)):
            items = self.currentGraph[index].copy(len(self))
            self.detourGraph.append(items)
        #marks a set of vertices in one of graph in detourgraph for being a service node
        for service in serviceList:
            self.detourGraph[service].isService = True
        #add an edges from between g' to g where edges from g' to g and g is a service node
        for index in range(len(self.currentGraph)):
            for edges in self.detourGraph[index].edge:
                if self.detourGraph[edges.v].isService:
                    self.detourGraph[index + numOfRoad].edge.append(edges)


    def quickestDetourPath(self, source, target):
        """
        Time complexity: O(E Log V), occurs for each Edges, rise method of heap is needed to perform
        to ensure new item fits in the heap. However, when additional graph is created,
        number of vertex is doubled  and number of edges is tripled atmost to a new graph which
        would lead to O(3E)*O(Log 2V) = O(3E)*O(Log V + log 2) which result O( E Log V).

        Space complexity: O(V), occurs when V space is needed to store V vertices in path and
        2V space is needed to store atmost all vertex in new graph into min heap
        which result O(V)

        where V is the total number of points/ nodes/ vertices in the road network
        and E is the total number of roads/edges in the road network.

        Error handle: nothing

        :param source: vertexID of start path
        :param target: vertexID of destination path
        :return: list[] contains list of path from source to target and time of that path

        Precondition: sources and target are vertex that in graph and graph is a non empty, simple graph.
        """
        # make sure source is integer
        source = int(source)
        # make sure target is integer
        target = int(target)
        # if new graph is empty or source out of range or target out of range then return [[],-1]
        if len(self.currentGraph) == 0 or source < 0 or source > len(self) - 1 or target < 0 or target > len(self) - 1:
            return [[], -1]
        # reset all vertex in currentGraph
        for vertex in self.detourGraph:
            vertex.reset()
        #when source is already a service node, perform dijakstra at original
        #copy of graph, g
        if self.detourGraph[source].isService:
            startPoint = source
        # when source is not a service node, perform dijakstra at
        # another copy of graph, g'
        else:
            startPoint = source + len(self.currentGraph)
        # set source vertices time as 0
        self.detourGraph[startPoint].time= 0
        # set origin path of source as target, such that when
        #backtracking take place, it shall end here.
        self.detourGraph[startPoint].originPath = target
        # let discover be minheap holds at most all vertices
        discover = MinHeap(len(self.detourGraph))
        # add source to discover
        discover.add(0, startPoint)
        # runs while discover has vertex
        while discover.hasVertexLeft():
            # get current minimum in discover
            item = discover.extract()
            # get vertex of discovered
            u = self.detourGraph[item[0]]
            # set u vertex has discovered
            u.visited = True
            # for each edges, choose to either add all vertex of current vertex or
            # update the distance from source to some target
            for edges in u.edge:
                v = self.detourGraph[edges.v]
                # when adjacent vertex is not visited yet
                if v.visited == False:
                    # when adjacent vertex is not discovered,
                    # add to discover
                    if v.discovered == False:
                        discover.add(item[1] + edges.w, edges.v)
                        v.update(item[1] + edges.w,item[0])
                    else:
                        # if current time taken is better than the vertex time taken
                        # and it is already discovered, update the time in discover
                        if v.time> item[1] + edges.w:
                            discover.update(edges.v, item[1] + edges.w)
                            v.update(item[1] + edges.w,item[0])
        # when origin path for target is none, surely, there
        # isnt a path from source, then return [[], -1]
        if self.detourGraph[target].originPath == None:
            return [[], -1]
        # get origin vertex of the target
        location = self.detourGraph[target].originPath
        # use for storing all path in reverse order
        path = [target]
        # while location is not at target yet
        while location != target:
            # add location to path
            # change location to origin vertex of a vertex
            path.append(location)
            location = self.detourGraph[location].originPath
        # let first pointer be first k item pointer
        fPointer = 0
        # let last pointer be last k item pointer
        lPointer = len(path) - 1
        # while first pointer has not reach last pointer
        # swap item on first pointer with item on last pointer at every iteration
        while fPointer <= lPointer:
            #ensure first k item does not goes out of range by changing its value back
            #to a vertex mirror to its
            if path[fPointer] >= len(self): path[fPointer] -= len(self)
            # ensure last k item does not goes out of range by changing its value back
            # to a vertex mirror to its
            if path[lPointer] >= len(self): path[lPointer] -= len(self)

            path[fPointer], path[lPointer] = path[lPointer], path[fPointer]

            fPointer += 1

            lPointer -= 1
        # to be return: path, time taken of path
        retVal = (path, self.detourGraph[target].time)

        return retVal

#   Actual codes for class is extracted from FIT1008, my_heap.py on week 12
#   https://lms.monash.edu/course/view.php?id=42395&section=19#19
#   Used by 28098552, Siew Ming Shern for assignment 1.
class MinHeap:

    def __init__(self, count):
        """
        :param count: number of item to be inserted to this heap
        Time Complexity:
            Best Case:  O(k)
            Worst Case: O(k)
            where k is kth value specify by user.
            Best Case == Worst Case, this function will always create a list for N item .
        Space Complexity:
            Best Case:  O(k)
            Worst Case: O(k)
            where N is Maximum number of item that can be inserted to this class.
            Best Case == Worst Case, this function will always create a list for N item which require 4 bytes each N items.
        """
        self.count = 0
        self.map = [None] * (count)
        self.array = [None] * (count + 1)

    def hasVertexLeft(self):
        """
        :return: boolean value, true if number of vertex in graph is not empty
                                false if number of vertex in graph is empty
        Time Complexity: O(1), O(1) time is needed to find if length of this heap is zero
        Space Complexity: O(1), O(1) space is needed to find if length of this heap is zero and
                                return boolean value
        """
        return len(self) > 0

    def __len__(self):
        """
        :return: self.count
        Time Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
            Best Case == Worst Case, this function is returning only number of item in This Heap excluding the first element,
            so that computation can be convienience.
        Space Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
            Best Case == Worst Case, this function will always return a byte for any number.
        """
        return self.count

    def swap(self, i, j):
        """
        :param i: index of an item in array
        :param j: index of an item in array
        :return: no return, it is void method.
        Time Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
            Best Case == Worst Case, this function is swapping position of two item in array.
        Space Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
            Best Case == Worst Case, no memory required for this operation since it does not create any memory for swapping.
        """

        self.map[self.array[i][0]], self.map[self.array[j][0]] = self.map[self.array[j][0]], self.map[self.array[i][0]]
        self.array[i], self.array[j] = self.array[j], self.array[i]

    def rise(self, k):
        """
        :param k: index of an item, mainly which is index of new item added to heap.
        :return: nothing to return, it is just swapping method for new item to fit in this heap.
        Time Complexity:
            Best Case:  O(1)
            Worst Case: O(log k)
            where k is number of non-empty item in list.

            Best Case: this function can terminate early if new child nodes is already bigger than parent nodes
            and if new child nodes is equal to parent nodes but already have bigger index value that it holds.
            Therefore, no swapping is required at all.
            Worst Case: when new child nodes is the smallest item of the heaps or the new child nodes is equal to all parent nodes
            but have a bigger index. Therfore, any time the new child nodes encounter any above cases, swapping is required at log n times
        Space Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
            Best Case == Worst Case, no memory required for this operation since it does not create any memory for swapping.
        """
        #rise item to parent nodes if parents is already bigger than current new nodes
        while k > 1 and self.array[self.getparent(k)][1] > self.array[k][1]:
            self.swap(k, self.getparent(k))
            k = self.getparent(k)
        # swapping parent with new nodes when it equals but have higher index value.
        while k > 1 and self.array[self.getparent(k)][1] == self.array[k][1] and self.array[self.getparent(k)][0] > self.array[k][0]:
            self.swap(k, self.getparent(k))
            k = self.getparent(k)

    def add(self, wordFreq, counter):
        """
        :param wordFreq: a value defines word count to be inserted in to array
        :param counter: index of that item in original list
        :return: nothing, just a void method
        Time Complexity:
            Best Case:  O(1)
            Worst Case: O(log k)
            where k is number of non-empty item in list.
            Best Case: this function can terminate early if new child nodes is already bigger than parent nodes
                and if new child nodes is equal to parent nodes but already have bigger index value that it holds.
                Therefore, no swapping is required at all.
            Worst Case: when new child nodes is the smallest item of the heaps or the new child nodes is equal to all parent nodes
                but have a bigger index. Therfore, any time the new child nodes encounter any above cases, swapping is reuired at log n times
        Space Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
            Best Case == Worst Case, no memory required for this operation since it does not create any memory for swapping.
        """
        newitem = [counter, wordFreq]
        #add item when current class length is less than heap array length
        if self.count + 1 < len(self.array):
            self.array[self.count + 1] = newitem
            self.count += 1
            self.map[newitem[0]] = self.count
            self.rise(self.count)

    def update(self,v,time):
        """
        :param v: vertexID
        :param time: new value for time
        :return: nothing
        Time Complexity: O(log v), when new child nodes is the smallest
        item of the heaps or the new child nodes is equal to all parent nodes
        but have a bigger index. Therefore, any time the new child nodes
        encounter any above cases, rising is required at log v times.

        Space Complexity: O(1), no memory required for this operation
        since it does not create any memory for swapping.

        where v is number of vertex in heap.
        """
        index = self.map[v]
        self.array[index][1] = time
        self.rise(index)

    def getmin(self):
        """
        :return: get minimum value of frequency of an array which always located at index 1.
        Time Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
        Best Case == Worst Case, access item of array which only tooks O(1).
        Space Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
        Best Case == Worst Case, access item of array which only tooks O(1). Therefore, a byte memory is all its needs.
        """
        return self.array[1][1]

    def extract(self):
        """
        :return:
        Time Complexity:
            Best Case:  O(log k)
            Worst Case: O(log k)
            where k is number of non-empty item in list.
            BestCase == Worst Case: Happens when new roots node is needed to swap all the way to other item to become a leaf node.
        Space Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
        Best Case == Worst Case, This operation is just swapping operation(in-place function). Therefore, no memory is required to
        perform it.
        """
        #store item to new variable, swap last child nodes with root then reduce length of array
        item = self.array[1]
        self.swap(1, self.count)
        self.array[self.count] = None
        self.count -= 1
        #sink method is called to adjust array to heap structure.
        self.sink(1)
        return item

    def sink(self, k):
        """
        :param k: index of an item in array
        :return: nothing, just a void method
        Time Complexity:
            Best Case:  O(log k)
            Worst Case: O(log k)
            where k is number of non-empty item in list.
            When last item is swapped with root, it always require comparing with child node and needs to sink all the way to become leaf node
            because in heap structure, child node in scallability is always bigger than any parent nodes of older root node
        Space Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
        Best Case == Worst Case, This operation is just swapping operation(in-place function) therefore, no memory is required to
        perform it.
        """
        #if array has left item, loops continute
        while self.getleft(k) <= self.count:
            child = self.getsmallestchild(k)
            #when parent nodes is already smaller than child node, break
            if self.array[k][1] < self.array[child][1]:
                break
            #when parent nodes is already equal to child node, then index value with higher is place at root
            elif self.array[k][1] == self.array[child][1] and self.array[k][0] < self.array[child][0]:
                break
            self.swap(child, k)
            k = child

    def getparent(self,k):
        """
        :param k: refers to index of item in array
        :return: computation of parent of k index in array.
        Time Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
        Best Case == Worst Case, this function just performing arithmetic operation.
        Space Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
        Best Case == Worst Case, for this function just performing arithmetic operation. Therefore, a byte memory is all its needs
        """
        return k//2

    def getleft(self,k):
        """
        :param k: refers to index of item in array
        :return: computation to get left child nodes
        Time Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
        Best Case == Worst Case, this function just performing arithmetic operation.
        Space Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
        Best Case == Worst Case, for this function just performing arithmetic operation. Therefore, a byte memory is all its needs
        """
        return 2*k

    def getright(self,k):
        """
        :param k: refers to index of item in array
        :return: computation to get right child nodes
        Time Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
        Best Case == Worst Case, this function just performing arithmetic operation.
        Space Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
        Best Case == Worst Case, for this function just performing arithmetic operation. Therefore, a byte memory is all its needs
        """
        return 2*k + 1

    def getsmallestchild(self, k):
        """
        :param k: index of item in array
        :return: either left index of k item or right item of k index
        Time Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
            Best Case == Worst Case, this function just comparing which child nodes is smaller when last nodes is swap with minimum
            after extract() method is called.
        Space Complexity:
            Best Case:  O(1)
            Worst Case: O(1)
            Best Case == Worst Case, This function only involve storing index into variable and accessing item in array.
            Therefore, two byte memory is all its needs resulting O(1)
        """
        left = self.getleft(k)
        right = self.getright(k)
        #check if there is length of class is equal to root with left item only or if left child is smaller, return left if true
        if left == len(self) or self.array[left][1] < self.array[right][1]:
            return left
        #check if there left child node is equal to right node and index value of left is bigger, return left if true
        elif self.array[left][1] == self.array[right][1] and self.array[left][0] < self.array[right][0]:
            return left
        #return right when we knew right node is assume smaller.
        else:
            return right


if __name__ == '__main__':
    # User Input
    filename = 'basicGraph.txt'
    cameraFile = 'camera.txt'
    tollFile = 'toll.txt'
    serviceFile = 'servicePoint.txt'
    startVertex = 4
    endVertex = 3
    # Driver Operation
    adjencyList = Graph()
    adjencyList.buildGraph(filename)
    adjencyList.augmentGraph(cameraFile, tollFile)
    adjencyList.addService(serviceFile)
    task1Result = adjencyList.quickestPath(startVertex, endVertex)
    task2Result = adjencyList.quickestSafePath(startVertex, endVertex)
    task3Result = adjencyList.quickestDetourPath(startVertex, endVertex)
    source = str(startVertex)
    sink = str(endVertex)
    print("---------------------------------------------------------------------")
    print("Enter the file name for the graph: " + filename)
    print("Enter the file name for camera nodes: " + cameraFile)
    print("Enter the file name for the toll roads: " + tollFile)
    print("Enter the file name for the service nodes: " + serviceFile)
    print("---------------------------------------------------------------------")
    print("Source node: " + source)
    print("Sink node: " + sink)
    print("---------------------------------------------------------------------")
    print("Quickest path:")
    if task1Result[1] == -1:
        print("No path exists")
        print("Time: 0 minute(s)")
    else:
        path = str(task1Result[0][0])
        for i in range(1, len(task1Result[0])):
            path += " --> " + str(task1Result[0][i])
        print(path)
        print("Time: " + str(task1Result[1]) + " minute(s)")
    print("---------------------------------------------------------------------")
    print("Safe quickest path:")
    if task2Result[1] == -1:
        print("No path exists")
        print("Time: 0 minute(s)")
    else:
        path = str(task2Result[0][0])
        for i in range(1, len(task2Result[0])):
            path += " --> " + str(task2Result[0][i])
        print(path)
        print("Time: " + str(task2Result[1]) + " minute(s)")
    print("---------------------------------------------------------------------")
    print("Quickest detour path:")
    if task3Result[1] == -1:
        print("No path exists")
        print("Time: 0 minute(s)")
    else:
        path = str(task3Result[0][0])
        for i in range(1, len(task3Result[0])):
            path += " --> " + str(task3Result[0][i])
        print(path)
        print("Time: " + str(task3Result[1]) + " minute(s)")
    print("---------------------------------------------------------------------")
    print("Program end")
