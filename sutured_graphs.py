"""
Python implementation of sutured graphs. Those are labelled trivalent ribbon
graphs. There are two possible labels for each vertex (+ or -) and edges are
labelled by two integers: the first is nonnegative and even, and the second is
odd if vertex labels differ and even if they are the same.

Those graphs encode a sutured handlebody. The trivalent ribbon graph encodes a
handlebody (by thickening it) together with a fixed parametrization of its
boundary (given by the ribbon structure). The labels on edges give the
Dehn-Thurston coordinates of the sutures, and the labels on vertices determine
which regions will be part of R_+ and which will be part of R_-.

Sutured graphs can be created using the SuturedGraph class (don't pass any
arguments, the vertices and edges should be added later). Then you use the
methods add_vertex or add_edge to build the graph.
"""

#-----------------------------------------------------------------------------#
#  Basic Classes
#-----------------------------------------------------------------------------#


class CyclicList(list):
    def __getitem__(self, index):
        if isinstance(index, int):
            return tuple.__getitem__(self, index % len(self))
        elif isinstance(index, slice):
            return tuple.__getitem__(self, index)

    def __setitem__(self, index):
        tuple.__setitem__(self, index % len(self))

    def succ(self, x):
        return self[(self.index(x)+1) % len(self)]

    def pred(self, x):
        return self[(self.index(x)-1) % len(self)]

    def insert(self, index, element):
        return tuple.insert(index % len(self), element)


class Edge(tuple):
    def __new__(cls, vertex1, vertex2, sutures, twist):
        # Checks if the number of intersections between sutures and the
        # associated compressing disk is even. This is necessary if we want to
        # split the complement of the sutures into two regions (R_+ and R_-).
        assert sutures % 2 == 0, "Number of sutures should be even."

        # Checks if twist satisfies the necessary properties
        if vertex1.label == vertex2.label:
            assert twist % 2 == 0, "Edge connecting vertices with same " + \
                "labels should have even twist."
        else:
            assert twist % 2 == 1, "Edge connecting vertices with " + \
                "different labels should have odd twist."

        return tuple.__new__(cls, (vertex1, vertex2, sutures, twist))

    def __eq__(self, other):
        # A pair of vertices can be incident to multiple edges with the same
        # parameters. Thus, to differentiate then we use id().
        return self is other

    def __ne__(self, other):
        return not (self is other)

    def __hash__(self):
        return id(self)

    def __call__(self, v):
        # If Edge gets called with one of its endpoints returns the other.
        # Otherwise raises an exception.
        if v == self[0]:
            return self[1]
        elif v == self[1]:
            return self[0]
        else:
            raise ValueError('Vertex is not an endpoint of this edge.')

    def __getitem__(self, index):
        if index == 0 or index == 1:
            return tuple.__getitem__(self, index)
        else:
            raise IndexError('An Edge only has two vertices.')

    def __repr__(self):
        return f'{self[0]} <--{self.twist}/{self.sutures}--> {self[1]}'

    @property
    def sutures(self):
        return tuple.__getitem__(self, 2)

    @property
    def twist(self):
        return tuple.__getitem__(self, 3)


class Vertex(tuple):
    '''
    Vertex is a tuple (name, label). Name should be a string or an int, and
    label should be '+' or '-'.
    '''
    def __new__(cls, name, label):
        assert isinstance(name, int | str), "Name should be string or integer."
        assert label == '+' or label == '-', "Label should be '+' or '-'."
        return tuple.__new__(cls, (name, label))

    @property
    def name(self):
        return self[0]

    @property
    def label(self):
        return self[1]

    def __eq__(self, other):
        # In case there is two vertices with the same name and label.
        return self is other

    def __ne__(self, other):
        return not (self is other)

    def __hash__(self):
        return id(self)

    def __repr__(self):
        if self.label == '+':
            return f'{self.name}\u208A'
        elif self.label == '-':
            return f'{self.name}\u208B'

#-----------------------------------------------------------------------------#
#  Main Class
#-----------------------------------------------------------------------------#


class SuturedGraph:
    def __init__(self):
        # Sets of vertices and edges. Can be acessed directly, but should
        # be changed only using the relevant methods (to keep them consistent).
        self.vertices = set()
        self.edges = set()

        # This dictionary has vertices as keys, and cyclic lists of all
        # corresponding incident edges as values. The order of the edges in the
        # cyclic list is the order given by the ribbon structure.
        self.incidence = dict()

    def __getitem__(self, name):
        '''
        Returns a vertex with name given (or None if there no such vertex).
        '''
        for vertex in self.vertices:
            if vertex.name == name:
                return vertex

    def add_vertex(self, name, sign):
        '''
        Adds a vertex with the indicated name and label to the graph.
        '''

        # Checks if there is a vertex with that name already
        assert self[name] is None, 'There is a vertex with that name already.'

        v = Vertex(name, sign)
        self.vertices.add(v)
        self.incidence[v] = []

    def add_edge(self, end1, end2, labels):
        '''
        This function takes 3 tuples. First and second tuples are of the form
        (vertex_name, edge_index) where the edge_index gives the position in
        the ribbon order of the given vertex where the edge should be inserted.

        This index works modulo the number of edges around the vertex. So if
        there are two edges incident to a vertex, an index = 0 (mod 2) would
        insert the edge before/after the two edges (the order is cyclic) and it
        would insert the edge in between the other edges if index = 1 (mod 2).

        The last tuple is of the form (sutures, twist) where those two numbers
        are integers, the former tells how many times the sutures intersect a
        compressing disk transversal to the edge, and the latter tells how much
        the sutures twist around the edge. The former should be even. The
        latter should be odd if the vertices have opposite labels and even
        otherwise.
        '''

        # Unpacking variables
        v_name, v_insertion = end1
        w_name, w_insertion = end2
        sutures, twist = labels

        # Vertices and edge references
        v = self.check_vertex(v_name)
        w = self.check_vertex(w_name)
        e = Edge(v, w, sutures, twist)

        # Update graph
        self.incidence[v].insert(v_insertion, e)
        self.incidence[w].insert(w_insertion, e)
        self.edges.add(e)

    def incident_to(self, v):
        '''
        Returns a CyclicList of edges incident to the given vertex. The cyclic
        order is given by the ribbon structure on the graph.
        '''
        return self.incidence[v]

    def check_vertex(self, v):
        '''
        Checks if the input is a vertex, if it is not, assumes the input is a
        vertex name, and returns the vertex with that name. (This function is
        manly used internally, to allow the user to input a vertex or a name.)
        '''

        if isinstance(v, Vertex):
            return v
        else:
            return self[v]

    def spanning_tree(self, root_name):
        '''
        Finds a spanning tree for the connected component of the vertex given.
        Input can be a vertex or a vertex name. It returns the set of edges in
        that spanning tree.
        '''

        # Gets the vertex if given a vertex name
        root = self.check_vertex(root_name)

        # Check if root is a vertex in the graph.
        assert root in self.vertices, 'Root is not in the graph.'

        # Below we have a standard depth-first search.
        # We add every edge that does not create a loop.

        # Variable that stores the final result
        spanning_tree = set()

        # Temporary variables used during the search.
        visited_vertices = {root}
        stack = [(e, e(root)) for e in self.incident_to(root)]

        while stack:
            e, v = stack.pop()

            if not v in visited_vertices:
                visited_vertices.add(v)
                spanning_tree.add(e)
                stack += [(e, e(v)) for e in self.incident_to(v)]

        return spanning_tree

    def valence(self, v_name):
        '''
        Count the number of ends incident to the vertex. (It counts loops twice
        and edges with different ends once.) The function takes as input the
        name of the vertex or the vertex itself.
        '''

        # Gets a vertex if given a vertex name
        v = self.check_vertex(v_name)

        # Since there can be loops, we have to separate into two cases.
        loops = [e for e in self.incident_to(v) if e[0] == e[1]]
        open_edges = [e for e in self.incident_to(v) if e[0] != e[1]]

        return 2*len(loops) + len(open_edges)

