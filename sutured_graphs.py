"""
Python implementation of sutured graphs. Those are labelled trivalent ribbon
graphs. There are two possible labels for each vertex (+ or -) and edges are
labelled by two integers: the first is nonnegative and even, and the second is
odd if incident vertex labels differ and even if they are the same.

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


from typing import Optional
from constructors import AAEdge, DDDVertex
from modules import MultiModule
from tensor import tensor


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

    def incident_to(self):
        '''Returns vertices incident to the Edge.'''
        return [self[0], self[1]]


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

        # Initialize Sutured Floer ranks
        self.SFH_rank = None
        self.SFH_ranks = None

    def __getitem__(self, name: str) -> Vertex:
        '''
        Returns a vertex with name given.
        '''
        for vertex in self.vertices:
            if vertex.name == name:
                return vertex

        raise IndexError('No vertex with that name.')

    def __call__(self, obj: str | Vertex | Edge, pos: int) -> \
            tuple[Vertex | Edge, int]:
        '''
        Given a vertex/edge end, returns the matched edge/vertex end.
        '''
        return self.incident_to(obj)[pos]

    #-------------------------------------------------------------------------#
    # Graph-editing methods
    #-------------------------------------------------------------------------#

    def add_vertex(self, name: str, sign: str) -> None:
        '''
        Adds a vertex with the indicated name and label to the graph.
        '''

        # Checks if there is a vertex with that name already
        try:
            self[name]
        except IndexError:
            v = Vertex(name, sign)
            self.vertices.add(v)
            self.incidence[v] = []
        else:
            raise IndexError('Graph already has a vertex with this name.')

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

        # Unpacking variables (if needed)
        if isinstance(end1, str):
            v = end1
            v_insertion = None
        else:
            v, v_insertion = end1
        if isinstance(end2, str):
            w_name = end2
            w_insertion = None
        else:
            w_name, w_insertion = end2

        sutures, twist = labels

        # Vertices and edge references
        v = self.get_vertex(v)
        w = self.get_vertex(w_name)
        e = Edge(v, w, sutures, twist)

        # In case insertion is not defined
        if v_insertion == None:
            v_insertion = len(self.incidence[v])
        if w_insertion == None:
            w_insertion = len(self.incidence[v])

        # Update graph
        self.incidence[v].insert(v_insertion, e)
        self.incidence[w].insert(w_insertion, e)
        self.edges.add(e)

    #-------------------------------------------------------------------------#
    # Local Properties
    #-------------------------------------------------------------------------#

    def incident_to(self, obj: str | Vertex | Edge) -> list[Edge | Vertex]:
        '''
        Returns the list of edges/vertices incident to the given vertex/edge.
        The order is given by the ribbon structure on the graph for vertices,
        and by the implied direction for edges.
        '''
        if isinstance(obj, Vertex | str):
            obj = self.get_vertex(obj)

            result = []
            for e in self.incidence[obj]:
                if e[0] == obj:
                    result.append((e, 0))
                elif e[1] == obj:
                    result.append((e, 1))

            return result

        elif isinstance(obj, Edge):
            result = []

            for v in obj.incident_to():
                i = self.incidence[v].index(obj)
                result.append((v, i))

            return result

    def get_vertex(self, v: Optional[str | Vertex]) -> Vertex:
        '''
        Function that normalizes vertex inputs. If the input is a vertex, it
        outputs it back, if it is a name it gives a vertex with that name back,
        if it exists, and gives a vertex in the graph otherwise.

        Mostly used internally to allow users to input names instead of
        vertices or to not make a vertex choice at all, when that is an option.
        '''

        if isinstance(v, Vertex):
            return v
        elif isinstance(v, str):
            return self[v]
        else:
            for v in self.vertices:
                return v

            # Only gets here if there are no vertices.
            raise IndexError('Graph is empty.')

    def valence(self, v: str | Vertex) -> int:
        '''
        Count the number of ends incident to the vertex. (It counts loops twice
        and edges with different ends once.) The function takes as input the
        name of the vertex or the vertex itself.
        '''

        return len(self.incident_to(v))

    def check_valence(self) -> bool:
        '''
        Checks if the graph is trivalent. This is necessary to compute SFH of
        the associated sutured handlebody.
        '''

        return all(self.valence(v) == 3 for v in self.vertices)

    def BS(self, a: str | Vertex | Edge) -> MultiModule:
        '''
        Returns the bordered sutured modules associated with a vertex or edge.
        '''

        if isinstance(a, str | Vertex):
            a = self.get_vertex(a)

            assert self.valence(a) == 3, 'Vertex must have valence 3.'
            e1, e2, e3 = self.incidence[a]

            M = DDDVertex(e1.sutures, e2.sutures, e3.sutures)
            return M if a.label == '+' else M.mirror()

        elif isinstance(a, Edge):
            return AAEdge(a.sutures, a.twist, a[0].label)

    #-------------------------------------------------------------------------#
    # Global Properties
    #-------------------------------------------------------------------------#

    def connected_components(self) -> set[set[Vertex]]:
        '''
        Returns the connected components as sets of vertices.
        '''
        components = set()

        # If one end of an edge is in a component, add both of its ends to it.
        # If the edge is not incident to any component, add a new component.
        for e in self.edges:
            for component in components:
                if e[0] in component or e[1] in component:
                    component.add(e[0])
                    component.add(e[1])
                    break
            else:
                components.add({e[0], e[1]})

        return components

    def number_components(self) -> int:
        return len(self.number_components())

    def is_connected(self) -> bool:
        return self.number_components() == 1

    def genus(self) -> int:
        '''
        Computes the genus of the graph.
        '''

        # Euler characteristic = number of components - genus.
        g = self.number_components() - len(self.vertices) + len(self.edges)

        return g

    def spanning_tree(self, root_name: str | Vertex = None) -> list[Edge]:
        '''
        Finds a spanning tree for the connected component of the vertex given.
        Input can be a vertex or a vertex name. It returns the set of edges in
        that spanning tree.
        '''

        # Get the root or pick a root if one was not choosen
        root = self.get_vertex(root_name)

        # Below we have a standard depth-first search.
        # We add every edge that does not create a loop.

        # Variable that stores the final result
        spanning_tree = []

        # Temporary variables used during the search.
        visited_vertices = {root}
        stack = [(e, e(root)) for e in self.incident_to(root)]

        while stack:
            e, v = stack.pop()

            if not v in visited_vertices:
                visited_vertices.add(v)
                spanning_tree.append(e)
                stack += [(e, e(v)) for e in self.incident_to(v)]

        return spanning_tree

    def fundamental_group_gens(self, root_name: str | Vertex = None) \
            -> list[Edge]:
        '''
        Chooses a spanning tree and returns the edges complementary to it.
        If the graph is connected those edges represent the generators of the
        fundamental group. If the graph is not connected, an error is raised.
        If needed, root_name lets you specify the root of the spanning tree.
        '''

        spanning_tree = self.spanning_tree(root_name=root_name)
        gens = self.edges.difference(set(spanning_tree))

        assert len(gens) == self.genus(), "Graph is not connected."
        return list(gens)

    #-------------------------------------------------------------------------#
    # Sutured Floer Invariants
    #-------------------------------------------------------------------------#

    @property
    def SFH_rank(self):
        '''Total rank of SFH of the Genus2Handlebody.'''
        return self.SFH()[0] if self._SFH_rank == None else self._SFH_rank

    @property
    def SFH_ranks(self):
        '''A dictionary that has as keys the spin_c gradings with nontrivial
        sutured Floer homology. The value associated with a key is the rank
        of SFH on the corresponding grading (with Z/2 coefficients).'''
        return self.SFH()[1] if self._SFH_ranks == None else self._SFH_ranks

    @SFH_rank.setter
    def SFH_rank(self, value):
        self._SFH_rank = value

    @SFH_ranks.setter
    def SFH_ranks(self, value):
        self._SFH_ranks = value

    def SFH(self, root_name: str | Vertex = None) -> MultiModule:
        '''
        Returns a dictionary {Spin^c grading: rank of SFH}.
        Assumes the sutured graph is connected.
        '''
        # Checks that graph is trivalent
        assert self.check_valence(), 'Graph is not trivalent.'

        # Get the root or pick a root if one was not choosen
        root = self.get_vertex(root_name)

        # Computes largest bordered invariant before having to take HHs.
        # In other words, it traverses a "graph" G such that all the vertices
        # and edges of self are included but some edges are not attached to
        # both of their ends. The result is required to be contractible. This
        # graph G has an associated bordered invariant that is computed in
        # the process (it is the M_left).

        # Initialize variables
        not_visited = self.vertices.union(self.edges)

        not_visited.remove(root)
        M_left = self.BS(root)
        left_slots = self.incident_to(root)

        while not_visited:
            # Pop the first vertex or edge that is not visited
            # Guarantees that we are not adding loops
            for i, (obj, pos) in enumerate(left_slots):
                if obj in not_visited:
                    left_index = i
                    left_obj, left_pos = obj, pos
                    left_slots.pop(i)
                    break

            M_right = self.BS(left_obj)
            right_slots = self.incident_to(left_obj)
            not_visited.remove(left_obj)

            # Pop the vertex or edge that is matched to left_slot
            for j, (obj, pos) in enumerate(right_slots):
                if self.incident_to(obj)[pos] == (left_obj, left_pos):
                    right_index = j
                    right_slots.pop(j)
                    break

            # Computes the bordered invariant by gluing left_slot
            # to the corresponing slot on the right.
            M_left = tensor(M_left, left_index, M_right, right_index)
            left_slots = left_slots + right_slots

        N = M_left

        # Identifies matched pair and apply HH to them.
        number_of_pairs = len(left_slots) // 2
        while left_slots:
            slot2 = left_slots.pop()
            slot2_index = len(left_slots)

            slot1_index = left_slots.index(self(*slot2))
            left_slots.pop(slot1_index)

            N = N.HH(slot1_index, slot2_index)

        SFH_ranks = {}
        for gen in N.generators:
            point = M_left.count_chords(gen)[number_of_pairs::]
            try:
                SFH_ranks[point] += 1
            except KeyError:
                SFH_ranks[point] = 1

        self.SFH_ranks = SFH_ranks
        self.SFH_rank = len(N.generators)

        return self.SFH_rank, self.SFH_ranks
