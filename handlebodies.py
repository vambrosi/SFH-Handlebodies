#-----------------------------------------------------------------------------#
# Sutured Handlebodies (given by a graph)
#-----------------------------------------------------------------------------#

from multimodules import MultiModule
from multimodules_constructors import (
    AlgebraHomology, CFDD_id, CFAA_id, PartialDehnTwist, SolidPairOfPants)
from box_tensor_product import box_tensor

import math
import csv
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------#
# Handlebody classes
#-----------------------------------------------------------------------------#

class Genus2Handlebody:
    def __init__(self, edge_labels):
        # TO DO: Still doesn't work when the edge_labels give a train track
        # that does not satisfy the triangle inequality.

        # Check that the geometric intersection number of each compressing disk
        # with the multicurve given by the sutures is even.
        
        if not all(points % 2 == 0 for points, _ in edge_labels):
            raise Exception('Number of intervals must be even.')

        # If the twist parameters are odd then the roles of R_- and R_+ in the
        # 2nd ball must be reversed for the gluing to work. The parameter
        # self.dual keeps track of that. Implicity in this description is the
        # fact that the twist parameters must have same have the same parity.
        # The next few lines also check if that condition is satisfied.

        if all(twist % 2 == 1 for _, twist in edge_labels):
            self.dual = True
        elif all(twist % 2 == 0 for _, twist in edge_labels):
            self.dual = False
        else:
            raise Exception('Invalid twist parameters.')

        self.edge_labels = edge_labels

        # Does not compute SFH on initialization. Only when you call the
        # function self.SFH() or when you call self.SFH_rank or self.SFH_ranks
        self.SFH_rank = None
        self.SFH_ranks = None            

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

    #-------------------------------------------------------------------------#
    # Computing SFH
    #-------------------------------------------------------------------------#

    def SFH(self, load=True, save=False, check=False):
        '''Computes the rank of SFH for all spin_c structures.'''
        
        # If check=True it loads precomputed values to compare with the
        # new computed values. For debugging purposes only.
        if check:
            try:
                self.SFH_ranks_from_csv()
                save = False
        
            except FileNotFoundError:    
                raise Exception('Precomputed values not found.')

        # Loads ranks if possible.
        elif load:
            try: 
                return self.SFH_ranks_from_csv() 
            except FileNotFoundError:
                print('Computing ranks. (Might take a while.)')
                pass

        # edge_labels gives the number of arcs in the arc diagram of each disk
        # and how much we have to rotate the arc diagram before gluing.
        
        edge_labels = [(points//2 - 1, (twist//2) % (points//2))
                       for points, twist in self.edge_labels]

        # Numbers of arcs for each parametrized disk in the 1st and 2nd ball
        arcs1 = [edge_labels[0][0], edge_labels[1][0], edge_labels[2][0]]
        arcs2 = [edge_labels[0][0], edge_labels[2][0], edge_labels[1][0]]

        # DDDModule for the 1st and 2nd ball, respectively.
        S1 = SolidPairOfPants(arcs1[0], arcs1[1], arcs1[2])
        
        if self.dual:
            S2 = SolidPairOfPants(arcs1[0], arcs1[1], arcs1[2]).dual()
        else:
            S2 = SolidPairOfPants(arcs2[0], arcs2[1], arcs2[2])

        AA_id = {} # making sure it doesn't compute twice.
        for n in arcs1:
            if not n in AA_id:
                AA_id[n] = CFAA_id(n)

        for n in arcs1:
            S1 = box_tensor(S1, 0, AA_id[n], 0)

        if self.dual:
            smallest_twist = {}
            
            for n in arcs1:
                if not n in smallest_twist:
                    smallest_twist[n] = AlgebraHomology(n).dual()
                S2 = box_tensor(S2, 0, smallest_twist[n], 1)

        else:
            for n in arcs2:
                S2 = box_tensor(S2, 0, AA_id[n], 0)
       
        DehnTwists = {} # making sure it doesn't compute twice again.
        for arcs, twist in edge_labels:
            if not (arcs, twist) in DehnTwists: 
                DehnTwists[(arcs, twist)] \
                = PartialDehnTwist(arcs, twist, AA_id[arcs])

        for arcs, twist in edge_labels:
            S1 = box_tensor(S1, 0, DehnTwists[(arcs, twist)], 0)

        # Glue two SolidPairOfPants
        M = box_tensor(S1, 0, S2, 0)
        if self.dual:
            N = M.HH(0, 2)
        else:
            N = M.HH(0, 3)
        N = N.HH(0, 1)
              
        # Creates a dictionary that associates the spin_c grading with the
        # rank of the corresponding sutured Floer homology. The spin_c grading
        # is given by the number of strands in the idempotents associated with
        # each generator. We only look at the idempotents corresponding to
        # disks 2 and 3 because they generate the homology group H_2(M).
        SFH_ranks = {}
        for gen in N.generators:
            point = M.count_chords(gen)[0:2]
            try:
                SFH_ranks[point] += 1
            except KeyError:
                SFH_ranks[point] = 1

        if check:
            if (self.SFH_rank == len(N.generators)
                    and self.SFH_ranks == SFH_ranks):
                print(f'SFH agrees with loaded file for {self.edge_labels}')
            else:
                print(f'SFH != from loaded file for {self.edge_labels}')
        else:
            self.SFH_rank = len(N.generators)
            self.SFH_ranks = SFH_ranks

        if save: self.SFH_ranks_to_csv()
        return self.SFH_rank, self.SFH_ranks

    #-------------------------------------------------------------------------#
    # Loading/Saving computations
    #-------------------------------------------------------------------------#

    def SFH_ranks_from_csv(self): # fix
        el = self.edge_labels
        name = f'{el[0][0]}-{el[0][1]} ' \
               + f'{el[1][0]}-{el[1][1]} ' \
               + f'{el[2][0]}-{el[2][1]}'

        with open(f'./SFHGenus2/{name}.csv', newline='') as csvfile:
            reader = csv.reader(csvfile, quoting=csv.QUOTE_MINIMAL)

            SFH_ranks = {}
            for row in reader:
                SFH_ranks[(int(row[0]), int(row[1]))] = int(row[2])

        self.SFH_ranks = SFH_ranks
        
        SFH_rank = 0
        for _, rank in SFH_ranks.items():
            SFH_rank += rank

        self.SFH_rank = SFH_rank

        return SFH_rank, SFH_ranks

    def SFH_ranks_to_csv(self):
        el = self.edge_labels
        name = f'{el[0][0]}-{el[0][1]} ' \
               + f'{el[1][0]}-{el[1][1]} ' \
               + f'{el[2][0]}-{el[2][1]}'

        with open(f'./SFHGenus2/{name}.csv', 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)

            for point, rank in self.SFH_ranks.items():
                writer.writerow([point[0], point[1], rank])

    #-------------------------------------------------------------------------#
    # Plotting SFH ranks
    #-------------------------------------------------------------------------#

    def SFH_plot(self, symmetric=True, **kwargs):
        points_by_rank = self.SFH_by_rank()

        el = [x//2 for x, _ in self.edge_labels]
        xlim = (-2, el[1] + 1)
        ylim = (-2, el[2] + 1)

        if symmetric:
            # Draw the plot in a triangular grid.
            fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(3,3), dpi=200)

            # Draw rank for points where SFH has nonzero rank
            s3 = math.sqrt(3)
            for rank in points_by_rank:
                ys = [-3*y - 1 for _, y in points_by_rank[rank]]
                xs = [s3*(2*x+y-1) for x, y in points_by_rank[rank]]
                
                latex_rank = f'${rank}$'
                m = max(points_by_rank)
                ax.scatter(xs, ys, marker='s', color='w', zorder=2)
                ax.scatter(xs, ys, marker=latex_rank, zorder=3, linewidth=0.1,
                           facecolor=plt.cm.inferno(math.sqrt((rank-1)/m)),
                           edgecolor='black')

            # Draw an triangular grid with diameter prescribed above.
            xs, ys = [], []
            for x in range(xlim[0] + 1, xlim[1]):
                for y in range(ylim[0] + 1, ylim[1]):
                    if ((el[1] + el[2] - el[0] - 4)/2 < x + y 
                        <= math.ceil((el[1] + el[2] + el[0] - 3)/2) + 1):
                        xs.append(s3 * (2*x+y-1))
                        ys.append(-3*y - 1)
            
            ax.scatter(xs, ys, marker='.', color='black',
                       s=0.5, alpha=0.5, zorder=1)

            xlim_tmp, ylim_tmp = xlim, ylim
            xlim = tuple(s3*(2*x+y-1) for x, y in zip(xlim_tmp, ylim_tmp))
            ylim = tuple(-3*y - 1 for y in ylim_tmp)

        else:
            # Draw the plot in a rectangular grid.
            fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(2,2), dpi=150)

            # Draw points_by_rank with nonzero rank
            for rank in points_by_rank:
                xs = [point[0] for point in points_by_rank[rank]]
                ys = [point[1] for point in points_by_rank[rank]]
                latex_rank = f'${rank}$'
                m = max(points_by_rank)
                ax.scatter(xs, ys, marker=latex_rank, zorder=3, linewidth=0.1,
                           facecolor=plt.cm.inferno(math.sqrt((rank-1)/m)),
                           edgecolor='black')

            # Draws a square grid with side given by sides.
            xs, ys = [], []
            for x in range(xlim[0]+1, xlim[1]):
                for y in range(ylim[0]+1, ylim[1]):
                    if not (x, y) in self.SFH_ranks:
                        xs.append(x)
                        ys.append(y)
            
            ax.scatter(xs, ys, marker='.', color='black', s=0.5, alpha=0.5)

        ax.set(**kwargs)
        ax.set_axis_off()
        ax.set(xlim=xlim, ylim=ylim)
        ax.set_aspect('equal')
        
        plt.show()
        plt.close()

    def SFH_by_rank(self):
        points_by_rank = {}
        for point, rank in self.SFH_ranks.items():
            try:
                points_by_rank[rank].append(point)
            except KeyError:
                points_by_rank[rank] = [point]
        return points_by_rank

class SuturedHandlebody: ## not working yet
    def __init__(self, graph_dict, vertices=None):
        if vertices == None:
            self.vertices = set()
        else:
            self.vertices = vertices
        
        self.edges = graph_dict

        for v in self.graph_dict:
            self.vertices.add(v)
            for w in self.graph[v]:
                self.vertices.add(w)

    def find_spanning_tree(self):
        tree = {}
        
        # pick any vertex and add to the tree
        v = next(iter(self.vertices))
        tree[v] = {}

        # start a stack of iterators

#-----------------------------------------------------------------------------#
# Test Functions
#-----------------------------------------------------------------------------#

def check_SFH_pretzel_knots():
    '''Consider the complement of the Seifert surface of the pretzel knot
    P(2r+1, 2s+1, 2t+1) given by the checkerboard coloring of the canonical
    projection. This function computes the rank of the SFH of the associated
    sutured manifold and it compares with the result in Example 8.4 of

    S. FRIEDL, A. JUHASZ, J. RASMUSSEN - 
    THE DECATEGORIFICATION OF SUTURED FLOER HOMOLOGY
    '''

    # Case where r,s,t > 0
    for r in range(1,4):
        for s in range(1,4):
            for t in range(1,4):
                M = Genus2Handlebody([(2*(t+r+1), 2*(t+r+1) - 1),
                                      (2*(r+s+1), 2*(r+s+1) - 1),
                                      (2*(s+t+1), 2*(s+t+1) - 1)])
                if M.SFH(load=False)[0] == (r+1)*(t+1) + (t+r+1)*s:
                    print(f'Works for r={r}, s={s}, t={t}')
                else:
                    print(f"Doesn't work for r={r}, s={s}, t={t}")
