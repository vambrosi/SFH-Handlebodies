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
        self.edge_labels = edge_labels
        self.SFH_rank = None
        self.SFH_ranks = None

    @property
    def SFH_rank(self):
        return self.SFH()[0] if self._SFH_rank == None else self._SFH_rank

    @property
    def SFH_ranks(self):
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

    def SFH(self, load=True, save=False, method='HH'):
        if load:
            try: 
                return self.SFH_ranks_from_csv() 
            except FileNotFoundError:
                print('Computing ranks. (May take a while.)')
                pass

        # Creates one SolidPairOfPants module (the two trimodules are equal) 
        # and take the box tensor products with the identity to get an 
        # AAA module.

        arcs = [label - 1 for label, _ in self.edge_labels]
        S1 = SolidPairOfPants(arcs[0], arcs[1], arcs[2])

        AA_id = {}
        for n in arcs:
            if not n in AA_id: # making sure it doesn't compute twice.
                AA_id[n] = CFAA_id(n) 

        for n in arcs:
            S1 = box_tensor(S1, 0, AA_id[n], 0)
        
        S2 = S1

        # edge_labels gives how many times R_+ intersects the compressing
        # disks and how much it is twisted along the handles.
        
        DehnTwists = {} # making sure it doesn't compute twice again.
        for n, k in self.edge_labels:
            if not (n, k) in DehnTwists: 
                DehnTwists[(n,k)] = PartialDehnTwist(n-1, k, AA_id[n-1])

        # Glue Dehn twists to one of the SolidPairOfPants
        for n, k in self.edge_labels:
            S2 = box_tensor(S2, 0, DehnTwists[(n,k)], 0)

        # Glue two SolidPairOfPants
        if method == 'HH':
            M = box_tensor(S1, 0, S2, 0, rename=True)
            N = M.HH(0, 2)
            N = N.HH(0, 1)

        elif method == 'box':
            N = box_tensor(S1, [0,1,2], S2, [0,1,2])
        
        self.SFH_rank = len(N.generators)
        
        SFH_ranks = {}
        for gen in N.generators:
            point = M.count_chords(gen)[0:2]
            try:
                SFH_ranks[point] += 1
            except KeyError:
                SFH_ranks[point] = 1

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
        
        self.SFH_rank = 0
        for _, rank in self.SFH_ranks.items():
            self.SFH_rank += rank

        return self.SFH_rank, self.SFH_ranks

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

    def SFH_plot(self, symmetric=True, diameter=None, **kwargs):
        points_by_rank = self.SFH_by_rank()
        
        if not points_by_rank:
            print('Trivial SFH, empty plot.')
            return

        xmax = max((point[0] for point in self.SFH_ranks))
        xmin = min((point[0] for point in self.SFH_ranks))
        
        ymax = max((point[1] for point in self.SFH_ranks))
        ymin = min((point[1] for point in self.SFH_ranks))

        xymax = max((point[0]+point[1] for point in self.SFH_ranks))
        xymin = min((point[0]+point[1] for point in self.SFH_ranks))

        xmid = (xmax + xmin)/2
        ymid = (ymax + ymin)/2
        xymid = (xymax + xymin)/2

        if diameter == None:
            radius = math.ceil(max(xmax-xmin, ymax-ymin, xymax-xymin, 3)/2) + 1
        elif int(diameter) < 0:
            raise Exception('diameter must be a positive integer.')
        else:
            radius = diameter//2 + 1

        xlim = (int(xmid - radius), int(xmid + radius))
        ylim = (int(ymid - radius), int(ymid + radius))

        if symmetric:
            # Draw the plot in a triangular grid.
            fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(3,3), dpi=150)

            # Draw rank for points where SFH has nonzero rank
            s3 = math.sqrt(3)
            for rank in points_by_rank:
                ys = [-3*y - 1 for _, y in points_by_rank[rank]]
                xs = [s3*(2*x+y-1) for x, y in points_by_rank[rank]]
                latex_rank = f'${rank}$'
                ax.scatter(xs, ys, marker='s', color='w', zorder=2)
                ax.scatter(xs, ys, marker=latex_rank, zorder=3)

            # Draw an triangular grid with diameter prescribed above.
            xs, ys = [], []
            for x in range(xlim[0] + 1, xlim[1]):
                for y in range(ylim[0] + 1, ylim[1]):
                    if (math.ceil(xymid - radius + 1)
                            <= x + y <= math.floor(xymid + radius - 1)):
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
                ax.scatter(xs, ys, marker=latex_rank)

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
        return ax

    def SFH_by_rank(self):
        points_by_rank = {}
        for point, rank in self.SFH_ranks.items():
            try:
                points_by_rank[rank].append(point)
            except KeyError:
                points_by_rank[rank] = [point]
        return points_by_rank

class SuturedHandlebody: # not working yet
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

def does_SFH_work():
    for r in range(1,4):
        for s in range(1,4):
            for t in range(1,4):
                M = Genus2Handlebody([(t+r+1,t), (r+s+1,r), (s+t+1,s)])
                if M.SFH(load=True, save=True)[0] == (r+1)*(t+1) + (t+r+1)*s:
                    print(f'Works for r={r}, s={s}, t={t}')
                else:
                    print(f"Doesn't work for r={r}, s={s}, t={t}")
