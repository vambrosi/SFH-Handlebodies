#-----------------------------------------------------------------------------#
# Sutured Handlebodies (given by a graph)
#-----------------------------------------------------------------------------#

from multimodules import MultiModule
from multimodules_constructors import (
    AlgebraHomology, CFDD_id, PartialDehnTwist, SolidPairOfPants)
from box_tensor_product import box_tensor

#from sage.plot.scatter_plot import scatter_plot
import matplotlib.pyplot as plt
import csv

class Genus2Handlebody:
    def __init__(self, edge_labels):
        self.edge_labels = edge_labels
        self.SFH_rank = None
        self.SFH_spinc = None

    def SFH(self, method='HH'):
        # Computes the AA_id modules for all n's.
        AA_id = []

        arcs = [n - 1 for n, _ in self.edge_labels]
        for n in arcs:
            DD_dual = CFDD_id(n).dual()
            A = AlgebraHomology(n)
            A_dual = A.dual()
            
            DA = box_tensor(DD_dual, 0, A, 0)   
            AA_id.append(box_tensor(DA, 0, A_dual, 1))

        # Creates one SolidPairOfPants module (the two trimodules are equal) 
        # and take the box tensor products with the identity to get an 
        # AAA module.

        S1 = SolidPairOfPants(arcs[0], arcs[1], arcs[2])
        S1 = box_tensor(S1, 0, AA_id[0], 0)
        S1 = box_tensor(S1, 0, AA_id[1], 0)
        S1 = box_tensor(S1, 0, AA_id[2], 0)

        S2 = S1

        # edge_labels gives how many times R_+ intersects the compressing disks
        # and how much it is twisted along the handles.
        DehnTwists = [PartialDehnTwist(n-1, k, AA_id[i]) for i, (n, k) 
                      in enumerate(self.edge_labels)]

        # Glue Dehn twists to one of the SolidPairOfPants
        S2 = box_tensor(S2, 0, DehnTwists[0], 0)
        S2 = box_tensor(S2, 0, DehnTwists[1], 0)
        S2 = box_tensor(S2, 0, DehnTwists[2], 0)

        # Glue two SolidPairOfPants
        if method == 'HH':
            M = box_tensor(S1, 0, S2, 0, rename=True)
            N = M.HH(0, 2)
            N = N.HH(0, 1)

        elif method == 'box':
            N = box_tensor(S1, [0,1,2], S2, [0,1,2])
        
        self.SFH_rank = len(N.generators)
        return N, M

    def does_SFH_work():
        for r in range(1,3):
            for s in range(1,4):
                for t in range(1,4):
                    H = Genus2Handlebody([(t+r+1,t), (r+s+1,r), (s+t+1,s)])
                    M = H.SFH()
                    if len(M.generators) != (r+1)*(t+1) + (t+r+1)*s:
                        return False
        return True

    def read_polytope_from_csv(self): # fix
        el = self.edge_labels
        name = f'{el[0][0]}-{el[0][1]} ' \
               + f'{el[1][0]}-{el[1][1]} ' \
               + f'{el[2][0]}-{el[2][1]}'

        with open(f'./SFHGenus2/{name}.csv', newline='') as csvfile:
            reader = csv.reader(csvfile, quoting=csv.QUOTE_MINIMAL)

            point_rank = {}
            for row in reader:
                point_rank[(int(row[0]), int(row[1]))] = int(row[2])

        return point_rank

    def SFH_polytope(self, load=True, save=True): 
        if load:
            if self.SFH_spinc != None:
                return self.SFH_spinc

            try:
                point_rank = self.read_polytope_from_csv()
                if point_rank:
                    self.SFH_spinc = point_rank
                    return point_rank
            except FileNotFoundError:
                point_rank = {}    
        else:
            point_rank = {}

        N, M = self.SFH()
        for gen in N.generators:
            point = M.count_chords(gen)[0:2]
            try:
                point_rank[point] += 1
            except KeyError:
                point_rank[point] = 1

        self.SFH_spinc = point_rank
        if save:
            el = self.edge_labels
            name = f'{el[0][0]}-{el[0][1]} ' \
                   + f'{el[1][0]}-{el[1][1]} ' \
                   + f'{el[2][0]}-{el[2][1]}'

            with open(f'./SFHGenus2/{name}.csv', 'w+', newline='') as csvfile:
                writer = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)

                for point, rank in point_rank.items():
                    writer.writerow([point[0], point[1], rank])

        return point_rank

    def SFH_by_rank(self):
        if self.SFH_spinc == None:
            self.SFH_polytope()

        points_by_rank = {}
        for point, rank in self.SFH_spinc.items():
            try:
                points_by_rank[rank].append(point)
            except KeyError:
                points_by_rank[rank] = [point]
        return points_by_rank

    def SFH_plot(self, **kwargs):
        if self.SFH_spinc == None:
            self.SFH_polytope()

        points_by_rank = self.SFH_by_rank()

        x_max = max((point[0] for point in self.SFH_spinc))
        x_min = min((point[0] for point in self.SFH_spinc))
        y_max = max((point[1] for point in self.SFH_spinc))
        y_min = min((point[1] for point in self.SFH_spinc))

        if not 'xlim' in kwargs:
            x_max = max((point[0] for point in self.SFH_spinc))
            x_min = min((point[0] for point in self.SFH_spinc))
            xlim = (x_min-1, x_max+1)
        else:
            xlim = kwargs['xlim']
        
        if not 'ylim' in kwargs:
            y_max = max((point[1] for point in self.SFH_spinc))
            y_min = min((point[1] for point in self.SFH_spinc))
            xlim = (y_min-1, y_max+1)
        else:
            ylim = kwargs['ylim']

        # plots = []
        # for rank in points_by_rank:
        #     plots.append(scatter_plot(points_by_rank[rank], 
        #                               axes=False, figsize=(2,2),
        #                               marker=f'${rank}$',
        #                               xmin=min_both-1, xmax=max_both+1,
        #                               ymin=min_both-1, ymax=max_both+1))

        # join_plot = sum(plots)
        # join_plot.grid(True)
        # return join_plot

        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(2,2), dpi=150)
        cm = plt.cm.get_cmap('tab10')

        for rank in points_by_rank:
            x = [point[0] for point in points_by_rank[rank]]
            y = [point[1] for point in points_by_rank[rank]]
            latex_rank = f'${rank}$'
            ax.scatter(x, y, marker=latex_rank, color=cm.colors[rank-1])

        x, y = [], []
        for xindex in range(xlim[0]+1, xlim[1]):
            for yindex in range(ylim[0]+1, ylim[1]):
                if not (xindex, yindex) in self.SFH_spinc:
                    x.append(xindex)
                    y.append(yindex)

        ax.scatter(x, y, marker='.', color='black', s=0.5, alpha=0.5)
        ax.set_axis_off()

        ax.set(**kwargs)
        ax.set(xlim=xlim, ylim=ylim)
        return ax

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
