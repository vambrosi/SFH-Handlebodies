#-----------------------------------------------------------------------------#
# Genus 2 Handlebodies
#-----------------------------------------------------------------------------#

from constructors import (AAEdge, DDDVertex, DDid)
from tensor import tensor

import math
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------#
# Handlebody class
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

    def SFH(self):
        '''Computes the rank of SFH for all spin_c structures.'''

        # edge_labels gives the number of arcs in the arc diagram of each disk
        # and how much we have to rotate the arc diagram before gluing.

        labels = [(sutures, twist % sutures)
                  for sutures, twist in self.edge_labels]

        # Numbers of arcs for each parametrized disk in the 1st and 2nd ball
        arcs1 = [labels[0][0], labels[1][0], labels[2][0]]
        arcs2 = [labels[0][0], labels[2][0], labels[1][0]]

        # DDDModule for the 1st ball
        S1 = DDDVertex(arcs1[0], arcs1[1], arcs1[2], '+')

        # Add twists
        for sutures, twist in labels:
            T = AAEdge(sutures, twist, '+')
            S1 = tensor(S1, 0, T, 0)

        if self.dual:
            # Code below gives the mirror image of S1
            S2 = DDDVertex(arcs1[0], arcs1[1], arcs1[2], '+').dual()
        else:
            S2 = DDDVertex(arcs2[0], arcs2[1], arcs2[2], '+')

        # Glue two vertices
        M = tensor(S1, 0, S2, 0)
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
            point = M.count_chords(gen)[0: 2]
            try:
                SFH_ranks[point] += 1
            except KeyError:
                SFH_ranks[point] = 1

        self.SFH_ranks = SFH_ranks
        self.SFH_rank = len(N.generators)

        return self.SFH_rank, self.SFH_ranks

    #-------------------------------------------------------------------------#
    # Plotting SFH ranks
    #-------------------------------------------------------------------------#

    def SFH_plot(self, symmetric=True, show=True, filename='', **kwargs):
        points_by_rank = self.SFH_by_rank()

        el = [x//2 for x, _ in self.edge_labels]
        xlim = (-2, el[1] + 1)
        ylim = (-2, el[2] + 1)

        if symmetric:
            # Draw the plot in a triangular grid.
            fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(3, 3), dpi=200)

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
            fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(2, 2), dpi=150)

            # Draw points_by_rank with nonzero rank
            m = max(points_by_rank)
            for rank in points_by_rank:
                xs = [point[0] for point in points_by_rank[rank]]
                ys = [point[1] for point in points_by_rank[rank]]
                latex_rank = f'${rank}$'
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

        if filename:
            plt.savefig(filename)
        if show:
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


class SolidTorus:
    def __init__(self, edge_label):
        assert edge_label[1] % 2 == 0, "Twist must be even."

        self.edge_label = edge_label

        # Does not compute SFH on initialization. Only when you call the
        # function self.SFH() or when you call self.SFH_rank or self.SFH_ranks
        self.SFH_rank = None
        self.SFH_ranks = None

    @property
    def SFH_rank(self):
        '''Total rank of SFH of the SolidTorus.'''
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

    def SFH_by_rank(self):
        points_by_rank = {}
        for point, rank in self.SFH_ranks.items():
            try:
                points_by_rank[rank].append(point)
            except KeyError:
                points_by_rank[rank] = [point]
        return points_by_rank

    def SFH(self):
        '''Computes SFH of the SolidTorus.'''

        sutures, twist = self.edge_label
        M = tensor(AAEdge(sutures, twist, '+'), 1, DDid(sutures), 0)
        N = M.HH(0, 1)

        SFH_ranks = {}
        for gen in N.generators:
            point = M.count_chords(gen)[0]
            try:
                SFH_ranks[point] += 1
            except KeyError:
                SFH_ranks[point] = 1

        self.SFH_ranks = SFH_ranks
        self.SFH_rank = len(N.generators)

        return self.SFH_rank, self.SFH_ranks

    def SFH_plot(self, **kwargs):
        points_by_rank = self.SFH_by_rank()

        el = self.edge_label[0] // 2
        xlim = (-2, el + 1)
        ylim = (-1, 1)

        # Draw the plot in a rectangular grid.
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(4, 2), dpi=150)

        # Draw points_by_rank with nonzero rank
        m = max(points_by_rank)
        for rank in points_by_rank:
            xs = points_by_rank[rank]
            ys = len(xs) * [0]
            latex_rank = f'${rank}$'
            ax.scatter(xs, ys, marker=latex_rank, zorder=3, linewidth=0.1,
                       facecolor=plt.cm.inferno(math.sqrt((rank-1)/m)),
                       edgecolor='black')

        # Draws a square grid with side given by sides.
        xs = []
        for x in range(xlim[0]+1, xlim[1]):
            if not x in self.SFH_ranks:
                xs.append(x)

        ys = len(xs) * [0]

        ax.scatter(xs, ys, marker='.', color='black', s=0.5, alpha=0.5)

        ax.set(**kwargs)
        ax.set_axis_off()
        ax.set(xlim=xlim, ylim=ylim)
        ax.set_aspect('equal')

        plt.show()
        plt.close()


#-----------------------------------------------------------------------------#
# Test Functions
#-----------------------------------------------------------------------------#


def check_SFH_pretzel_knots(rlim: int, slim: int, tlim: int, show=False) -> None:
    '''Consider the complement of the Seifert surface of the pretzel knot
    P(2r+1, 2s+1, 2t+1) given by the checkerboard coloring of the canonical
    projection. This function computes the rank of the SFH of the associated
    sutured manifold and it compares with the result in Example 8.4 of

    S. FRIEDL, A. JUHASZ, J. RASMUSSEN -
    THE DECATEGORIFICATION OF SUTURED FLOER HOMOLOGY
    '''

    # Case where r,s,t > 0
    for r in range(1, rlim + 1):
        for s in range(1, slim + 1):
            for t in range(1, tlim + 1):
                M = Genus2Handlebody([(2*(t+r+1), 2*(t+r+1) - 1),
                                      (2*(r+s+1), 2*(r+s+1) - 1),
                                      (2*(s+t+1), 2*(s+t+1) - 1)])
                if M.SFH(load=False)[0] == (r+1)*(t+1) + (t+r+1)*s:
                    print(f'Works for r={r}, s={s}, t={t}')
                else:
                    print(f"Doesn't work for r={r}, s={s}, t={t}")

                if show:
                    M.SFH_plot()
