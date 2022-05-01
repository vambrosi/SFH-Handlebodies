#-----------------------------------------------------------------------------#
# Sutured Graph MultiModules Constructors.
#-----------------------------------------------------------------------------#
from functools import lru_cache
from basics import shift_ones
from bitsets import (left_rotate, next_one, opposite, complement, singletons,
                     singletons_complement, d_complement, shift_ones_d)

from homology import H, I
from modules import MultiModule
from tensor import tensor


#-----------------------------------------------------------------------------#
# Bimodules
#-----------------------------------------------------------------------------#

@lru_cache(maxsize=64)
def DDid(sutures):
    arcs = sutures//2 - 1
    action_types = [('D', 'left'), ('D', 'left')]

    generators = {}
    for S in range(1 << arcs):
        T = d_complement(S, arcs)
        generators[(S, T)] = [I(S), I(T)]

    arrows = {gen: {} for gen in generators}
    for (S, T) in generators:
        for U, _ in shift_ones(S, 1, arcs):
            V = d_complement(U, arcs)
            arrows[(S, T)][(U, V)] = [[H(S, U), H(T, V)]]

    return MultiModule(generators=generators, arrows=arrows,
                       action_types=action_types)

@lru_cache(maxsize=64)
def AAHomology(sutures):
    # Module parameter
    arcs = sutures // 2 - 1

    # Adding generators
    gens = {}
    for S in range(1 << arcs):
        for U in range(S, 1 << arcs):
            if H(S, U) != 0:
                gens[(S, U)] = [I(S), I(U)]

    # Adding arrows
    arrows = {gen: {} for gen in gens}

    for S, T1 in gens:
        for T2, U in gens:
            if T1 == T2:
                if H(S, U) != 0:
                    arrows[(T1, U)][(S, U)] = [[[H(S, T1)], []]]
                    arrows[(S, T1)][(S, U)] = [[[], [H(T1, U)]]]

        arrows[(S, T1)][(S, T1)].append([[I(S)], []])

    return MultiModule(generators=gens, arrows=arrows,
                       action_types=[('A', 'left'), ('A', 'right')])

@lru_cache(maxsize=64)
def DDEdge(sutures, twist, sign):
    assert sutures > 0 and sutures % 2 == 0

    # Trivial cases
    if twist % sutures == 0:
        if sign == '+':
            return DDid(sutures)
        elif sign == '-':
            return DDid(sutures).dual()
    elif twist % sutures == 1:
        id = DDid(sutures)
        if sign == '+':
            M = tensor(id, 1, AAHomology(sutures).dual(), 0)
            return tensor(M, 1, id.dual(), 0)
        elif sign == '-':
            M = tensor(id.dual(), 1, AAHomology(sutures).dual(), 1)
        return tensor(M, 1, id, 0)

    # Module parameters
    n = sutures // 2
    k = (twist % sutures) // 2
    is_odd = twist % 2

    # Adding generators
    gens = {}

    # Generators are parametrized by tuples (U,V,u,v)
    # u and v will be singleton sets, for simplicity.

    u0 = 1 << n - k

    for A in range(1 << n - k - 1):
        for B in range(1 << k - 1):
            U0 = ((B << n - k) | A) << 1

            # Tuples where u = {n-k}
            U = U0
            u = u0
            V, v = twist_complement(U, u, n, k)

            gens[(U, V, u, v)] = [I(U >> 1), I(V >> 1)]

            # Tuples where u != {n-k}
            U = U0 + u0

            for u in singletons_complement(U, n):
                V, v = twist_complement(U, u, n, k)
                gens[(U, V, u, v)] = [I(U >> 1), I(V >> 1)]

    # Adding arrows
    arrows = {gen: {} for gen in gens}

    for (U1, V1, u1, v1) in gens:
        for (U2, V2, u2, v2) in gens:
            if (U1 < U2 and V1 <= V2) or (V1 < V2 and U1 <= U2):
                if V1 == V2 and next_one(U2, u2) == u1:
                    arrows[(U1, V1, u1, v1)][(U2, V2, u2, v2)] = \
                        [[H(U1 >> 1, U2 >> 1), H(V1 >> 1, V2 >> 1)]]
                elif U1 == U2 and next_one(V2, v2) == v1:
                    arrows[(U1, V1, u1, v1)][(U2, V2, u2, v2)] = \
                        [[H(U1 >> 1, U2 >> 1), H(V1 >> 1, V2 >> 1)]]
                elif u1 == u2 and twist_flip_check(U1, U2, V1, V2):
                    arrows[(U1, V1, u1, v1)][(U2, V2, u2, v2)] = \
                        [[H(U1 >> 1, U2 >> 1), H(V1 >> 1, V2 >> 1)]]

    M = MultiModule(generators=gens, arrows=arrows,
                    action_types=[('D', 'left'), ('D', 'left')])

    if sign == '+':
        if is_odd:
            M = tensor(M, 1, AAHomology(sutures).dual(), 0)
            return tensor(M, 1, DDid(sutures).dual(), 0)
        else:
            return M
    elif sign == '-':
        if is_odd:
            M = tensor(AAHomology(sutures).dual(), 0, M, 1)
            return tensor(DDid(sutures).dual(), 1, M, 0)
        else:
            return M.dual()

@lru_cache(maxsize=64)
def AAEdge(sutures, twist, sign):
    if sign == '+':
        M = DDEdge(sutures, twist, '-')
        M = tensor(AAHomology(sutures), 0, M, 0)
        if twist % 2 == 0:
            return tensor(M, 1, AAHomology(sutures).dual(), 1)
        else:
            return tensor(M, 1, AAHomology(sutures).dual(), 0)
    if sign == '-':
        M = DDEdge(sutures, twist, '+')
        M = tensor(AAHomology(sutures), 1, M, 0)
        if twist % 2 == 0:
            return tensor(M, 1, AAHomology(sutures).dual(), 0)
        else:
            return tensor(M, 1, AAHomology(sutures).dual(), 1)

#-----------------------------------------------------------------------------#
# Trimodule
#-----------------------------------------------------------------------------#

@lru_cache(maxsize=64)
class DDDVertex(MultiModule):
    def __init__(self, sutures1, sutures2, sutures3):
        self.action_types = [('D', 'left'), ('D', 'left'), ('D', 'left')]
        self.arcs = [sutures1//2 - 1, sutures2//2 - 1, sutures3//2 - 1]
        self.alpha_arcs = [(self.arcs[0] - self.arcs[1] + self.arcs[2])//2,
                           (self.arcs[1] - self.arcs[2] + self.arcs[0])//2,
                           (self.arcs[2] - self.arcs[0] + self.arcs[1])//2]

        for arcs in self.alpha_arcs:
            if arcs < 0:
                raise Exception('Invalid parameters.')

        mid_arcs = [1 << arcs for arcs in self.alpha_arcs]

        # Initializing generators and arrows in the MultiModule
        self.generators = {}
        self.arrows = {}

        # Case where R_+ has two triangular regions.
        if int(self.arcs[0] + self.arcs[1] + self.arcs[2]) % 2:

            for S1 in range(mid_arcs[0]):
                U2 = (d_complement(S1, self.alpha_arcs[0])
                      << self.alpha_arcs[2] << 1)

                for T1 in range(mid_arcs[1]):
                    S2 = (d_complement(T1, self.alpha_arcs[1])
                          << self.alpha_arcs[0] << 1)

                    for U1 in range(mid_arcs[2]):
                        T2 = (d_complement(U1, self.alpha_arcs[2])
                              << self.alpha_arcs[1] << 1)

                        S = S1 + mid_arcs[0] + S2
                        T = T1 + mid_arcs[1] + T2
                        U = U1 + U2

                        self.generators[(S, T, U)] = [I(S), I(T), I(U)]
                        self.arrows[(S, T, U)] = {}

                        # Octogonal regions between disk 2 and 0
                        for (A1, C2) in shift_ones_d(S1, self.alpha_arcs[0]):
                            A = A1 + mid_arcs[0] + S2
                            C = U1 + (C2 << self.alpha_arcs[2] << 1)
                            self.arrows[(S, T, U)][(A, T, C)] = \
                                [[H(S, A), I(T), H(U, C)]]

                        # Octogonal regions between disk 0 and 1
                        for (B1, A2) in shift_ones_d(T1, self.alpha_arcs[1]):
                            A = S1 + mid_arcs[0] + \
                                (A2 << self.alpha_arcs[0] << 1)
                            B = B1 + mid_arcs[1] + T2
                            self.arrows[(S, T, U)][(A, B, U)] = \
                                [[H(S, A), H(T, B), I(U)]]

                        # Octogonal regions between disk 1 and 2
                        for (C1, B2) in shift_ones_d(U1, self.alpha_arcs[2]):
                            B = T1 + mid_arcs[1] + \
                                (B2 << self.alpha_arcs[1] << 1)
                            C = C1 + U2
                            self.arrows[(S, T, U)][(S, B, C)] = \
                                [[I(S), H(T, B), H(U, C)]]

                        S = S1 + S2
                        T = T1 + mid_arcs[1] + T2
                        U = U1 + mid_arcs[2] + U2

                        self.generators[(S, T, U)] = [I(S), I(T), I(U)]
                        self.arrows[(S, T, U)] = {}

                        # Octogonal regions between disk 2 and 0
                        for (A1, C2) in shift_ones_d(S1, self.alpha_arcs[0]):
                            A = A1 + S2
                            C = U1 + mid_arcs[2] + \
                                (C2 << self.alpha_arcs[2] << 1)
                            self.arrows[(S, T, U)][(A, T, C)] = \
                                [[H(S, A), I(T), H(U, C)]]

                        # Octogonal regions between disk 0 and 1
                        for (B1, A2) in shift_ones_d(T1, self.alpha_arcs[1]):
                            A = S1 + (A2 << self.alpha_arcs[0] << 1)
                            B = B1 + mid_arcs[1] + T2
                            self.arrows[(S, T, U)][(A, B, U)] = \
                                [[H(S, A), H(T, B), I(U)]]

                        # Octogonal regions between disk 1 and 2
                        for (C1, B2) in shift_ones_d(U1, self.alpha_arcs[2]):
                            B = T1 + mid_arcs[1] + \
                                (B2 << self.alpha_arcs[1] << 1)
                            C = C1 + mid_arcs[2] + U2
                            self.arrows[(S, T, U)][(S, B, C)] = \
                                [[I(S), H(T, B), H(U, C)]]

                        S = S1 + mid_arcs[0] + S2
                        T = T1 + T2
                        U = U1 + mid_arcs[2] + U2

                        self.generators[(S, T, U)] = [I(S), I(T), I(U)]
                        self.arrows[(S, T, U)] = {}

                        # Octogonal regions between disk 2 and 0
                        for (A1, C2) in shift_ones_d(S1, self.alpha_arcs[0]):
                            A = A1 + mid_arcs[0] + S2
                            C = U1 + mid_arcs[2] + \
                                (C2 << self.alpha_arcs[2] << 1)
                            self.arrows[(S, T, U)][(A, T, C)] = \
                                [[H(S, A), I(T), H(U, C)]]

                        # Octogonal regions between disk 0 and 1
                        for (B1, A2) in shift_ones_d(T1, self.alpha_arcs[1]):
                            A = S1 + mid_arcs[0] + \
                                (A2 << self.alpha_arcs[0] << 1)
                            B = B1 + T2
                            self.arrows[(S, T, U)][(A, B, U)] = \
                                [[H(S, A), H(T, B), I(U)]]

                        # Octogonal regions between disk 1 and 2
                        for (C1, B2) in shift_ones_d(U1, self.alpha_arcs[2]):
                            B = T1 + (B2 << self.alpha_arcs[1] << 1)
                            C = C1 + mid_arcs[2] + U2
                            self.arrows[(S, T, U)][(S, B, C)] = \
                                [[I(S), H(T, B), H(U, C)]]

            # Arrows corresponding to the middle regions.
            last_arcs = [arcs >> 1 for arcs in mid_arcs]

            # Middle octogonal region between disk 2 and 0
            for S1 in range(last_arcs[0]):
                for T1 in range(mid_arcs[1]):
                    for U1 in range(mid_arcs[2]):
                        S2 = d_complement(T1, self.alpha_arcs[1])
                        T2 = d_complement(U1, self.alpha_arcs[2])
                        U2 = d_complement(S1, self.alpha_arcs[0]-1)

                        S = S1 + last_arcs[0] + (S2 << self.alpha_arcs[0] << 1)
                        T = T1 + mid_arcs[1] + (T2 << self.alpha_arcs[1] << 1)
                        U = U1 + mid_arcs[2] + (U2 << self.alpha_arcs[2] << 2)

                        A = S + last_arcs[0]
                        C = U + mid_arcs[2]

                        self.arrows[(S, T, U)][(A, T, C)] \
                            = [[H(S, A), I(T), H(U, C)]]

            # Middle octogonal region between disk 0 and 1
            for S1 in range(mid_arcs[0]):
                for T1 in range(last_arcs[1]):
                    for U1 in range(mid_arcs[2]):
                        S2 = d_complement(T1, self.alpha_arcs[1]-1)
                        T2 = d_complement(U1, self.alpha_arcs[2])
                        U2 = d_complement(S1, self.alpha_arcs[0])

                        S = S1 + mid_arcs[0] + (S2 << self.alpha_arcs[0] << 2)
                        T = T1 + last_arcs[1] + (T2 << self.alpha_arcs[1] << 1)
                        U = U1 + mid_arcs[2] + (U2 << self.alpha_arcs[2] << 1)

                        A = S + mid_arcs[0]
                        B = T + last_arcs[1]

                        self.arrows[(S, T, U)][(A, B, U)] \
                            = [[H(S, A), H(T, B), I(U)]]

            # Middle octogonal region between disk 1 and 2
            for S1 in range(mid_arcs[0]):
                for T1 in range(mid_arcs[1]):
                    for U1 in range(last_arcs[2]):
                        S2 = d_complement(T1, self.alpha_arcs[1])
                        T2 = d_complement(U1, self.alpha_arcs[2]-1)
                        U2 = d_complement(S1, self.alpha_arcs[0])

                        S = S1 + mid_arcs[0] + (S2 << self.alpha_arcs[0] << 1)
                        T = T1 + mid_arcs[1] + (T2 << self.alpha_arcs[1] << 2)
                        U = U1 + last_arcs[2] + (U2 << self.alpha_arcs[2] << 1)

                        B = T + mid_arcs[1]
                        C = U + last_arcs[2]

                        self.arrows[(S, T, U)][(S, B, C)] \
                            = [[I(S), H(T, B), H(U, C)]]

        # Case where R_+ has one triangular region.
        else:
            for S1 in range(mid_arcs[0]):
                for T1 in range(mid_arcs[1]):
                    for U1 in range(mid_arcs[2]):
                        S2 = (d_complement(T1, self.alpha_arcs[1])
                              << self.alpha_arcs[0])
                        T2 = (d_complement(U1, self.alpha_arcs[2])
                              << self.alpha_arcs[1])
                        U2 = (d_complement(S1, self.alpha_arcs[0])
                              << self.alpha_arcs[2])

                        S = S1 + S2
                        T = T1 + T2
                        U = U1 + U2

                        self.generators[(S, T, U)] = [I(S), I(T), I(U)]
                        self.arrows[(S, T, U)] = {}

                        # Octogonal regions between disk 2 and 0
                        for (A1, C2) in shift_ones_d(S1, self.alpha_arcs[0]):
                            A = A1 + S2
                            C = U1 + (C2 << self.alpha_arcs[2])
                            self.arrows[(S, T, U)][(A, T, C)] = \
                                [[H(S, A), I(T), H(U, C)]]

                        # Octogonal regions between disk 0 and 1
                        for (B1, A2) in shift_ones_d(T1, self.alpha_arcs[1]):
                            A = S1 + (A2 << self.alpha_arcs[0])
                            B = B1 + T2
                            self.arrows[(S, T, U)][(A, B, U)] = \
                                [[H(S, A), H(T, B), I(U)]]

                        # Octogonal regions between disk 1 and 2
                        for (C1, B2) in shift_ones_d(U1, self.alpha_arcs[2]):
                            B = T1 + (B2 << self.alpha_arcs[1])
                            C = C1 + U2
                            self.arrows[(S, T, U)][(S, B, C)] = \
                                [[I(S), H(T, B), H(U, C)]]

            # Arrows corresponding to the middle region.
            last_arcs = [arcs >> 1 for arcs in mid_arcs]

            for S1 in range(last_arcs[0]):
                U2 = d_complement(S1, self.alpha_arcs[0]-1)

                for T1 in range(last_arcs[1]):
                    S2 = d_complement(T1, self.alpha_arcs[1]-1)

                    for U1 in range(last_arcs[2]):
                        T2 = d_complement(U1, self.alpha_arcs[2]-1)

                        S = S1 + last_arcs[0] + (S2 << self.alpha_arcs[0] << 1)
                        T = T1 + last_arcs[1] + (T2 << self.alpha_arcs[1] << 1)
                        U = U1 + last_arcs[2] + (U2 << self.alpha_arcs[2] << 1)

                        A = S + last_arcs[0]
                        B = T + last_arcs[1]
                        C = U + last_arcs[2]

                        self.arrows[(S, T, U)][(A, B, C)] \
                            = [[H(S, A), H(T, B), H(U, C)]]

        # Rank corresponding to each grading
        self.ranks = {}
        for gen in self.generators:
            b = self.count_chords(gen)
            try:
                self.ranks[b] += 1
            except KeyError:
                self.ranks[b] = 1


#-----------------------------------------------------------------------------#
# Auxiliary functions
#-----------------------------------------------------------------------------#

def twist_complement(U, u, n, k):
    V = complement(left_rotate(opposite(U + u, n), n-k, n), n)
    v = 1 << ((n-k+1-u.bit_length()) % n)
    return V, v


def twist_flip_check(U1, U2, V1, V2):
    U_xor = U1 ^ U2
    V_xor = V1 ^ V2

    U_start = next(singletons(U_xor)).bit_length()
    V_start = next(singletons(V_xor)).bit_length()

    return (U_xor >> U_start - 1) == 0b11 and (V_xor >> V_start - 1) == 0b11
