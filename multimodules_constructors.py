#-----------------------------------------------------------------------------#
# Main MultiModules constructors.
#-----------------------------------------------------------------------------#
from basics import (d_complement, shift_ones, shift_ones_d,
                    singletons, reverse, crop)

from strands_algebra_homology import H, I
from multimodules import MultiModule
from box_tensor_product import box_tensor

#-----------------------------------------------------------------------------#
# Strands Algebra related modules
#-----------------------------------------------------------------------------#

# Upper case letters represent sets (which correspond to idempotents).

# This file contains most of the nontrivial constructions. It will not make
# any sense unless you saw the Heegaard diagrams of the relevant manifolds.
# All manifolds here are partially sutured balls with two or three parametrized
# disks on their boundaries.

def AlgebraHomology(n):
    '''Returns homology of the strands algebra as an AA Module.'''
    action_types = [('A','left'), ('A','right')]

    generators = {}
    for S in range(1 << n):
        for U in range(S, 1 << n):
            if H(S, U) != 0:
                generators[(S, U)] = [I(S), I(U)]

    arrows = {gen: {} for gen in generators}

    for S in range(1 << n):
        for U in range(S, 1 << n):
            if H(S, U) != 0:
                for T in range(S, U+1):
                    gen1 = H(S, T)
                    gen2 = H(T, U)

                    if gen1 != 0 and gen2 != 0:
                        arrows[(S, T)][(S,U)] = [[[], [gen2]]]
                        arrows[(T, U)][(S,U)] = [[[gen1], []]]

        arrows[(S,S)][(S,S)].append([[], [I(S)]])

    return MultiModule(generators=generators, arrows=arrows,
                       action_types=action_types)

def CFDD_id(arcs):
    action_types = [('D','left'), ('D','left')]

    generators = {}
    for S in range(1 << arcs):
        T = d_complement(S, arcs)
        generators[(S,T)] = [I(S), I(T)]

    arrows = {gen: {} for gen in generators}
    for (S, T) in generators:
        for U, _ in shift_ones(S, 1, arcs):
            V = d_complement(U, arcs)
            arrows[(S,T)][(U,V)] = [[H(S,U), H(T,V)]]

    return MultiModule(generators=generators, arrows=arrows,
                       action_types=action_types)

def CFAA_id(arcs):
    A = AlgebraHomology(arcs)
    A_dual = A.dual()

    DA = box_tensor(CFDD_id(arcs).dual(), 0, A, 0)
    AA_id = box_tensor(DA, 0, A_dual, 1)

    return AA_id

#-----------------------------------------------------------------------------#
# Handlebody Modules
#-----------------------------------------------------------------------------#

def PartialDehnTwist(arcs, power, AA_id=None):
    action_types = [('D','left'), ('D','left')]

    # If power is 0 the partial Dehn twist is the identity.
    if power == 0 or arcs == 0:
        return CFDD_id(arcs)

    # Adding generators
    generators = {}
    n = arcs - 1
    for S in range(1 << n):
        T = d_complement(S, n)
        U = (1<<n) + S
        V = (1<<n) + T

        generators[(U, T)] = [I(U), I(T)]
        for k in singletons(U):
            generators[(U - k, V)] = [I(U - k), I(V)]

        
    arrows = {gen: {} for gen in generators}

    # Sums of the rightmost octogonal domains. (Sums of B_k's k<n.)
    for T in range(1 << n):
        V = (1<<n) + T
        for shift in range(1, n):
            for (W, b) in shift_ones(V, shift, n):
                U = d_complement(W + b, n) + (1<<n)
                arrows[(U, V)][(U, W)] = [[I(U), H(V,W)]]

    # Sums of rightmost domains that include B_n. 
    # (Sums of B_k's including B_n.)
    for V in range(1, 1 << n):
        U = (1<<n) + d_complement(V, n)
        W = (1<<n) + V - (1 << V.bit_length()-1)
        arrows[(U,V)][(U, W)] = [[I(U), H(V,W)]]

    # Sums of leftmost domains. (sums of A_k's)
    for S in range(1 << n):
        T = d_complement(S, n)
        U = (1<<n) + S
        V = (1<<n) + T

        for a in singletons(U):
            for b in singletons(U):
                if b > a:
                    A = U - a
                    B = U - b
                    alg_elem = H(B,A)
                    if not alg_elem.is_zero(do_full_check=False):
                        arrows[(B,V)][(A, V)] = [[alg_elem, I(V)]]

    # Sums of domains crossing over (A_{n-1-k} + B_k for k<n.)
    for (A, V) in generators:
        for (W, b) in shift_ones(V, 1, n):
            a = reverse(b << 1, n)
            if A | a == A:
                B = A + a
                arrows[(A,V)][(B, W)] = [[H(A,B), H(V,W)]]

    # In case you precomputed AA_id
    if not isinstance(AA_id, MultiModule):
        AA_id = CFAA_id(arcs)

    DD = MultiModule(generators=generators, arrows=arrows,
                     action_types=action_types)
    
    DA = box_tensor(DD, 1, AA_id, 0)
    for i in range(power - 1):
        DD = box_tensor(DA, 1, DD, 0)
        DD.cancel_all_arrows()

    return DD

class SolidPairOfPants(MultiModule):
    def __init__(self, arcs1, arcs2, arcs3):
        self.action_types = [('D','left'), ('D','left'), ('D','left')]
        self.arcs = [arcs1, arcs2, arcs3]
        self.alpha_arcs = [(self.arcs[0] - self.arcs[1] + self.arcs[2])//2, \
                           (self.arcs[1] - self.arcs[2] + self.arcs[0])//2, \
                           (self.arcs[2] - self.arcs[0] + self.arcs[1])//2]

        for arcs in self.alpha_arcs:
            if arcs < 0: raise Exception('Invalid parameters.')

        mid_arcs = [1 << arcs for arcs in self.alpha_arcs]

        # Initializing generators and arrows in the MultiModule
        self.generators = {}
        self.arrows = {}

        ## Case where R_+ has two triangular regions.
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

                        self.generators[(S,T,U)] = [I(S),I(T),I(U)]
                        self.arrows[(S,T,U)] = {}

                        # Octogonal regions between disk 2 and 0
                        for (A1, C2) in shift_ones_d(S1, self.alpha_arcs[0]):
                            A = A1 + mid_arcs[0] + S2
                            C = U1 + (C2<<self.alpha_arcs[2]<<1)
                            self.arrows[(S,T,U)][(A,T,C)] = \
                                [[H(S,A), I(T), H(U,C)]]

                        # Octogonal regions between disk 0 and 1
                        for (B1, A2) in shift_ones_d(T1, self.alpha_arcs[1]):
                            A = S1 + mid_arcs[0] + (A2<<self.alpha_arcs[0]<<1)
                            B = B1 + mid_arcs[1] + T2
                            self.arrows[(S,T,U)][(A,B,U)] = \
                                [[H(S,A), H(T,B), I(U)]]

                        # Octogonal regions between disk 1 and 2
                        for (C1, B2) in shift_ones_d(U1, self.alpha_arcs[2]):
                            B = T1 + mid_arcs[1] + (B2<<self.alpha_arcs[1]<<1)
                            C = C1 + U2
                            self.arrows[(S,T,U)][(S,B,C)] = \
                                [[I(S), H(T,B), H(U,C)]]

                        S = S1 + S2
                        T = T1 + mid_arcs[1] + T2
                        U = U1 + mid_arcs[2] + U2

                        self.generators[(S,T,U)] = [I(S),I(T),I(U)]
                        self.arrows[(S,T,U)] = {}

                        # Octogonal regions between disk 2 and 0
                        for (A1, C2) in shift_ones_d(S1, self.alpha_arcs[0]):
                            A = A1 + S2
                            C = U1 + mid_arcs[2] + (C2<<self.alpha_arcs[2]<<1)
                            self.arrows[(S,T,U)][(A,T,C)] = \
                                [[H(S,A), I(T), H(U,C)]]

                        # Octogonal regions between disk 0 and 1
                        for (B1, A2) in shift_ones_d(T1, self.alpha_arcs[1]):
                            A = S1 + (A2<<self.alpha_arcs[0]<<1)
                            B = B1 + mid_arcs[1] + T2
                            self.arrows[(S,T,U)][(A,B,U)] = \
                                [[H(S,A), H(T,B), I(U)]]

                        # Octogonal regions between disk 1 and 2
                        for (C1, B2) in shift_ones_d(U1, self.alpha_arcs[2]):
                            B = T1 + mid_arcs[1] + (B2<<self.alpha_arcs[1]<<1)
                            C = C1 + mid_arcs[2] + U2
                            self.arrows[(S,T,U)][(S,B,C)] = \
                                [[I(S), H(T,B), H(U,C)]]

                        S = S1 + mid_arcs[0] + S2
                        T = T1 + T2
                        U = U1 + mid_arcs[2] + U2

                        self.generators[(S,T,U)] = [I(S),I(T),I(U)]
                        self.arrows[(S,T,U)] = {}

                        # Octogonal regions between disk 2 and 0
                        for (A1, C2) in shift_ones_d(S1, self.alpha_arcs[0]):
                            A = A1 + mid_arcs[0] + S2
                            C = U1 + mid_arcs[2] + (C2<<self.alpha_arcs[2]<<1)
                            self.arrows[(S,T,U)][(A,T,C)] = \
                                [[H(S,A), I(T), H(U,C)]]

                        # Octogonal regions between disk 0 and 1
                        for (B1, A2) in shift_ones_d(T1, self.alpha_arcs[1]):
                            A = S1 + mid_arcs[0] + (A2<<self.alpha_arcs[0]<<1)
                            B = B1 + T2
                            self.arrows[(S,T,U)][(A,B,U)] = \
                                [[H(S,A), H(T,B), I(U)]]

                        # Octogonal regions between disk 1 and 2
                        for (C1, B2) in shift_ones_d(U1, self.alpha_arcs[2]):
                            B = T1 + (B2<<self.alpha_arcs[1]<<1)
                            C = C1 + mid_arcs[2] + U2
                            self.arrows[(S,T,U)][(S,B,C)] = \
                                [[I(S), H(T,B), H(U,C)]]

            # Arrows corresponding to the middle regions.
            last_arcs = [arcs >> 1 for arcs in mid_arcs]

            # Middle octogonal region between disk 2 and 0
            for S1 in range(last_arcs[0]):
                for T1 in range(mid_arcs[1]):
                    for U1 in range(mid_arcs[2]):
                        S2 = d_complement(T1, self.alpha_arcs[1])
                        T2 = d_complement(U1, self.alpha_arcs[2])
                        U2 = d_complement(S1, self.alpha_arcs[0]-1)

                        S = S1 + last_arcs[0] + (S2<<self.alpha_arcs[0]<<1)
                        T = T1 + mid_arcs[1] + (T2<<self.alpha_arcs[1]<<1)
                        U = U1 + mid_arcs[2] + (U2<<self.alpha_arcs[2]<<2)

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

                        S = S1 + mid_arcs[0] + (S2<<self.alpha_arcs[0]<<2)
                        T = T1 + last_arcs[1] + (T2<<self.alpha_arcs[1]<<1)
                        U = U1 + mid_arcs[2] + (U2<<self.alpha_arcs[2]<<1)

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

                        S = S1 + mid_arcs[0] + (S2<<self.alpha_arcs[0]<<1)
                        T = T1 + mid_arcs[1] + (T2<<self.alpha_arcs[1]<<2)
                        U = U1 + last_arcs[2] + (U2<<self.alpha_arcs[2]<<1)

                        B = T + mid_arcs[1]
                        C = U + last_arcs[2]

                        self.arrows[(S, T, U)][(S, B, C)] \
                            = [[I(S), H(T, B), H(U, C)]]

        ## Case where R_+ has one triangular region.
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
                        self.arrows[(S,T,U)] = {}

                        # Octogonal regions between disk 2 and 0
                        for (A1, C2) in shift_ones_d(S1, self.alpha_arcs[0]):
                            A = A1 + S2
                            C = U1 + (C2<<self.alpha_arcs[2])
                            self.arrows[(S,T,U)][(A,T,C)] = \
                                [[H(S,A), I(T), H(U,C)]]

                        # Octogonal regions between disk 0 and 1
                        for (B1, A2) in shift_ones_d(T1, self.alpha_arcs[1]):
                            A = S1 + (A2<<self.alpha_arcs[0])
                            B = B1 + T2
                            self.arrows[(S,T,U)][(A,B,U)] = \
                                [[H(S,A), H(T,B), I(U)]]

                        # Octogonal regions between disk 1 and 2
                        for (C1, B2) in shift_ones_d(U1, self.alpha_arcs[2]):
                            B = T1 + (B2<<self.alpha_arcs[1])
                            C = C1 + U2
                            self.arrows[(S,T,U)][(S,B,C)] = \
                                [[I(S), H(T,B), H(U,C)]]

            # Arrows corresponding to the middle region.
            last_arcs = [arcs >> 1 for arcs in mid_arcs]

            for S1 in range(last_arcs[0]):
                U2 = d_complement(S1, self.alpha_arcs[0]-1)

                for T1 in range(last_arcs[1]):
                    S2 = d_complement(T1, self.alpha_arcs[1]-1)

                    for U1 in range(last_arcs[2]):
                        T2 = d_complement(U1, self.alpha_arcs[2]-1)

                        S = S1 + last_arcs[0] + (S2<<self.alpha_arcs[0]<<1)
                        T = T1 + last_arcs[1] + (T2<<self.alpha_arcs[1]<<1)
                        U = U1 + last_arcs[2] + (U2<<self.alpha_arcs[2]<<1)

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