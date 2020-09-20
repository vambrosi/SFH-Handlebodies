#-----------------------------------------------------------------------------#
# Main MultiModules constructors.
#-----------------------------------------------------------------------------#
from basics import d_complement, shift_ones, singletons, reverse, crop

from strands_algebra_homology import H, I
from multimodules import MultiModule
from box_tensor_product import box_tensor

#-----------------------------------------------------------------------------#
# Strands Algebra related modules
#-----------------------------------------------------------------------------#

# Upper case letters represent sets (which correspond to idempotents).
# Lower case letters represent elements of the corresponding set.

def AlgebraHomology(n):
    '''Returns homology of the strands algebra as an AA Module.'''
    action_types = [('A','left'), ('A','right')]

    generators = {}
    for S in range(1 << n):
        for U in range(S, 1 << n):
            if not H(S, U).is_zero():
                generators[(S, U)] = [I(S), I(U)]

    arrows = {gen: {} for gen in generators}

    for S in range(1 << n):
        for U in range(S, 1 << n):
            if not H(S, U).is_zero():
                for T in range(S, U+1):
                    gen1 = H(S, T)
                    gen2 = H(T, U)

                    if not gen1.is_zero():
                        arrows[(S, T)][(S,U)] = [[[], [gen2]]]
                        
                    if not gen2.is_zero():
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

    DA = box_tensor(CFDD_id(arcs).dual(), 0, 0, A)
    AA_id = box_tensor(DA, 0, 1, A_dual)
    AA_id.cancel_all_arrows()

    return AA_id

#-----------------------------------------------------------------------------#
# Handlebody Modules
#-----------------------------------------------------------------------------#

def PartialDehnTwist(arcs, power, AA_id=None):
    action_types = [('D','left'), ('D','left')]

    # If power is 0 the partial Dehn twist is the identity.
    if not power:
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
        # TO DO: less repetition on the for loops (to make it faster).
        self.action_types = [('D','left'), ('D','left'), ('D','left')]
        self.arcs = [arcs1, arcs2, arcs3]
        self.alpha_arcs = [(arcs1 - arcs2 + arcs3)//2, \
                           (arcs2 - arcs3 + arcs1)//2, \
                           (arcs3 - arcs1 + arcs2)//2]

        for arcs in self.alpha_arcs:
            if arcs < 0:
                # Temporary until I add the nongeneric case.
                raise Exception('Invalid parameters.')

        if (arcs1 + arcs2 + arcs3) % 2:
            # Temporary while this case is not ready.
            raise Exception('Invalid parameters.')

            # Assumes you want the triangular region to be in R_+
            # TO DO: add an extra parameter in case you want it to be in R_-

            # Adding generators 

        else:
            # Adding generators
            self.generators = {}
            for S in range(1<<arcs1):
                V = d_complement(crop(S, self.alpha_arcs[0]+1, arcs1), 
                                 self.alpha_arcs[1])
                W = d_complement(crop(S,1,self.alpha_arcs[0]),
                                 self.alpha_arcs[0]) << self.alpha_arcs[2]

                for i in range(1<<self.alpha_arcs[2]):
                    T = V + (i << self.alpha_arcs[1])
                    U = d_complement(i, self.alpha_arcs[2]) + W
                    self.generators[(S,T,U)] = [I(S), I(T), I(U)]

            # Adding arrows
            self.arrows = {gen:{} for gen in self.generators}
            for (S, T, U) in self.generators:
                # Octogonal regions between disk 2 and 0
                for (V, b) in shift_ones(S, 1, self.alpha_arcs[0]):
                    a = reverse(b,self.alpha_arcs[0])<<self.alpha_arcs[2] >> 1
                    self.arrows[(S, T, U)][(V, T, U + a)] \
                        = [[H(S, V), I(T), H(U, U+a)]]

                # Octogonal regions between disk 0 and 1
                for (V, b) in shift_ones(T, 1, self.alpha_arcs[1]):
                    a = reverse(b,self.alpha_arcs[1])<<self.alpha_arcs[0] >> 1
                    self.arrows[(S, T, U)][(S + a, V, U)] \
                        = [[H(S, S+a), H(T, V), I(U)]]

                # Octogonal regions between disk 1 and 2
                for (V, b) in shift_ones(U, 1, self.alpha_arcs[2]):
                    a = reverse(b,self.alpha_arcs[2])<<self.alpha_arcs[1] >> 1
                    self.arrows[(S, T, U)][(S, T + a, V)] \
                        = [[I(S), H(T, T+a), H(U, V)]]

            # Arrows corresponding to the middle region (with 12 sides).

            # a1 a2 a3 represent the arcs that must occupied in the 
            a1 = 1<<self.alpha_arcs[0] >> 1
            a2 = 1<<self.alpha_arcs[1] >> 1
            a3 = 1<<self.alpha_arcs[2] >> 1

            for S1 in range(a1):
                for T1 in range(a2):
                    for U1 in range(a3):
                        S2 = d_complement(T1, self.alpha_arcs[1]-1)
                        T2 = d_complement(U1, self.alpha_arcs[2]-1)
                        U2 = d_complement(S1, self.alpha_arcs[0]-1)

                        S = S1 + a1 + (S2 << self.alpha_arcs[0]+1)
                        T = T1 + a2 + (T2 << self.alpha_arcs[1]+1)
                        U = U1 + a3 + (U2 << self.alpha_arcs[2]+1)

                        self.arrows[(S, T, U)][(S+a1, T+a2, U+a3)] \
                            = [[H(S, S+a1), H(T, T+a2), H(U, U+a3)]]

        # Rank corresponding to each grading
        self.ranks = {}
        for gen in self.generators:
            b = self.count_chords(gen)
            try:
                self.ranks[b] += 1
            except KeyError:
                self.ranks[b] = 1