#-----------------------------------------------------------------------------#
# Sutured Graph MultiModules Constructors.
#-----------------------------------------------------------------------------#
from matplotlib.pyplot import box
from bitsets import (left_rotate, next_one, next_zero, opposite, complement,
                        singletons, singletons_complement)

from strands_algebra_homology import H, I
from multimodules import MultiModule
from box_tensor_product import box_tensor

# Temporary fix, final version should be here
from multimodules_constructors import CFDD_id as DDid

#-----------------------------------------------------------------------------#
# Bimodules
#-----------------------------------------------------------------------------#

def AAHomology(sutures):
    ## Module parameter
    arcs = sutures // 2 - 1

    ## Adding generators
    gens = {}
    for S in range(1 << arcs):
        for U in range(S, 1 << arcs):
            if H(S, U) != 0:
                gens[(S, U)] = [I(S), I(U)]

    ## Adding arrows
    arrows = {gen: {} for gen in gens}

    for S in range(1 << arcs):
        for U in range(S, 1 << arcs):
            if H(S, U) != 0:
                for T in range(S, U+1):
                    gen1 = H(S, T)
                    gen2 = H(T, U)

                    if gen1 != 0 and gen2 != 0:
                        arrows[(S, T)][(S, U)] = [[[], [gen2]]]
                        arrows[(T, U)][(S, U)] = [[[gen1], []]]

        arrows[(S, S)][(S, S)].append([[], [I(S)]])

    return MultiModule(generators=gens, arrows=arrows,
                       action_types=[('A', 'left'), ('A', 'right')])


def DDEdge(sutures, twist, sign):
    assert sutures>0 and sutures % 2 == 0

    ## Trivial cases
    if twist % sutures == 0:
        if sign == '+':
            return DDid(sutures//2 - 1)
        elif sign == '-':
            return DDid(sutures//2 - 1).dual()
    elif twist % sutures == 1:
        id = DDid(sutures//2 - 1)
        if sign == '+':
            M = box_tensor(id, 1, AAHomology(sutures).dual(), 0)
            return box_tensor(M, 1, id.dual(), 0)
        elif sign == '-':
            M = box_tensor(id.dual(), 1, AAHomology(sutures).dual(), 1)
        return box_tensor(M, 1, id, 0)

    ## Module parameters
    n = sutures // 2
    k = (twist % sutures) // 2
    is_odd = twist % 2

    ## Adding generators
    gens = {}

    # Generators are parametrized by tuples (U,V,u,v)
    # u and v will be singleton sets, for simplicity.

    u0 = 1 << n - k

    for A in range(1 << n - k - 1):
        for B in range(1 << k - 1):
            U0 = ((B << n - k) | A ) << 1

            # Tuples where u = {n-k}
            U = U0
            u = u0
            V, v = twist_complement(U, u, n, k)

            gens[(U,V,u,v)] = [I(U >> 1), I(V >> 1)]

            # Tuples where u != {n-k}
            U = U0 + u0

            for u in singletons_complement(U, n):
                V, v = twist_complement(U, u, n, k)
                gens[(U,V,u,v)] = [I(U >> 1), I(V >> 1)]

    ## Adding arrows
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
            M = box_tensor(M, 1, AAHomology(sutures).dual(), 0)
            return box_tensor(M, 1, DDid(sutures//2 - 1).dual(), 0)
        else:
            return M
    elif sign == '-':
        if is_odd:
            M = box_tensor(AAHomology(sutures).dual(), 0, M, 1)
            return box_tensor(DDid(sutures//2 - 1).dual(), 1, M, 0)
        else:
            return M.dual()


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