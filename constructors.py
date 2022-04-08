#-----------------------------------------------------------------------------#
# Sutured Graph MultiModules Constructors.
#-----------------------------------------------------------------------------#
from bitsets import (left_rotate, next_zero, opposite, complement,
                        singletons,     singletons_complement)

from strands_algebra_homology import H, I
from multimodules import MultiModule

def DDEdge(sutures, twist):
    # Twist must be even (odd case will be added later)
    assert twist % 2 == 0

    ## Module parameters
    n = sutures // 2
    k = twist // 2

    ## Adding generators
    gens = {}

    # Generators are parametrized by tuples (U,V,u,v)
    # u and v will be singleton sets, for simplicity.

    u0 = 1 << n - k
    v0 = 1

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
                if V1 == V2 and next_zero(V1, v1) == v2:
                    arrows[(U1, V1, u1, v1)][(U2, V2, u2, v2)] = \
                        [[H(U1 >> 1, U2 >> 1), H(V1 >> 1, V2 >> 1)]]
                elif U1 == U2 and next_zero(U1, u1) == u2:
                    arrows[(U1, V1, u1, v1)][(U2, V2, u2, v2)] = \
                        [[H(U1 >> 1, U2 >> 1), H(V1 >> 1, V2 >> 1)]]
                elif u1 == u2 and twist_flip_check(U1, U2, V1, V2):
                    arrows[(U1, V1, u1, v1)][(U2, V2, u2, v2)] = \
                        [[H(U1 >> 1, U2 >> 1), H(V1 >> 1, V2 >> 1)]]

    return MultiModule(generators=gens, arrows=arrows,
                       action_types=[('D', 'left'), ('D', 'left')])

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