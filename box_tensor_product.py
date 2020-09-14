#-----------------------------------------------------------------------------#
# Box tensor product
#-----------------------------------------------------------------------------#
from basics import drop
from multimodules import MultiModule

def box_tensor(M, boundariesM, N, boundariesN, rename=False, simplify=True):
    ### Doesn't work for more than one boundary ###

    '''
    Returns the box tensor product of the MultiModules M and N. bM and bN are 
    integers that denote the components that will be glued. The code assumes 
    that the bM action on M is on the right and that the bN action on N is on 
    the left.
    '''

    # In case there is just one boundary
    if isinstance(boundariesM, int) and isinstance(boundariesN, int):
        boundariesM = [boundariesM]
        boundariesN = [boundariesN]

    # The resulting module has actions carried over from to the ones in M 
    # and N (we remove the glued boundary components from the list).
    action_types = (drop(M.action_types, boundariesM) 
                    + drop(N.action_types, boundariesN))

    # Generators are given by tensor products over the idempotent ring, so
    # they are pairs of generators (gen1, gen2) where the right idempotent
    # of gen1 agrees with the left idempotent of gen2.

    boundaries = list(zip(boundariesM, boundariesN))

    generators = {}
    for genM, idpM in M.generators.items():
        for genN, idpN in N.generators.items():
            if all(idpM[bM] == idpN[bN] for bM, bN in boundaries):
                generators[(genM, genN)] = (drop(idpM, boundariesM)
                                            + drop(idpN, boundariesN))

    arrows = {gen: {} for gen in generators}
    P = MultiModule(generators=generators, arrows=arrows,
                    action_types=action_types)

    for bM, bN in boundaries:
        if M.action_types[bM][0] == 'D' and N.action_types[bN][0] == 'A':
            for (genM1, genN1) in generators:
                for genN2, label in N.arrows_out(genN1):
                    paths, ends = M.paths(genM1, bM, label[bN])

                    for path, genM2 in zip(paths, ends):
                        composite = M.compose_arrows(genM1, path, boundariesM)
                        if not composite is 0:
                            P.add_arrow((genM1, genN1), (genM2, genN2),
                                        composite + drop(label, boundariesN))

        elif M.action_types[bM][0] == 'A' and N.action_types[bN][0] == 'D':
            for (genM1, genN1) in generators:
                for genM2, label in M.arrows_out(genM1):
                    paths, ends = N.paths(genN1, bN, label[bM])

                    for path, genN2 in zip(paths, ends):
                        composite = N.compose_arrows(genN1, path, boundariesN)
                        if not composite is 0:
                            P.add_arrow((genM1, genN1), (genM2, genN2),
                                        drop(label, boundariesM) + composite)

    if rename:
        P.rename_generators(permanently=True, format_type='alphabet')
    if simplify:
        P.cancel_all_arrows()

    return P
