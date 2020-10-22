from basics import drop
from strands_algebra_homology import H, I

import math
import csv

class MultiModule:

    def __init__(self, generators={}, arrows={}, action_types=[]):
        self.generators = generators
        self.arrows = arrows
        self.action_types = action_types

        self.ranks = {}
        for gen in self.generators:
            b = self.count_chords(gen)
            try:
                self.ranks[b] += 1
            except KeyError:
                self.ranks[b] = 1

    def __eq__(self, other):
        self.generators == other.generators
        self.arrows == other.arrows
        self.action_types == other.action_types

    def __repr__(self):
        act = [a for (a,b) in self.action_types]
        act = ''.join(act)
        return f'{act}Module with {len(self.generators)} generators'

    #-------------------------------------------------------------------------#
    # User-iterface Functions
    #-------------------------------------------------------------------------#

    def rename_generators(self, permanently=True, format_type='idempotent'):

        if format_type == 'idempotent':
            translation = {ord(','): ord('|'), ord('('): None,
                           ord(')'): None, ord('['): None, ord(']'): None}
            renaming = {gen: str(idp).translate(translation)
                        for gen, idp in self.generators.items()}

        elif format_type == 'alphabet':
            alphabet = [chr(ord('a') + i) for i in range(26)]
            for j in range(len(self.generators)//26):
                alphabet += [chr(ord('a') + i) + str(j) for i in range(26)]
            
            renaming = {gen: alphabet[i] for i, gen 
                        in enumerate(self.generators)}

        generators = {renaming[gen]: self.generators[gen]
                      for gen in self.generators}
        arrows = {renaming[gen]: {} for gen in self.generators}
        
        for gen1 in self.arrows:
            for gen2 in self.arrows[gen1]:
                self.arrows[gen1][gen2]
                arrows[renaming[gen1]][renaming[gen2]] \
                    = self.arrows[gen1][gen2]

        if permanently:
            self.generators = generators
            self.arrows = arrows
        else:
            return generators, arrows

    def show(self, rename=True, format_type='alphabet',
             what='all', chords=None):
        at = self.action_types

        # If you want to rename the generators only for viewing purposes.
        if rename:
            generators, arrows = self.rename_generators(
                                    permanently=False, format_type=format_type)
        else:
            generators = self.generators
            arrows = self.arrows

        # If you want to see only to see generators that have idempotents in
        # A(n1; n2; ...; nm) where strands = [n1, ..., nm].

        if what=='all' or what=='generators':
            overlap = len(self.generators) != len(generators)
            print(f'Generators ({len(self.generators)}, overlap {overlap})')
            
            for gen, idp in generators.items():
                print(f'{gen}\t', idp)
            print()

        if what=='all' or what=='operations':
            print('Operations')
            for start_gen, gen_dict in arrows.items():
                for end_gen, labels in gen_dict.items():
                    for label in labels:
                        op_string = '\u03BC('

                        for k, alg_gen in enumerate(label):              
                            if at[k] == ('A', 'left'):
                                if alg_gen == []:
                                    op_string += '\u2219 , '
                                else:
                                    op_string += f'{alg_gen}, '

                        op_string += f'\u2329{start_gen}\u232A, '

                        for k, alg_gen in enumerate(label):              
                            if at[k] == ('A', 'right'):
                                if alg_gen == []:
                                    op_string += '\u2219 , '
                                else:
                                    op_string += f'{alg_gen}, '

                        op_string = op_string[:-2] + ') = '

                        for k, alg_gen in enumerate(label):              
                            if at[k] == ('D', 'left'):
                                if alg_gen is None or alg_gen.is_idempotent():
                                    op_string += f'\u2219|'
                                else:    
                                    op_string += f'{alg_gen}|'

                        op_string = op_string[:-1] + f' \u2329{end_gen}\u232A '

                        for k, alg_gen in enumerate(label):              
                            if at[k] == ('D', 'right') and not alg_gen is None:
                                if alg_gen is None or alg_gen.is_idempotent():
                                    op_string += f'\u2219|'
                                else:    
                                    op_string += f'{alg_gen}|'

                        op_string = op_string[:-1]

                        print(op_string)
            print()

    #-------------------------------------------------------------------------#
    # Gradings, submodules, and related functions
    #-------------------------------------------------------------------------#

    def count_chords(self, gen):
        return tuple(map(lambda x: x.count_chords(), self.generators[gen]))

    def spinc_submodule(self, chords):
        gens = {gen: self.generators[gen] for gen in self.generators
                if self.count_chords(gen) == chords}
        arrows = {gen: self.arrows[gen] for gen in gens}
        at = self.action_types
        return MultiModule(generators=gens, arrows=arrows, action_types=at)

    def spinc_gradings(self):
        nonempty_gradings = set()
        for gen in self.generators:
            nonempty_gradings.add(self.count_chords(gen))
        return nonempty_gradings

    def gradings_by_rank(self):
        points_by_rank = {}
        for point, rank in self.ranks.items():
            try:
                points_by_rank[rank].append(point)
            except KeyError:
                points_by_rank[rank] = [point]
        return points_by_rank

    #-------------------------------------------------------------------------#
    # Related MultiModules
    #-------------------------------------------------------------------------#

    def dual(self):
        action_types = []

        for (D_A, left_right) in self.action_types:
            if left_right == 'left':
                action_types.append((D_A, 'right'))
            elif left_right == 'right':
                action_types.append((D_A, 'left'))

        gens = self.generators
        arrows = {gen: {} for gen in self.generators}

        for start_gen in self.arrows:
            for end_gen in self.arrows[start_gen]:
                arrows[end_gen][start_gen] = self.arrows[start_gen][end_gen]

        return MultiModule(generators=gens, arrows=arrows,
                           action_types=action_types)

    def HH(self, b1, b2, max_iterations=30, simplify=True):
        if self.action_types[b1][0] == 'D':
            b1, b2 = b2, b1

        action_types = drop(self.action_types, [b1, b2])

        generators = {}
        for gen, idp in self.generators.items():
            if idp[b1] == idp[b2]:
                generators[gen] = drop(idp, [b1, b2])

        arrows = {gen: {} for gen in generators}
        M = MultiModule(generators=generators, arrows=arrows,
                        action_types=action_types)

        for start_gen in generators:
            stack = []

            # First operation
            for end_gen, label in self.arrows_out(start_gen):
                if not label[b1]:
                    if label[b2].is_idempotent():
                        M.add_arrow(start_gen, end_gen, drop(label, [b1, b2]))
                    else:
                        stack.append([label, end_gen])

            # Subsequent operations
            while stack:
                path = stack.pop()
                end_gen = path.pop()

                if len(stack) > max_iterations:
                    raise Exception('Loop timeout')

                A_sequence = []
                D_sequence = []

                for arrow in path:
                    A_sequence.extend(arrow[b1])
                    D_sequence.append(arrow[b2])

                D_sequence = D_sequence[len(A_sequence):]

                for new_end, label in self.arrows_out(end_gen):

                    if label[b1] == D_sequence and label[b2].is_idempotent():
                        new_label = self.compose_arrows(start_gen,
                                                        path + [label], 
                                                        [b1, b2])
                        if not new_label is 0 and new_end in M.generators:
                            M.add_arrow(start_gen, new_end, new_label)

                    if (not label[b2].is_idempotent() 
                        and label[b1] == D_sequence[:len(label[b1])]):
                        stack.append(path + [label, new_end])

        if simplify:
            M.cancel_all_arrows()
        return M
        
    #-------------------------------------------------------------------------#
    # Arrow manipulation
    #-------------------------------------------------------------------------#

    def add_arrow(self, gen1, gen2, label):
        # If label is on list remove it. (Coefficients are mod 2.)
        # If ValueError the label is not on the list, so it adds the label.
        try:
            self.arrows[gen1][gen2].remove(label)
        except ValueError:
            self.arrows[gen1][gen2].append(label)
        except KeyError:
            if gen1 in self.generators and gen2 in self.generators:
                self.arrows[gen1][gen2] = [label]

    def arrows_into(self, end_gen):
        for start_gen in self.arrows:
            try:
                for label in self.arrows[start_gen][end_gen]:
                    yield start_gen, label
            except KeyError:
                pass

    def arrows_out(self, start_gen):
        for end_gen, labels in self.arrows[start_gen].items():
            for label in labels:
                yield end_gen, label

    def compose_arrows(self, start_gen, path, omit):
        at = drop(self.action_types, omit)

        try:
            last = drop(path.pop(), omit)
        except:
            result = []
            idempotents = drop(self.generators[start_gen], omit)
            for action_type, idp in zip(at, idempotents):
                if action_type[0] == 'A':
                    result.append([])
                else:
                    result.append(idp)
            return result.copy()

        while path:
            result = []
            next_last = drop(path.pop(), omit)

            for (DA, lr), arrow_gen1, arrow_gen2 in zip(at, next_last, last):
                if DA == 'A':
                    if lr == 'left':
                        result.append(arrow_gen2 + arrow_gen1)
                    else:
                        result.append(arrow_gen1 + arrow_gen2)
                else:
                    if lr == 'left':
                        result.append(arrow_gen1*arrow_gen2)
                    else:
                        result.append(arrow_gen2*arrow_gen1)

            last = result

        if any(alg_gen is 0 for alg_gen in last):
            return 0
        else:
            return last.copy()

    #-------------------------------------------------------------------------#
    # Cancelation lemma related functions
    #-------------------------------------------------------------------------#

    def find_cancelable_arrow(self):
        at = self.action_types
        for start_gen in self.arrows:
            for end_gen in self.arrows[start_gen]:
                if len(self.arrows[start_gen][end_gen]) <= 1:
                    is_cancelable_arrow = True
                    try:
                        label = self.arrows[start_gen][end_gen][0]
                        for action_type, alg_elems in zip(at, label):
                            if action_type[0] == 'A' and alg_elems:
                                is_cancelable_arrow = False
                                break
                            if (action_type[0] == 'D'
                                    and not alg_elems.is_idempotent()):
                                is_cancelable_arrow = False
                                break
                        if is_cancelable_arrow:
                            return start_gen, end_gen, label
                    except IndexError:
                        return start_gen, end_gen, label

    def remove_generator(self, gen):
        del self.generators[gen]
        del self.arrows[gen]

        for other_gen in self.generators:
            try:
                del self.arrows[other_gen][gen]
            except:
                pass

    def cancel_arrow(self, start_gen, end_gen, label):
        '''Removes one arrow. (According to the Cancelation Lemma.)'''

        # Adds arrows for all nontrivial compositions
        for gen1, label1 in self.arrows_into(end_gen):
            for gen2, label2 in self.arrows_out(start_gen):
                composition = self.compose_arrows(gen1, [label1, label2], [])
                if (not composition is 0
                        and gen1 != start_gen
                        and gen2 != end_gen):
                    self.add_arrow(gen1, gen2, composition)

        # Removes both endpoints of the cancelable arrow and also removes
        # all arrows going into and out of those two generators
        self.remove_generator(start_gen)
        self.remove_generator(end_gen)

    def cancel_all_arrows(self):
        while True:
            try:
                start_gen, end_gen, label = self.find_cancelable_arrow()
                self.cancel_arrow(start_gen, end_gen, label)
            except:
                break

    def is_arrows_consistent(self):
        for gen1 in self.arrows:
            if not gen1 in self.generators: return False
            for gen2 in self.arrows[gen1]:
                if gen2 not in self.generators: return False
        return True

    #-------------------------------------------------------------------------#
    # Box-Tensor related functions
    #-------------------------------------------------------------------------#

    def paths(self, start_gen, component, label):
        if self.action_types[component][0] == 'A':
            raise Exception('Must use a D boundary component.')

        if self.action_types[component][1] == 'right':
            label.reverse()

        paths = []
        path = []
        ends = []
        end_gen = start_gen

        stack = [self.arrows_out(start_gen)]

        while stack:
            try:
                alg_elem = label[len(stack) - 1]
                end_gen, end_label = next(stack[-1])

                if end_label[component] == alg_elem:
                    path.append(end_label)
                    stack.append(self.arrows_out(end_gen))

            except IndexError:
                paths.append(path.copy())
                ends.append(end_gen)
                
                try:
                    path.pop()
                    stack.pop()
                except IndexError:
                    break

            except StopIteration:
                try:
                    path.pop()
                    stack.pop()
                except IndexError:
                    break

        return paths, ends