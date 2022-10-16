#-----------------------------------------------------------------------------#
# Strands Algebra Homology
#-----------------------------------------------------------------------------#
from basics import *
from typing import NamedTuple


class HomologyGenerator(NamedTuple):
    '''
    This class encodes the generators of the homology of A(infty,k).

    Generators of A(infty, k) are tuples (S, T, phi) where S and T are finite
    subsets of the positive integers and phi: S -> T is a bijection. Since no
    two strands cross in a nonzero element h of the homology, phi is just the
    unique increasing map between S and T. Thus, h is completely determined by
    its 'left' and 'right' sets (S and T).

    To minimize the memory used subsets of the positive integers will be
    represented by nonnegative integers.

    Example: 5=(101) in base 2, so 5 represents the set {1, 3}.

    Using this convention HomologyGenerator can be used as follows.

    Examples------------------------------------------------------------------

    1) HomologyGenerator(3,5) represents the increasing map 1 2 -> 1 3

    2) HomologyGenerator does not reject 'invalid' inputs like

        HomologyGenerator(3,4)

    which has sets of different size. However,

        HomologyGenerator(3,4).is_zero() returns true
        and HomologyGenerator(3,4) == 0 is also true

    The same holds if we have strands veering downward

        HomologyGenerator(5,3) == 0 is true (5=(101) and 3=(011))

    or if that generator is a boundary (which happens if phi(i) >= j is true
    for some i < j). As an example,

        HomologyGenerator(3,6) == 0 is true (3=(011)_2 and 6=(110)_2)

    3) Products of generators are given by (S,T)*(T,U) = (S,U) with the caveat
    that we have to check if (S,U) is a coboundary.

        H(3,5)*H(5,9) == H(3,9) is true

    Here we used the simplified notation H(n,m) = HomologyGenerator(n,m).

    4) To make multiplication faster, we assume that each term is nonzero
    on the product.

        H(5,3)*H(3,9) == H(5,9) is true

    However, we do check if the result is a boundary.

        H(3,5)*H(5,6) returns 0.

    '''

    left_set: int
    right_set: int

    def __eq__(self, other):
        if self is 0 or self.is_zero():
            return other is 0 or other.is_zero()
        elif other is 0 or other.is_zero():
            return self is 0 or self.is_zero()
        else:
            return self.left_set == other.left_set \
                and self.right_set == other.right_set

    def __repr__(self):
        if self.left_set == 0:
            return '\u2205'
        elif self.is_zero():
            return '0'
        elif self.is_idempotent():
            endpoints = str(list_subset(self.left_set))[1:-1].replace(',', '')
            return f'({endpoints})'
        else:
            left = str(list_subset(self.left_set))[1:-1].replace(',', '')
            right = str(list_subset(self.right_set))[1:-1].replace(',', '')
            return f'({left} \u21A6 {right})'

    def __mul__(self, other):
        try:
            if self.right_set != other.left_set:
                return 0
            else:
                gen = HomologyGenerator(self.left_set, other.right_set)
                if gen.is_zero(do_full_check=False):
                    return 0
                else:
                    return gen

        except AttributeError:
            if other % 2:
                return self
            else:
                return 0

    def __rmul__(self, other):
        # Defines product by scalar on the left.
        try:
            if other % 2:
                return self
            else:
                return 0
        except:
            raise Exception('Multiplication not defined!')

    def __len__(self):
        if size(self.left_set) == size(self.right_set):
            return size(self.left_set)
        else:
            raise Exception('Not well defined.')

    def is_zero(self, do_full_check=True):
        '''The variables same_size and upward_veering are there to check if the
        generator is actually in the strands algebra; overlap checks if the
        generator is a boundary. More details in LOT Bimodules paper.'''

        left = list_subset(self.left_set)
        right = list_subset(self.right_set)

        overlap = any(a <= b for a, b in zip(left[1:], right[:-1]))

        if do_full_check:
            same_size = len(left) == len(right)
            upward_veering = all(a <= b for a, b in zip(left, right))
            return not same_size or not upward_veering or overlap

        return overlap

    def is_idempotent(self):
        return self.left_set == self.right_set

    def count_chords(self):
        return size(self.left_set)

#-----------------------------------------------------------------------------#
# Functions to simplify the notation. They will be used in the other files.
#-----------------------------------------------------------------------------#


def I(subset):
    return HomologyGenerator(subset, subset)


def H(left_set, right_set):
    return HomologyGenerator(left_set, right_set)
