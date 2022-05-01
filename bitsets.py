#-----------------------------------------------------------------------------#
# Bitset Functions
#-----------------------------------------------------------------------------#

from numpy import single


def left_rotate(bitset, shift, max_bits):
    '''
    Bitwise left rotate bitset by shift wrapping around max_bits.

    left_rotate(0b010110, 3, 6) --> 0b110010
    '''
    all_ones = (1 << max_bits) - 1

    return ((bitset << (shift % max_bits)) & all_ones) | \
        ((bitset & all_ones) >> ((max_bits - shift) % max_bits))

def reverse(bitset, max_bits):
    '''
    Reverses order of bits in bitset. First bit goes to max_bits and so on.

    reverse(0b010110, 6) --> 0b011010
    '''
    r = 0
    b = bitset
    n = max_bits

    for i in range(n):
        r = (r << 1) | (b & 1)
        b >>= 1

    return r

def opposite(bitset, max_bits):
    return (reverse(bitset >> 1, max_bits - 1) << 1) | (bitset & 1)

def complement(bitset, max_bits):
    '''
    Complement of a bitset, i.e. swaps 1s and 0s until max_bits.

    complement(0b010110, 6) --> 0b101001
    '''
    return ((1 << max_bits) - 1) ^ bitset

def singletons(n0):
    '''Yields the powers of two present in the binary expansion of a number.
    Equivalently, it yields the one elements subsets of n0 (viewed as a set).'''
    n = n0

    while n:
        b = n & (~n+1)
        yield b
        n ^= b

def singletons_complement(n0, m):
    '''Yields the one elements subsets of the complement of n0 in {1,...,m}.'''
    n = (1 << m) + ~n0

    while n:
        b = n & (~n+1)
        yield b
        n ^= b

def next_one(n0, start):
    n = n0 & ~((start << 1) - 1)
    return n & (~n+1)

def d_complement(S, m):
    return reverse((1 << m) + ~S, m)

def crop(n, first, last):
    return (n & (1 << last) - (1 << first-1)) >> first-1

def size(n):
    '''Return sum of the digits in the binary expansion.
    Equivalently, it gives the number of elements in the subset.'''
    k = 0
    for _ in singletons(n):
        k += 1
    return k

def shift_ones_d(S, m):
    a = singletons(crop(S, 1, m) + (1 << m))
    b = next(a)

    while True:
        try:
            c = next(a)
            d = c - b

            if size(d) > 1:
                T = S + (b << 1) - b
                yield (T, d_complement(T, m))

            b = c

        except StopIteration:
            break