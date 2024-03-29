#-----------------------------------------------------------------------------#
# Just to make the code easier to read
#-----------------------------------------------------------------------------#

def drop(some_list, indices):
    '''Returns a copy of some_list without the entry some_list[index]'''
    new_list = [elem for k, elem in enumerate(some_list) if not k in indices]
    return new_list

#-----------------------------------------------------------------------------#
# Sets as nonpositive integers
#-----------------------------------------------------------------------------#
# We can represent subsets of the positive integers using nonnegative integers.
# Example: 5 = (101)_2 represents the set {1,3}.
# The functions below define some basic set operations using integers.


def singletons(n):
    '''Yields the powers of two present in the binary expansion of a number.
    Equivalently, it yields the one elements subsets of n (viewed as a set).'''

    while n:
        b = n & (~n+1)
        yield b
        n ^= b


def singletons_complement(n, m):
    '''Yields the one elements subsets of the complement of n in {1,...,m}.'''
    n = (1 << m) + ~n

    while n:
        b = n & (~n+1)
        yield b
        n ^= b


def elements(n):
    '''Yields the positions of the ones present in the binary expansion of a
    number (counting from 1). Equivalently, it yields the elements of a subset
    of the natural numbers.'''

    while n:
        b = n & (~n+1)
        yield b.bit_length()
        n ^= b


def contains(n, k):
    return n & k == k


def d_complement(S, m):
    return reverse((1 << m) + ~S, m)


def crop(n, first, last):
    return (n & (1 << last) - (1 << first-1)) >> first-1


def shift_ones(n, k, m):
    '''
    Find a 1 in the binary expansion and shift that 1 by k digits if there are
    only zeros in between, and the ending position is at or before the m digit.
    It yields the number with the shifted 1 and the power of two corresponding
    to the initial position of the shifted 1.

    Then, starts again with the next 1 in the expansion.

    Ex: shift_one(0b1001_0011, 2, 8) yields
        (0b1001_1001, 0b0000_0010), (0b1100_0011, 0b0000_1000), and stops.
    '''
    a = singletons(crop(n, 1, m) + (1 << m))
    b = next(a)

    while True:
        try:
            c = next(a)
            d = c - b

            if size(d) > k:
                yield (n + (b << k) - b, b)

            b = c

        except StopIteration:
            break


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


def singleton_next_not(n, m):
    '''Yields the powers of two present in the binary expansion that are
    followed by a zero (disregarding digits after the mth bit).

    In terms of subsets, it yield the singleton subsets {n} of a set T such
    that n+1 is not in T.

    Example: singleton_next_not(0b1010_1101, 6) = [0b0000_0001, 0b0000_1000]
    (0b0010_0000 is not included because the following 0 is in the 7th bit)'''

    m = 1 << m
    n &= m - 1
    m = (m >> 1) & n
    while n > m:
        b = n & (~n+1)
        if not b << 1 & n:
            yield b
        n ^= b


def size(n):
    '''Return sum of the digits in the binary expansion.
    Equivalently, it gives the number of elements in the subset.'''
    k = 0
    for _ in elements(n):
        k += 1
    return k


def reverse(n, m):
    return sum(1 << (m - i) for i in elements(crop(n, 1, m)))


def list_subset(n):
    return list(elements(n))


def differ_by_two(S, T):
    ''' This function answer the question: is the cardinality of the symmetric
    difference of those two sets equal to two?'''
    symmetric_difference = S ^ T

    size = 0
    for _ in singletons(symmetric_difference):
        size += 1
        if size > 2:
            return False

    if size < 2:
        return False
    else:
        return True
