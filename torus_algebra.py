#-----------------------------------------------------------------------------#
# Torus Algebra (for testing)
#-----------------------------------------------------------------------------#
# We will use the set = numbers convention for consistency

from basics import contains

class TorusAlgebraGenerator:  
    def __init__(self, left_set, right_set):
        self.left_set = left_set
        self.right_set = right_set

    def __eq__(self, other):
        if self is 0:
            return other is 0
        elif other is 0:
            return self is 0
        else:
            return self.left_set == other.left_set \
                and self.right_set == other.right_set

    def __mul__(self, other):
        try:
            if self.right_set != other.left_set:
                if (self.is_idempotent() 
                        and contains(self.right_set, other.left_set)):
                    return other
                elif (other.is_idempotent()
                        and contains(other.left_set, self.right_set)):
                    return self
                else:
                    return 0
            else:
                return TorusAlgebraGenerator(self.left_set, other.right_set)
        
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

    def __repr__(self):
        if self.left_set == 5:
            return 'i'
        elif self.left_set == 10:
            return 'j'
        else:
            first = self.left_set.bit_length()
            last = self.right_set.bit_length() - 1
            if first == 1 and last == 3:
                return '123'
            elif last == first:
                return str(first)
            else:
                return f'{first}{last}'

    def is_idempotent(self):
        return self.left_set == self.right_set