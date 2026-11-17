"""Exact integer supercell matrices and coset enumeration.

Only the 3x3 integer matrix and exact rational arithmetic are used; no
floating point value is ever produced.  A matrix S maps the source fractional
column f to lattice-fractional column y through ``S y = f + n`` for an integer
n.  The D = |det S| classes of n in Z^3 / S Z^3 each contain exactly one y in
the half-open cube [0,1)^3; this module enumerates those representatives by a
finite exact search.
"""
from fractions import Fraction
from math import ceil, floor, gcd
from .exact_numbers import read_integer
from .rational_matrix import determinant3

MAX_ENTRY = 3
MAX_DETERMINANT = 32


def read_supercell_matrix(value):
    if type(value) not in (list, tuple) or len(value) != 3:
        raise ValueError("supercell matrix must have three rows")
    rows = []
    for row in value:
        if type(row) not in (list, tuple) or len(row) != 3:
            raise ValueError("supercell matrix rows must have length three")
        rows.append(tuple(read_integer(x, "matrix entry") for x in row))
    matrix = tuple(rows)
    if any(abs(x) > MAX_ENTRY for row in matrix for x in row):
        raise ValueError("matrix entry exceeds the finite bound")
    determinant = determinant3(matrix)
    if determinant == 0 or abs(determinant) > MAX_DETERMINANT:
        raise ValueError("matrix determinant is outside the finite bound")
    return matrix


def integer_adjugate(matrix):
    def minor(row, column):
        rows = [i for i in range(3) if i != row]
        columns = [j for j in range(3) if j != column]
        return (matrix[rows[0]][columns[0]] * matrix[rows[1]][columns[1]]
                - matrix[rows[0]][columns[1]] * matrix[rows[1]][columns[0]])
    cofactors = tuple(tuple((-1) ** (i + j) * minor(i, j) for j in range(3))
                      for i in range(3))
    return tuple(tuple(cofactors[j][i] for j in range(3)) for i in range(3))


def _matvec_integer(matrix, vector):
    return tuple(sum(matrix[i][j] * vector[j] for j in range(3))
                 for i in range(3))


def _lcm(values):
    result = 1
    for value in values:
        result = result * value // gcd(result, value)
    return result


def translation_box(matrix, fractional):
    """Integer box that provably contains every admissible source n.

    For ``y`` in [0,1)^3 each component ``(S y)_i`` lies in the closed interval
    from the sum of the negative entries of row ``i`` to the sum of its positive
    entries.  Hence every integer component of ``n = S y - f`` lies between
    ``ceil(lo_i - f_i)`` and ``floor(hi_i - f_i)``.  The bounds are deliberately
    closed (a superset) because a component can reach either end.
    """
    lower, upper = [], []
    for i, row in enumerate(matrix):
        low = sum(min(0, x) for x in row) - fractional[i]
        high = sum(max(0, x) for x in row) - fractional[i]
        lower.append(ceil(low))
        upper.append(floor(high))
    return tuple(zip(lower, upper))


def coset_key(matrix, translation):
    """Key of ``translation`` in Z^3 / S Z^3, independent of source site."""
    determinant = determinant3(matrix)
    magnitude = abs(determinant)
    sign = 1 if determinant > 0 else -1
    adjugate = integer_adjugate(matrix)
    product = _matvec_integer(adjugate, translation)
    return tuple((sign * product[i]) % magnitude for i in range(3))


def site_representatives(matrix, fractional):
    """Return the D representatives (n, y) of ``S y = fractional + n``.

    ``n`` is the integer triple, ``y`` lies in [0,1)^3 and the n values hit
    every class of Z^3 / S Z^3 exactly once.
    """
    determinant = determinant3(matrix)
    magnitude = abs(determinant)
    sign = 1 if determinant > 0 else -1
    adjugate = integer_adjugate(matrix)
    denominator = _lcm([x.denominator for x in fractional])
    scaled = tuple(x.numerator * (denominator // x.denominator)
                   for x in fractional)
    offset = _matvec_integer(adjugate, scaled)
    target = denominator * magnitude
    boxes = translation_box(matrix, fractional)
    representatives = []
    for n0 in range(boxes[0][0], boxes[0][1] + 1):
        for n1 in range(boxes[1][0], boxes[1][1] + 1):
            for n2 in range(boxes[2][0], boxes[2][1] + 1):
                translation = (n0, n1, n2)
                product = _matvec_integer(adjugate, translation)
                numerator = tuple(sign * (offset[i] + denominator * product[i])
                                  for i in range(3))
                if all(0 <= x < target for x in numerator):
                    y = tuple(Fraction(x, target) for x in numerator)
                    representatives.append((translation, y))
    if len(representatives) != magnitude:
        raise ValueError("supercell representative count differs from determinant")
    keys = {coset_key(matrix, n) for n, _ in representatives}
    if len(keys) != magnitude:
        raise ValueError("supercell representatives repeat a coset")
    representatives.sort(key=lambda item: item[0])
    return tuple(representatives)


def integer_fiber(matrix, start, end, offset):
    """Translation T with ``S T = offset + start - end``, or None.

    ``S T = S y_v - (f_v + e) + e + s - (f_u + s) + ...`` reduces exactly to
    ``offset + start - end`` by the defining equations of the replicas, so a
    unique integral solution exists because the fiber is a coset of Z^3.
    """
    determinant = determinant3(matrix)
    adjugate = integer_adjugate(matrix)
    rhs = tuple(offset[i] + start[i] - end[i] for i in range(3))
    numerator = _matvec_integer(adjugate, rhs)
    if determinant < 0:
        numerator = tuple(-x for x in numerator)
        denominator = -determinant
    else:
        denominator = determinant
    if any(x % denominator for x in numerator):
        return None
    return tuple(x // denominator for x in numerator)
