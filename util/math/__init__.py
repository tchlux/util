# Make the "Fraction" class available here.
from util.math.fraction import Fraction

# Useful for checking near-equal cases when roundoff error is involved.
SMALL = 1.4901161193847656*10**(-8) 
#       ^^ SQRT(EPSILON( 0.0_REAL64 ))

# Given a list of lists, flatten it to one dimension.
def flatten(l): return [v for row in l for v in row]

# Given a list of lists, transpose it and return.
def transpose(l): return [list(row) for row in zip(*l)]

# Return a boolean "is_numeric"
def is_numeric(obj):
    try:
        abs((.3*obj + 1*obj) - .3*obj)
        return True
    except: return False

# Function for performing absolute difference between numbers (or vectors).
def abs_diff(v1, v2):
    if hasattr(v1, "__iter__"):
        try: return sum(abs(v1 - v2))
        except: return int(v1 == v2)
    try:    return abs(v1 - v2)
    except: return int(v1 == v2)

# Return all unique prime factors of a number in sorted order.
def primes(num):
    factors = {}
    candidate = 2
    while candidate**2 <= num:
        while (num % candidate) == 0:
            factors[candidate] = factors.get(candidate,0) + 1
            num //= candidate
        candidate += 1
    if num > 1: factors[num] = factors.get(num,0) + 1
    return sorted((k,factors[k]) for k in factors)

# Check if a variable equals none
def is_none(v): return type(v) == type(None)

# Function for calculating (n choose k), more efficiently than
# using the raw factorials (takes advantage of cancellation).
# 
# Parameters:
#   n -- positive integer
#   k -- positive integer less than or equal to "n"
# 
# Example:
#   > choose(10,7)
#    120
# 
# Notes:
#   Raises generic exception for negative input integers.
#   Undefined behavior for non-integer inputs.
def choose(n, k):
    if (n < 0) or (k < 0): raise(Exception("Both 'n' and 'k' should be positive numbers."))
    if (n == k) or (k == 0): return 1
    numerator = 1
    denominator = 1
    for i in range(n,max(n-k,k),-1): numerator *= i
    for i in range(1,min(n-k,k)+1):  denominator *= i
    return numerator // denominator

