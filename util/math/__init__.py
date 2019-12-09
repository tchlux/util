# Make the "Fraction" class available here.
from util.math.fraction import Fraction
from util.math.points import mesh, chebyshev, polynomials, \
    polynomial_indices, fekete_indices, fekete_points
from util.math.pairs import pair_to_index, index_to_pair, \
    num_from_pairs, pairwise_distance

# Useful for checking near-equal cases when roundoff error is involved.
SMALL = 2**(-26)
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

# Function for performing absolute difference between numbers (or
# vectors). Falls back to equality for things that don't support
# difference, absolute value, or sums.
def abs_diff(v1, v2):
    if hasattr(v1, "__iter__"):
        try: return sum(abs(v1 - v2))
        except: return int(v1 == v2)
    try:    return abs(v1 - v2)
    except: return int(v1 == v2)

# Function for computing the product of a list of numbers.
def product(sequence):
    v = 1
    for number in sequence: v *= number
    return v

# Compute the probability that 1 specific value is picked from a
# population of size `N` given a set number of `samples`.
def pick(samples, N):
    return 1 - product( ((N-i-1)/(N-i) for i in range(samples)) )

# Return True if a number is prime, False otherwise.
def is_prime(n):
    for i in range(2,int(n**(1/2))+1):
        if (not (n%i)): return False
    return True

# Return all unique prime factors of a number in sorted order.
def primes(n):
    factors = {}
    candidate = 2
    while candidate**2 <= n:
        while (n % candidate) == 0:
            factors[candidate] = factors.get(candidate,0) + 1
            n //= candidate
        candidate += 1
    # Store the last remainder if no more primes factors were found.
    if (n > 1): factors[n] = factors.get(n,0) + 1
    # Sort all prime factors by their index.
    return sorted((k,factors[k]) for k in factors)

# Return all primes up to a given number.
def primes_up_to(n):
    if   (n <= 0): return []
    elif (n == 1): return [1]
    elif (n == 2): return [1,2]
    prime_numbers = [1,2]
    for i in range(3, n+1):
        if is_prime(i): prime_numbers.append( i )
    return prime_numbers

# Function for calculating (n choose k), more efficiently than
# using the raw factorials (takes advantage of cancellation).
# 
# Parameters:
#   n -- positive integer
#   k -- positive integer less than or equal to "n"
# 
# Example:
#   > choose(10,7)
#   120
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

