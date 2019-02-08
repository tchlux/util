# Useful for checking near-equal cases when roundoff error is involved.
SMALL_NUMBER = 1.4901161193847656*10**(-8) 
#              ^^ SQRT(EPSILON( 0.0_REAL64 ))

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
    if num > 1: factors.add(num)
    return sorted((k,factors[k]) for k in factors)

# Check if a variable equals none
def is_none(v): return type(v) == type(None)

