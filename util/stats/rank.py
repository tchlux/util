from util.stats import *

import os, pickle, time
from itertools import combinations

# Pre:  "indexable" is a list type, "value" is an indexable whose
#       elements have comparison operators, "key" is a function to use
#       for sorting
# Post: A binary search over the indexable is performed and "value" is
#       inserted in a sorted order left of anything with equal key,
#       "key" defaults to sorting by direct element value comparison.
def insert_sorted(indexable, value, key=lambda i:i):
    low = 0
    high = len(indexable)
    index = (low + high) // 2
    while high != low:
        if key(indexable[index]) >= key(value):   high = index
        elif key(indexable[index]) < key(value):  low = index + 1
        index = (low + high) // 2
    indexable.insert(index, value)

# Pre:  "numbers" is a python list of numbers
# Post: All of the values in "numbers" multiplied together
def product(numbers):
    value = 1
    for n in numbers:
        value *= n
    return value

# Pre:  "indexable" is an object that can be index accessed
#       "indexable" must be in sorted order
#       "value" is a value that can be compared with values in "indexable"
# Post: The number of elements in "indexable" that are greater than "value"
def count_greater(indexable, value):
    low = 0
    high = len(indexable)
    index = (low + high) // 2    
    while high != low:
        if   indexable[index] >  value:
            high = index
        elif indexable[index] <= value:  
            low  = index + 1
        else:
            string = "High: %i\tLow: %i\tIndex: %i\n"%(high,low,index)
            string += "Indexable:\n%s"%str(indexable)
            raise(Exception("Encountered invalid number in "+
                            "comparison.\n\n"+string))
        index = (low + high) // 2
    return len(indexable) - index

# Pre:  "interest" is an iterable of length > 0
#       "other_data" is a list of iterables each of length > 0
#       "rank" is a 0-indexed position in an ordering, 0'th being best
#       "order" is 'True' if order matters for data set values.
# Post: The probability that if random single elements are selected
#       from each of the sets in "other_data" and "interest", that the
#       value selected from "interest" will be in position "rank" from
#       the smallest value
def rank_probability(interest, other_data, rank=0, order=False):
    # Check basic assumptions that should hold at the start of this function.
    assert(len(interest) > 0)
    assert(len(other_data) > 0)
    assert(sum(map(len,other_data)) > 0)
    # Convert negative ranks into their positive counterparts.
    while (rank < 0): rank += len(other_data) + 1
    # Prepare all the data in tuple format (less memory) make sure
    # that all lists in other data are sorted for improved runtime
    interest = sorted(interest)
    other_data = [sorted(d) for d in other_data]

    # If order matters, cycle through the data once computing
    # likelihood and exit.
    if order:
        length = len(interest)
        if not all(l == length for l in map(len,other_data)):
            raise(Exception("All provided value series must have same length."))
        # Initialize holder for rank occurrence.
        rank_count = 0
        for group in zip(interest, *other_data):
            interest_value = group[0]
            order = sorted(group)
            rank_count += int(order.index(interest_value) == rank)
        # return the probability that "interest" takes "rank" among peers.
        return rank_count / length

    # If order does *not* matter, then use combinatorics to determine
    # the probability that any value in "interest" is the lowest of
    # all values.

    # Initialize list for holding all probabilities    
    probabilities = []
    # Calculate the total different selections possible
    total = product(map(len,other_data))
    # Cycle trough and calculate a probability that each value in
    # "interest" has exactly "rank" greater than it

    for value in interest:
        # The probability that value from one list is better than current value
        greater = [count_greater(d, value) for d in other_data]
        # Handle the special case of rank 0  (no fancy operations necessary).
        if rank == 0:
            num_outcomes = product(greater)
        # Handle the special case of rank -1 (no fancy operations necessary).
        elif rank == len(other_data):
            num_outcomes = product([(len(d) - g) for (d,g) in zip(other_data,greater)])
        # Handle the general case, where we mix lessers and greaters.
        else:
            # Initialize the number greater to 0
            num_outcomes = 0
            # Get (in terms of indices) all combinations of the
            # other data sets that could be less than value
            for s in combinations(range(len(other_data)), rank):
                # Get the substitued "lesser" lengths of the selected indices
                lesser = [len(other_data[index]) - greater[index] for index in s]
                # Temporarily remove those elements from the list of lengths
                temp = [greater.pop(index) for index in s[::-1]][::-1]
                # Count up the number of possible unique sets
                num_outcomes += product(lesser + greater)
                # Reinsert the temporarily popped values back
                for i,index in enumerate(s):
                    greater.insert(index,temp[i])
        # Append the probability for this value to the list
        probabilities.append( num_outcomes / total )
    # Return the average of all probabilities (assuming each value is
    # equally likely to be selected from interest set)
    return sum(probabilities) / len(probabilities)

# Pre:  "interest" is an iterable of length > 0
#       "other_data" is a list of iterables each of length > 0
#       ** The minimum values at each index in "interest" and the
#          sub-sequences of "other_data" MUST NOT BE 0
# Post: The ratio defined in More and Wild in "Benchmarking
#       Derivative-Free Optimization Algorithms" (2009) as a
#       subsequent to performance profiles and data profiles.
def performance_ratio(interest, other_data):
    absolute_mins = [min(o) for o in zip(*(other_data+[interest]))]
    return [ i / m for i,m in zip(interest, absolute_mins) ]

# Pre:  "performance" is a list of performance ratios (numbers) where
#       smaller values denote better performances
#       "to_beat" is the performance ratio that we want to measure
#       against and guage percentage passed
# Post: The percentage of performances that were better than "to_beat"
def performance_profile(performance, to_beat):
    return sum(p <= to_beat for p in performance) / len(performance)

# Pre:  "performance" is a 2D python list of floats of performances
#       [iterations [trials]], where lower is better.
#       "best" is a 1D python list of floats of the absolute minimum
#       value that could be seen in trials for a given iteration
#       "to_beat" is a list of percentages that represent the maximum
#       distance from convergence allowed after n iterations for an
#       algorithm to be considered "passing"
# Post: The data profile for the performances given is returned in
#       terms of a dictionary as defined below
#       {max tolerance: [percentage passed per iteration]}
def data_profile(performances, best, to_beat=[10**(-i) for i in (1,3,5,7)]):
    profile = {t:[] for t in to_beat}
    for trials,abs_min in zip(performances,best):
        for t in to_beat:
            passed = sum( performances[0][i] - trials[i] >= 
                          (1 - t)*(performances[0][i] - abs_min) 
                          for i in range(len(trials)) ) / len(trials)
            profile[t].append( passed )
    return profile
