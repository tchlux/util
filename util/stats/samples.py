from util.stats import *

# =====================================================
#      Approximating the desired number of samples     
# =====================================================


# Compute the confidence of "num_samples" samples having less than or
# equal to "max_error" at exactly 1/2. This is achieved by factoring
# the expression (sum_{i=j}^k choose(n,i) / 2^n) to reduce operations.
def _half_confidence(num_samples, max_error):
    from util.math import Fraction, choose
    # Calculate the maximum number of steps away from the mid-point
    # that we can take *before* violating the error condition.
    odd = num_samples % 2
    start = (num_samples+1) // 2
    #         initial_error + steps * step_error <= max_error
    #   (max_error - initial_error) / step_error >= steps
    #  (max_error - initial_error) * num_samples >= steps
    # knowing that                    step_error  = 1 / num_samples
    min_error = odd * Fraction(1, (2*num_samples))    
    steps = int((max_error - min_error) * num_samples)
    # Put steps into allowable bounds.
    steps = max(0, min(num_samples - start, steps))
    # Handle two cases where there is no functioning confidence bound.
    if (odd) and (max_error < min_error): return Fraction()
    # Compute the fraction.
    numerator   = 1
    denominator = 2**(num_samples)
    # First compute the constant multiple outside of the sum.
    for i in range(start+steps+1, num_samples+1): numerator   *= i
    for i in range(2,             start+1):       denominator *= i
    # Compute the sum of the inner parts of the distribution.
    total = 0
    for i in range(start, start+steps+1):
        v = 1
        for j in range(i+1,             start+steps+1): v *= j
        for j in range(1+num_samples-i, start+1):       v *= j
        if (v > 1): total += v * (2 if ((i != start) or odd) else 1)
    # Perform the final numerator update.
    if (total > 0): numerator *= total
    # Return the (reduced) fraction.
    return Fraction(numerator, denominator)


# Given a list of numbers, return True if the given values provide an
# estimate to the underlying distribution with confidence bounded error.
def samples(size=None, error=None, confidence=None, at=None):
    # Determine what to calculate based on what was provided.
    from util.math import choose, Fraction
    if   (size is None):       to_calculate = "samples"
    elif (error is None):      to_calculate = "error"
    elif (confidence is None): to_calculate = "confidence"
    else:                      to_calculate = "verify"
    # Default evaluation point is at (1/2), where the error is greatest.
    if (at is None): at = Fraction(1, 2)
    else:            at = Fraction(at)
    # Set the default values for other things that were not provided.
    if type(error) == type(None):      error      = Fraction(10, 100)
    if type(confidence) == type(None): confidence = Fraction(95, 100)
    # Convert error and confidence to fraction types if necessary.
    if not type(error)      == Fraction: error      = Fraction(error)
    if not type(confidence) == Fraction: confidence = Fraction(confidence)
    # If the user provided something with a length, use that number.
    if hasattr(size, "__len__"): size = len(size)
    # \sum_{i=0}^n choose(n, i) * ( at^i (1-at)^(n-i) )
    if (size is not None):
        # Compute the probability of any given observed EDF value.
        prob = lambda i: choose(size, i) * (at**i * (1-at)**(size-i))
        # If we are calculating the confidence or verifying, compute confidence.
        if to_calculate in {"confidence", "verify"}:
            if (at == 1/2): conf = _half_confidence(size, error)
            else:
                conf = Fraction()
                steps = 0
                # Sum those probabilities that are closer than "error" distance.
                for i in range(size+1):
                    p = Fraction(i, size)
                    if (abs(p - at) <= error):
                        steps += 1
                        conf += prob(i)
            # Return the total confidence.
            if to_calculate == "confidence": return float(conf)
            else:                            return conf >= confidence
        elif to_calculate == "error":
            # Store the "contained" outcomes by "allowed error".
            error = Fraction()
            contained = Fraction()
            # Sort the percentiles by their distance from "at".
            i_p = sorted(enumerate(Fraction(i,size,reduce=False)
                                   for i in range(size+1)),
                         key=lambda ip: abs(ip[1]-at))
            # Cycle through percentiles, starting closest to "at" and moving out.
            for step in range(len(i_p)):
                # If this step has the same probability as the last, skip.
                if (i_p[step][1] == i_p[step-1][1]): continue
                i, p = i_p[step]
                # Compute the amount of data contained by this step away.
                next_contained = contained + prob(i)
                # If the distance from "at" is the same for two steps, take two.
                if (step+1 < len(i_p)) and (abs(at-i_p[step][1]) == abs(at-i_p[step+1][1])):
                    next_contained += prob( i_p[step+1][0] )
                # Only update the "allowed error" if confidence is maintained.
                if next_contained < confidence:
                    contained = next_contained
                    error = abs(i_p[step][1] - at)
                else: break
            return float(error)
    else:
        # Compute the number of samples required.
        size, step = 2**10, 2**9
        # print("Desired ----------------")
        # print("error:      ",error)
        # print("confidence: ",confidence)
        # for size in range(2, 500):
        #     conf_below = samples(size-1, error=error, at=at)
        #     conf_at = samples(size, error=error, at=at)
        #     print("", "size: ",size, float(f"{conf_below:.2e}"), float(f"{conf_at:.2e}"))
        # exit()

        under, over = 0, None
        # We have the right size when any smaller size is not passable.
        conf_below = samples(size=size-1, error=error, at=at)
        conf_at = samples(size=size, error=error, at=at)
        print("", "size: ",size, float(f"{conf_below:.2e}"), float(f"{conf_at:.2e}"))
        while not (conf_below < confidence <= conf_at):
            if conf_at < confidence:
                # Update "under". Scale up if we haven't found "over".
                under = max(under, size)
                if (over == None): step *= 2
                # Take the step.
                size += step
            else:
                # Update "over". Take step. Scale down.
                over = min(over if (over != None) else float('inf'), size)
                size = size - step
                step = step // 2
            # Recompute the confidence at and below this step size.
            conf_below = samples(size-1, error=error, at=at)
            conf_at = samples(size, error=error, at=at)
            print("", "size: ",size, float(f"{conf_below:.2e}"), float(f"{conf_at:.2e}"))
            # Correct for strange sample size error that can happen bc
            # of alignment of "at" and one of the sample values.
            if conf_at < conf_below:
                size = size-1
                conf_at = conf_below
                conf_below = samples(size-1, error=error, at=at)
                print("", "size: ",size, float(f"{conf_below:.2e}"), float(f"{conf_at:.2e}"))
        # Return the computed best sample size.
        return size


