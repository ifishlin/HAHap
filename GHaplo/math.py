__author__ = 'ifish'

from math import factorial, log, ceil, pow, exp

def combination(n,k):
    numerator=factorial(n)
    denominator=(factorial(k)*factorial(n-k))
    answer=numerator/denominator
    return answer

def likelihood(n, n1, n2):
    mp1 = 0.49
    mp2 = 0.49
    mp3 = 0.02
    n3 = n - n1 - n2
    p = n1 * log(mp1, 2)  + n2 * log(mp2, 2)  + n3 * log(mp3, 2)   + log(combination(n, n1),2)  + log(combination(n - n1, n2), 2)
    return p

def prior(ref, h1, h2):
    h11_p = 0.98 if h1[0] == ref[0] else 0.02
    h12_p = 0.98 if h1[1] == ref[1] else 0.02
    h21_p = 0.98 if h2[0] == ref[0] else 0.02
    h22_p = 0.98 if h2[1] == ref[1] else 0.02
    return log(h11_p * h12_p * h21_p * h22_p, 2)

def calc_posterior(n, n1, n2, ref, h1, h2):
    return likelihood(n ,n1 ,n2) + prior(ref, h1, h2)

def score(n, n1, n2):
    mp1 = 0.49
    mp2 = 0.49
    mp3 = 0.02
    n3 = n - n1 -n2
    return n1 * log(mp1, 2) + n2 * log(mp2, 2) + n3 * log(mp3, 2) + log(combination(n, n1),2) + log(combination(n - n1, n2), 2)

def normalized_factor(n):
    mp1 = 0.49
    mp2 = 0.49
    mp3 = 0.02
    n3 = round(n * mp3)
    if((n-n3)%2 == 0):
        n1 = (n - n3)/2
        n2 = n1
    else:
        n1 = ceil((n - n3)/2)
        n2 = n - n3 - n1
    return n1 * log(mp1, 2) + n2 * log(mp2, 2) + n3 * log(mp3, 2) + log(combination(n, n1), 2) + log(combination(n - n1, n2), 2)

def sigmoid(x, alpha):
    return 1 / (1 + exp(alpha*(-x + 0.5)))

def normalized_score_scale_to_max(n, n1, n2, orignal_n, max_sample_size, alpha):
    """scale should put it here"""
    p = score(n, n1, n2) - normalized_factor(n)
    return log(sigmoid(orignal_n/max_sample_size, alpha) * pow(2, p), 2)

