import random


def prime_test(N, k):
    # This is main function, that is connected to the Test button. You don't need to touch it.
    return fermat(N, k), miller_rabin(N, k)


# Time complexity of n^3
def mod_exp(x, y, N):
    # When we are finished with recursion
    if y == 0:
        return 1
    z = mod_exp(x, y//2, N)

    # Begin returning values from the bottom up until we get our final value
    if y % 2 == 0:
        return z ** 2 % N
    else:
        return x * z ** 2 % N


# Time complexity of n^2 with integer k which is size n bits
def fprobability(k):
    return 1 - (1 / 2 ** k)


# Time complexity of n^2 with integer k which is size n bits
def mprobability(k):
    return 1 - (1 / 4 ** k)


# Time complexity of kn^3 because we call mod_exp (which is order n^3) k times
def fermat(N, k):

    # if N is an even number return composite
    if N % 2 == 0:
        return 'composite'

    # Run through the number of tests
    while k > 0:

        # select a random base and compute modular exponentiation
        rand_num = random.randint(1, N - 1)
        num = mod_exp(rand_num, N - 1, N)

        # If the result of the modular exponentiation is not 1- we know the number is composite
        if num != 1:
            return 'composite'
        k -= 1

    # Return prime if we didn't return composite
    return 'prime'


# Time complexity of kn^2 because we call mod_exp (which is order n^3) k times. A
def miller_rabin(N, k):

    # If number is even we know it is composite
    if N % 2 == 0:
        return 'composite'

    # Set exponent equal to test value minus one
    exponent = N - 1

    while k > 0:

        # Select random base, compute modular exponentiation to see if it passes the first test
        rand_num = random.randint(1, N - 1)
        num = mod_exp(rand_num, exponent, N)
        if num == 1:

            # if N passes the first test call helper function to run miller rabin test until it returns
            encounter = False
            is_prime = miller_rabin_helper(rand_num, N - 1, N, encounter)
            if not is_prime:
                return 'composite'
        else:
            return 'composite'
        k -= 1

    return 'prime'


# n^3 apply max rule would be n^3 + (n-1)/2
def miller_rabin_helper(rand_num, exponent, test_num, encounter):

    # If the exponent is odd we want to move on to test the next base
    if exponent % 2 != 0:
        return True

    result = mod_exp(rand_num, exponent, test_num)

    if result == test_num - 1:
        encounter = True

    # if the result of mod_exp is not one and it is also not test_num - 1 return false if we have not seen a result
    # be test_num - 1
    if result != 1 and result != test_num - 1:
        if not encounter:
            return False
        else:
            return True

    else:
        return miller_rabin_helper(rand_num, exponent/2, test_num, encounter)

