import numpy as np

def kernel(r, h):
    """
    The standard GADGET Kernel, with 

    + r the interparticle separation
    + h the smoothing length of the particle.
    """
    factor = r/h
    factor2 = factor * factor
    prefactor = 8/(np.pi * h * h * h)

    if factor <= 0.5:
        poly = 1 - 6 * factor2 + 6 * factor2 * factor
    elif factor <= 1:
        one_minus_factor = 1 - factor
        poly = 2 * one_minus_factor * one_minus_factor * one_minus_factor
    else:
        poly = 0.

    return prefactor * poly



