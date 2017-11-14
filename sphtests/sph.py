from numpy import pi

def kernel(r, h):
    """
    The standard GADGET Kernel, with 

    + r the interparticle separation
    + h the smoothing length of the particle.
    """
    factor = r/h
    factor2 = factor * factor
    prefactor = 1/(0.375 * h)

    if factor <= 0.5:
        poly = 1 - 6 * factor2 + 6 * factor2 * factor
    elif factor <= 1:
        one_minus_factor = 1 - factor
        poly = 2 * one_minus_factor * one_minus_factor * one_minus_factor
    else:
        poly = 0.

    return prefactor * poly


def separations(radius, radii):
    """
    Finds the separation between all in the radii list and the radius that is
    supplied.
    """
    return [abs(r - radius) for r in radii]


def diff(x, y):
    """
    Finds the sum of the absolute difference between x and y.
    """
    return sum([abs(i-j) for i, j in zip(x, y)])

