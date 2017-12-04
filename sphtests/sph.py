from numpy import pi, exp, sqrt


def gadget_kernel(r, h):
    """
    The standard GADGET Kernel, with 

    + r the interparticle separation
    + h the smoothing length of the particle.
    """
    factor = r/h
    factor2 = factor * factor
    prefactor = 4/(3 * h)

    if factor <= 0.5:
        poly = 1 - 6 * factor2 + 6 * factor2 * factor
    elif factor <= 1:
        one_minus_factor = 1 - factor
        poly = 2 * one_minus_factor * one_minus_factor * one_minus_factor
    else:
        poly = 0.

    return prefactor * poly


def cubic_kernel(r, h):
    """
    The cubic kernel from Price, 2012, with

    + r the interparticle separation
    + h the smoothing length of the particle
    """
    factor = r/h
    prefactor = 2/(3 * h)

    if factor < 1:
        poly = 0.25 * (2 - factor)**3 - (1 - factor)**3
    elif factor < 2:
        poly = 0.25 * (2 - factor)**3
    else:
        poly = 0.

    return prefactor * poly


def quintic_kernel(r, h):
    """
    The quntic kernel from Price, 2012, with

    + r the interparticle separation
    + h the smoothing length of the particle
    """
    q = r/h
    prefactor = 1/(120 * h)

    if q < 1:
        poly = (3 - q)**5 - 6 * (2 - q)**5 + 15 * (1 - q)**5
    elif q < 2:
        poly = (3 - q)**5 - 6 * (2 - q)**5
    elif q < 3:
        poly = (3 - q)**5
    else:
        poly = 0

    return poly * prefactor


def gaussian_kernel(r, h):
    """
    A perfect kernel.

    + r the interparticle separation
    + h the smoothing length of the particle
    """
    sigma = h/2
    prefactor = 1/sqrt(2 * pi * sigma)
    exponential = exp(-0.5 * (r / sigma)**2)

    return prefactor * exponential


def tophat_kernel(r, h):
    """
    A tophat kernel. Possibly the worst kernel you can have.

    + r the interparticle separation
    + h the smoothing length of the particle
    """
    prefactor = 1/(2 * h)

    if r/h < 1:
        return prefactor
    else:
        return 0


def triangle_kernel(r, h):
    """
    Triangle kernel.

    + r the interparticle separation
    + h the smoothing length of the particle
    """
    prefactor = 1/h

    if r/h < 1:
        return (1 - r/h) * prefactor
    else:
        return 0


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

