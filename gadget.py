import numpy as np
from sph import kernel

def density(r, h, masses=None):
    """
    Calculates the SPH density from a set of particle separations.

    + r are the particle separations.
    + h are the smoothing lengths.
    + masses are the particle masses. If not given, we assume the
      particles are all equally massive with M=1.
    """

    weights = map(kernel, r, h)
    if masses is not None:
        density = sum([m*w for m, w in masses, weights])
    else:
        density = sum(weights)

    return density


