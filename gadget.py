import numpy as np
from scipy.optimize import root
from sph import kernel

def density(r, h, masses=None):
    """
    Calculates the SPH density from a set of particle separations.

    + r are the particle separations.
    + h is the smoothing lengths.
    + masses are the particle masses. If not given, we assume the
      particles are all equally massive with M=1.
    """
    
    def kernel_at_h(radius): return kernel(radius, h)

    weights = map(kernel_at_h, r)

    if masses is not None:
        density = sum([m*w for m, w in zip(masses, weights)])
    else:
        density = sum(weights)

    return density


def h(r, initial=1., mass=1, eta=0.84, tol=0.01, masses=None):
    """
    Calculates the smoothing length for a particle.

    + r are the interparticle separations
    + initial is the initial guess for the smoothing length
    + mass is the mass of the particle which is begin smoothed
    + eta is the smoothing length in terms of the mean interparticle separation
    + tol is the tolerance to find h to
    + masses are the masses of the other particles
    """
    def to_reduce(this_h):
        dens = density(r, this_h, masses)
        
        return (this_h / (mass * eta)) * dens - 1

    fitted = root(
        to_reduce,
        initial,
        tol=tol,
    )

    return fitted.x[0]


def gas_pressure(density, internal_energy, gamma=4./3.):
    """
    The gas pressure according to GADGET2, i.e.

            P = (gamma - 1) * rho * u

    + gamma has an initial value of 4/3
    """

    g_minus_1 = gamma - 1
    return g_minus_1 * density * internal_energy

