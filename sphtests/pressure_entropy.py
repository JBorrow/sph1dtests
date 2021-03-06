from sphtests.sph import gadget_kernel
from scipy.optimize import root


def A(energy, density, gamma=4./3.):
    """
    The adiabat for the particle.

    + energy is the internal energy of the particle
    + density is the density of the particle
    """
    g_minus_1 = gamma - 1.
    return energy * g_minus_1 / (density**g_minus_1) 


def pressure(r, A, h, kernel, gamma=4./3., masses=None):
    """
    Calculates the pressure-entropy smoothed pressure (note this is _not_ the
    gas pressure) based on the:

    + r, the interparticle separations,
    + A, the adiabats of the particles,
    + h, the smoothing length of the particle being considered,
    + gamma of the gas,
    + masses, the masses of the other particles. If none, assumed to all be 1,
    + kernel, a callable with arguents (r, h). Defaults to GADGET.
    """
    
    one_over_gamma = 1./gamma

    def kernel_at_h(radius): return kernel(radius, h)
    def A_one_over_gamma(A): return A**(one_over_gamma)

    weights = map(kernel_at_h, r)
    As = map(A_one_over_gamma, A)

    if masses is not None:
        P = sum([m*w*a for m, w, a in zip(masses, weights, As)])
    else:
        P = sum([w*a for w, a in zip(weights, As)])

    return P**gamma


def smoothed_density(A, P, gamma=4./3.):
    """
    The smoothed density of the particle, according to Pressure-Entropy
    SPH. This is _not_ the physical density.

    + A is the adiabat of the particle
    + P is the pressure of the particle
    + gamma of the gas.
    """
    return (P/A)**(1/gamma)


def A_reduced(r, A, h, energy, initial, index, kernel, gamma=4./3., tol=None, masses=None):
    """
    Calculates the optimum value of A for the particle, with

    + r the interparticle sepataions,
    + A the adiabats of all of the relevant particles,
    + h the smoothing length of the particle,
    + energy the internal energy of the particle,
    + initial, the initial value of A,
    + index, the index in A that corresponds to _this_ particle.
      If the particle is not in the array, set index to be -1.
    + gamma of the gas,
    + tol the tolerance to find A to,
    + masses are the masses of all of the particles. Assumed to be 1.
    """
    gamma_minus_1 = gamma - 1

    def to_reduce(this_A):
        # Has some lovely side effects, sorry about that...
        if index != -1:
            A[index] = this_A

        P_sum = pressure(r, A, h, kernel, gamma, masses)
        P_calc = this_A * ((energy/this_A) * gamma_minus_1)**(gamma/gamma_minus_1)

        return 1 - P_sum/P_calc

    if tol is not None:
        fitted = root(
            to_reduce,
            initial,
            tol=tol,
        )
    else:
        fitted = root(
            to_reduce,
            initial,
        )

    if index != -1:
        A[index] = fitted.x[0]

    return fitted.x[0], A

