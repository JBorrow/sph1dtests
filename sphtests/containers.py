from sphtests import gadget, pressure_entropy, sph


class GadgetData(object):
    """
    Container for the GADGET data and required methods to calculate smoothed
    densities and pressures, based on the GADGET routines.
    """
    def __init__(
            self,
            positions,
            energies=None,
            adiabats=None,
            eta=5,
            silent=False,
            gamma=4./3.
        ):
        """
        Please specify one of _either_ energies or adiabats.

        After creation, the properties can be accessed through:

        + GadgetData.smoothing_lengths,
        + GadgetData.densities,
        + GadgetData.pressures.
        """

        if adiabats is None and energies is None:
            raise AttributeError(
                "You have not provided adiabats or internal energies."
            )
        elif adiabats is not None and energies is not None:
            raise AttributeError(
                "You have provided _both_ adiabats and energies, please " +\
                "Only provide one."
            )

        self.positions = positions
        self.energies = energies
        self.adiabats = adiabats
        self.eta = eta
        self.kernel = sph.kernel
        
        if not silent: print("Calculating smoothing lengths")
        self.smoothing_lengths = self.calculate_smoothing_lengths(
            self.positions,
            eta=self.eta
        )

        if not silent: print("Calculating densities")
        self.densities = self.calculate_densities(
            self.positions,
            self.smoothing_lengths
        )

        if not silent: print("Calculating pressures")
        if adiabats is not None:
            self.pressures = self.calculate_pressures_adiabats(
                self.densities,
                self.adiabats
            )

            # We now need to update internal energies to match.
            self.energies = self.calculate_energies(
                self.adiabats,
                self.densities
            )
            
        else: # We must be using internal energies
            self.pressures = self.calculate_pressures(
                self.densities,
                self.energies
            )


        return



    def calculate_densities(self, positions, smoothing_lengths):
        """
        Calculates the density at all of the particle positions.
        
        + h, the smoothing length at that position
        + positions, the positions of the other particles.
        """

        sep_between_all = [sph.separations(r, positions) for r in positions]

        return list(
            map(
                gadget.density,
                sep_between_all,
                smoothing_lengths
            )
        )


    def calculate_smoothing_lengths(self, positions, initial=1., eta=0.2, tol=None):
        """
        Calculates all of the smoothing lengths for all of the particles given
        in positions.

        Assumes they all have mass 1.
        """

        def this_h(sep): return gadget.h(sep, eta=eta, tol=tol)

        return list(
            map(
                this_h,
                [sph.separations(r, positions) for r in positions]
            )
        )


    def calculate_pressures(self, densities, energies):
        """
        Calculates the pressures for all of the particles given their
        densities and energies.
        """

        return list(map(gadget.gas_pressure, densities, energies))


    def calculate_pressures_adiabats(self, densities, adiabats, gamma=4./3.):
        """
        Calculates the pressures for all of the particles given their
        densities and adiabats.
        """

        return list(map(gadget.gas_pressure_adiabat, densities, adiabats))


    def calculate_energies(self, adiabats, densities, gamma=4./3.):
        """
        Calculates the internal energies of the particles given their
        densities and adiabats using the thermodynamical relation:
        
                u = A rho^(gamma - 1)/(gamma - 1)

        """
        def u_at_gamma(a, rho): return gadget.internal_energy(a, rho, gamma)

        return list(
            map(
                u_at_gamma,
                adiabats,
                densities
            )
        )


class PressureEntropyData(object):
    """
    Container object for Pressure-Entropy.
    """
    def __init__(
            self,
            positions,
            energies=None,
            adiabats=None,
            eta=5,
            silent=False,
            gamma=4./3.
        ):
        """
        + positions are the positions of the particles,
        + energies/adiabats:
            You must provide either of the following:
                - energies are the internal energies,
                - adiabats are the initial adiabats,
        + eta is the kernel-eta,
        + silent is a boolean, if True it enables printing of information.
        + gamma is the ratio of specific heats of the gas.

        This class makes available the following properties:
        + positions,
        + energies,
        + gadget, the GadgetData object that is associated with 
          the data you have input,
        + densities, the TSPH (Traditional SPH) densities of the particles,
        + smoothing_lenghts, the TSPH smoothing lengths of the particles,
        + adiabats, the PESPH (Pressure-Entropy SPH) minimised adiabats
          associated with each of the particles,
        + smootehd_pressures, the PESPH smoothed pressures of the particles,
        + smoothed_densities, the PESPH smoothed densities of the particles.
        """
        
        if adiabats is None and energies is None:
            raise AttributeError(
                "You have not provided adiabats or internal energies."
            )
        elif adiabats is not None and energies is not None:
            raise AttributeError(
                "You have provided _both_ adiabats and energies, please " +\
                "Only provide one."
            )

        self.silent = silent
        self.positions = positions

        if not silent: print("Grabbing the GadgetData object")
        self.gadget = GadgetData(
            positions,
            energies,
            adiabats,
            eta,
            silent,
            gamma
        )

        self.densities = self.gadget.densities
        self.smoothing_lengths = self.gadget.smoothing_lengths
        self.energies = self.gadget.energies
        self.adiabats = self.gadget.adiabats

        if not silent: print("Starting Pressure-Entropy calculation")
        if adiabats is None:
            # Sort out your adiabats/energies.
            self.adiabats = self.calculate_adiabats(
                self.gadget.energies,
                self.gadget.densities,
            )

        if not silent: print("Minimising to find values of A")
        self.adiabats = self.minimise_A(
            self.adiabats,
            self.gadget.positions,
            self.gadget.smoothing_lengths,
            self.gadget.energies
        )

        if not silent: print("Calculating smoothed pressures")
        self.smoothed_pressures = self.pressures(
            self.gadget.positions,
            self.adiabats,
            self.gadget.smoothing_lengths,
        )

        if not silent: print("Calculating smoothed densities")
        self.smoothed_densities = self.density_twiddle(
            self.adiabats,
            self.smoothed_pressures
        )

        return


    def calculate_adiabats(self, energies, densities):
        """
        Calculates the ininitial Adiabats based on the traditional SPH
        calculation of the denstiy.

        + energies are the internal energies
        + densities are the densities of the particles.
        """

        return list(map(pressure_entropy.A, energies, densities))


    def pressures(self, r, A, h):
        """
        Calculates the smoothed pressures according to Pressure-Entropy,
        at the positions of each of the particles.

        + r are the particle positions,
        + A are the particle Adiabats,
        + h are the smoothing lengths of the particles from GADGETSPH.
        """

        sep_between_all = [sph.separations(this_r, r) for this_r in r]
        
        def p_given_A(separations, this_h):
            return pressure_entropy.pressure(separations, A, this_h)

        return list(
            map(
                p_given_A,
                sep_between_all,
                h
            )
        )


    def density_twiddle(self, A, P, gamma=4./3.):
        """
        The smoothed density of the particles.

        + A are the adiabats,
        + P are the pressures,
        + gamma of the gas.
        """
        def smoothed_at_gamma(a, p): return pressure_entropy.smoothed_density(a, p, gamma)
        return list(
            map(
                smoothed_at_gamma, A, P
            )
        )


    def minimise_A(self, A, r, h, energies, tol=1e-7):
        """
        Finds the equlibrium value of A using the Pressure-Entropy SPH
        technique.

        + A are the adiabats for each particle
        + r are the positions of each particle
        + h are the smoothing lenghts of each particle
        + energies are the internal energies of each particles
        + tol is the tolerance between iterations.
        """

        difference = tol + 1
        old = A.copy()
        new = old.copy()

        # As each particle's A depends on each other, we must iterate until
        # convergence in this lazy way.
        while difference > tol:
            # We iterate over each particle and update its A to be the
            # Equilibrium given the values of its neighbors
            for index, (this_A, this_h, this_u) in enumerate(zip(old, h, energies)):
                separations = sph.separations(r[index], r)

                this_A, new = pressure_entropy.A_reduced(
                    separations, new, this_h, this_u, this_A, index
                )
                
            difference = sph.diff(old, new)
            old = new.copy()

            if not self.silent: print("Difference: {}".format(difference))

        return old

