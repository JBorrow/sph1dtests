from sphtests import gadget, pressure_entropy, sph


class GadgetData(object):
    """
    Container for the GADGET data and required methods to calculate smoothed
    densities and pressures, based on the GADGET routines.
    """
    def __init__(self, positions, energies, eta=0.1, silent=False):
        """
        After creation, the properties can be accessed through:

        + GadgetData.smoothing_lengths,
        + GadgetData.densities,
        + GadgetData.pressures.
        """
        self.positions = positions
        self.energies = energies
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


    def calculate_smoothing_lengths(self, positions, initial=1., eta=0.2, tol=0.01):
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


class PressureEntropyData(object):
    """
    Container object for Pressure-Entropy.
    """
    def __init__(self, positions, energies, eta=0.1, silent=False):
        self.silent = silent

        if not silent: print("Grabbing the GadgetData object")
        self.gadget = GadgetData(positions, energies, eta, silent)

        if not silent: print("Starting Pressure-Entropy calculation")
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

        return


    def calculate_adiabats(self, energies, densities):
        """
        Calculates the ininitial Adiabats based on the traditional SPH
        calculation of the denstiy.
        """

        return list(map(pressure_entropy.A, energies, densities))


    def pressures(self, r, A, h):
        """
        Calculates the smoothed pressures according to Pressure-Entropy,
        at the positions of each of the particles.
        """

        sep_between_all = [sph.separations(r, positions) for r in positions]
        
        def p_given_A(separations, this_h):
            return pressure_entropy.pressure(separations, A, this_h)

        return list(
            map(
                p_given_A,
                sep_between_all,
                h
            )
        )


    def minimise_A(self, A, r, h, energies, tol=0.01):
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

        # As each particle's A depends on each other, we must iterate until
        # convergence in this lazy way.
        while difference > tol:
            new = old.copy()
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

