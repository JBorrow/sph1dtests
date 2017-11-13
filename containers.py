import gadget
import sph

class GadgetData(object):
    """
    Container for the GADGET data and required methods to calculate smoothed
    densities and pressures, based on the GADGET routines.
    """
    def __init__(self, positions, energies):
        """
        After creation, the properties can be accessed through:

        + GadgetData.smoothing_lengths,
        + GadgetData.densities,
        + GadgetData.pressures.
        """
        self.positions = positions
        self.energies = energies
        self.kernel = sph.kernel
        
        self.smoothing_lengths = self.calculate_smoothing_lengths(self.positions)
        self.densities = self.calculate_densities(self.positions, self.smoothing_lengths)
        self.pressures = self.calculate_pressures(self.densities, self.energies)


        return


    def separations(self, radii, radius):
        """
        Calculates the interparticle separations.
        """

        return [r - radius for r in radii]
    

    def calculate_densities(self, positions, smoothing_lengths):
        """
        Calculates the density at all of the particle positions.
        
        + h, the smoothing length at that position
        + positions, the positions of the other particles.
        """

        sep_between_all = [separations(r, positions) for r in positions]

        return map(gadget.density, sep_between_all, smoothing_lengths)


    def calculate_smoothing_lengths(self, positions, initial=1., eta=0.84, tol=0.01):
        """
        Calculates all of the smoothing lengths for all of the particles given
        in positions.

        Assumes they all have mass 1.
        """

        return map(gadget.h, [separations(r, positions) for r in positions])


    def calculate_pressures(self, densities, energies):
        """
        Calculates the pressures for all of the particles given their
        densities and energies.
        """

        return map(gadget.gas_pressure, densities, energies)

