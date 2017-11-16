SPH Tests
=========

Several SPH tests that look at the diferences between tradiational and 'Pressure-Entropy' SPH.

This module and the notebooks included require `python3`.

To use the API objects, you can do the following:

```python
from sphtests import PressureEntropyData

data = PressureEntropyData(
    positions, # A list of positions in 1D,
    energies=E, # A list of internal energies of your particles,
    eta=eta, # The kernel eta (for the GADGET2 Kernel -- try 50. Default is 5.)
)

# You can then access the TSPH properties as follows:

tsph_data = data.gadget

# This object contains:
densities = tsph.densities
pressures = tsph.pressures
hsml = tsph.smoothing_lengths

# Then, the Pressure Entropy stuff can be accessed as follows:
adiabats = data.adiabats
smoothed_P = data.smoothed_pressures
smoothed_rho = data.smoothed_densities

# Note that you can also specify adiabats rather than energy as follows:

data_adiabat = PressureEntropyData(
    positions,
    adiabats=A,
    eta=eta,
)
```

Users can also choose their own Kernels. There are two Kernels available with
this module; the one used in GADGET-2 (re-normalized for 1D use), or the C_2
kernel used in ANARCHY. They can be used as follows:

```python
import sphtests.sph as sph

data = PressureEntropyData(
    ...
    kernel=sph.gadget_kernel,
)

# or

data = PressureEntropyData(
    ...
    kernel=sph.ANARCHY_kernel,
)
```

Users can also define their own kernels. Note that their functions should have
units of 1/length, and have the following structure:

```python
def my_kernel(r, h):
    """
    + r is the interparticle separation,
    + h is the smoothing length,
    """
    return ...
```
