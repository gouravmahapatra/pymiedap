"""
Physical and spectroscopic constants used in CK-distribution computations.

All quantities are in CGS / HITRAN-compatible units unless otherwise noted.
"""

# Second radiation constant: h * c / k_B  [cm K]
C2 = 1.4387769

# Loschmidt number: number density of ideal gas at STP [molecules cm⁻³]
LOSCH = 2.6867774e19

# Boltzmann constant [J K⁻¹]
K_B = 1.380649e-23

# Avogadro constant [mol⁻¹]
N_A = 6.02214076e23

# HITRAN reference temperature [K]
T_REF = 296.0

# Standard atmosphere pressure [bar]
P_REF = 1.01325

# Conversion factor: bar → atm (HAPI expects pressure in atm)
BAR_TO_ATM = 0.986923

# Conversion factor: Pa → bar
PA_TO_BAR = 1.0e-5

# Gravity on Earth [m s⁻²]
G_EARTH = 9.807

# Gravity on Venus [m s⁻²]  (Ignatiev et al. value used in original code)
G_VENUS = 8.87

# Gravity on Mars [m s⁻²]
G_MARS = 3.721
