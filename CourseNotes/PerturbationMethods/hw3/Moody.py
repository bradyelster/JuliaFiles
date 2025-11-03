import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Define functions for calculating friction factor
def colebrook(Re, eD):
    """Solves the Colebrook equation for friction factor using numerical methods."""
    def f_colebrook(f):
        return 1 / np.sqrt(f) + 2 * np.log10(eD / 3.7 + 2.51 / (Re * np.sqrt(f)))
    return fsolve(f_colebrook, 0.01)[0]

def laminar_friction_factor(Re):
    """Calculates the friction factor for laminar flow."""
    return 64 / Re

# Generate Data
Re_values = np.logspace(3, 8, 500)
eD_values = [0, 0.0001, 0.001, 0.01, 0.02, 0.04, 0.05]

for eD in eD_values:
    f_values = []
    for Re in Re_values:
        if Re < 2300:
            f_values.append(laminar_friction_factor(Re))
        else:
            f_values.append(colebrook(Re, eD))
    plt.loglog(Re_values, f_values, label=f'e/D = {eD}')

# Plot the diagram
plt.xlabel('Reynolds Number (Re)')
plt.ylabel('Friction Factor (f)')
plt.title('Moody Diagram')
plt.legend()
plt.grid(True, which='both', linestyle='--')
plt.show()