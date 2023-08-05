import numpy as np

# Beam properties
L = 10.0  # Length of the beam in meters
E = 200e9  # Young's modulus in N/m^2
I = 1e-6  # Moment of inertia in m^4

# Applied point loads
loads = [
    {'magnitude': 10.0, 'position': 2.0},
    {'magnitude': 20.0, 'position': 4.0},
    {'magnitude': 15.0, 'position': 7.0},
    {'magnitude': 25.0, 'position': 9.0}
]

# Bending moment equation due to point loads
def bending_moment(x):
    moment = 0
    for load in loads:
        distance = x - load['position']
        moment += load['magnitude'] * distance
    return moment

# Deflection equation using single integration
def deflection(x):
    return bending_moment(x) * x * (L - x) / (6 * E * I)

# Calculate deflection at different points along the beam
x_values = np.linspace(0, L, 100)
deflections = [deflection(x) for x in x_values]

# Print the calculated deflections
for x, d in zip(x_values, deflections):
    print(f"At x = {x:.2f} m, Deflection = {d:.6f} m")
print(max(deflections))