import numpy as np

# Beam properties
L = round(738/1000, 2)  # Length of the beam
E = 190e9  # Young's modulus
I = 4.22e-5  # Moment of inertia (example value)

# Applied point loads
loads = [
    {'magnitude': 10133.28, 'position': round(134/1000, 2)},
    {'magnitude': -11183.66, 'position': round(269/1000, 2)},
    {'magnitude': -11183.66, 'position': round(469/1000, 2)},
    {'magnitude': 10133.28, 'position': round(604/1000, 2)}
]

# Bending moment equation due to point loads
def bending_moment(x):
    moment = 0
    for load in loads:
        distance = x - load['position']
        moment += load['magnitude'] * distance
    return moment

# Curvature equation
def curvature(x):
    return bending_moment(x) / (E * I)

# Deflection equation using cumulative integration
def deflection(x):
    deflection_integral = np.zeros_like(x)
    for i in range(1, len(x)):
        deflection_integral[i] = deflection_integral[i-1] + np.trapz(curvature(x[:i+1]), x[:i+1])
    return deflection_integral

# Calculate deflection at different points along the beam
x_values = np.linspace(0, L, 100)
deflections = deflection(x_values)

# Print the calculated deflections
for x, d in zip(x_values, deflections):
    print(f"At x = {x:.2f} m, Deflection = {d:.6f} m")
print(f"Maximum Deflection: {max(deflections):.6f} m")