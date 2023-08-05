import numpy as np


L = round(786/1000,2) 
E = 190e9  
I = 6.12e-6  

# Applied point loads
loads = [

    {'magnitude': -4093.40, 'position': round(293/1000,2)},
    {'magnitude': 100 * 9.81, 'position': round(723/1000,2)}
]

# Bending moment equation due to point loads
def bending_moment(x):
    moment = 0
    for load in loads:
        distance = x - load['position']
        moment += load['magnitude'] * distance
    return moment

# curvature equation
def curvature(x):
    return bending_moment(x) / (E * I)

# integrate to deflection
def deflection(x):
    deflection_integral = np.zeros_like(x)
    for i in range(1, len(x)):
        deflection_integral[i] = deflection_integral[i-1] + np.trapz(curvature(x[:i+1]), x[:i+1])
    return deflection_integral

# calculate deflection
x_values = np.linspace(0, L, 100)
deflections = deflection(x_values)


for x, d in zip(x_values, deflections):
    print(f"At x = {x:.2f} m, Deflection = {d:.6f} m")
print(f"Maximum Deflection: {max(deflections):.6f} m")
