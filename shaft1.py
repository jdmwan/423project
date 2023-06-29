import math
import numpy as np
pi = math.pi
maxyield = 579e6

d4 = 24.4/2/1000 #mm
d5 = 24.4/2/1000 #mm
r6 = (84 + 2/3)/2/1000
# F_t = 183.0765

px = 52.85*1000
py = 24.02*1000
T = 559.324
Fr = 100 * 9.81
g6dist = round(100/1000,2)

crankDist = round(730/1000,2)
shaftL = round(800/1000,2)



px = 52.85*1000
py = 19.24*1000
print(px)
# reactions y 
sumFy = Fr - py
sumMy = g6dist* Fr - (crankDist*py)
by = np.array([sumFy,sumMy])
ay= np.array([[-1,-1] , [0, -shaftL]])
yr = np.linalg.solve(ay,by)
Mgear6y = g6dist * yr[0]
Mcranky = (crankDist -g6dist) * (yr[0] + Fr) + Mgear6y
Mbearingy = (shaftL -crankDist) * (yr[0]-py +Fr) +Mcranky
print(yr)
print(Mgear6y)
print(Mcranky)

print(Mbearingy)
#reactions x
sumFx = -px
sumMx = - (crankDist*px)
bx = np.array([sumFx,sumMx])
ax = np.array([[-1,-1] , [0, -shaftL]])
xr = np.linalg.solve(ax,bx)
print(xr)
Mgear6x = g6dist * xr[0]
Mcrankx = (crankDist -g6dist) * (xr[0]) + Mgear6x
Mbearingx = (shaftL -crankDist) * (xr[0] -px) +Mcrankx
print(Mgear6x)
print(Mcrankx)

print(Mbearingx)

M = math.sqrt((Mcrankx**2) + (Mcranky**2))
print(M)
Msquared16 = (M**2) * (16 ** 2)
Tsquared16 = (T**2) * (16**2)
d = ((Msquared16 + Tsquared16)/((pi**2) *(maxyield**2))) ** (1/6)
print(d)