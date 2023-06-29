import math
import numpy as np
pi = math.pi
maxyield = 1000e6
h = 4/1000
r = 60/1000
Lm = 437400
d6 = 122/2 #mm
d7 = 122 #mm
# F_t = 183.0765
F_t = 183.0765/2*1000
Torque = F_t *d6 /1000
# print(Torque)
# Lm = 300000 old
Fr = 76.94/2*1000
g6dist = round(20/1000,2)
g7dist = round(580/1000,2)
crankDist = round(200/1000,2)
crankDist2 = round(400/1000,2)
shaftL = round(600/1000,2)
alpha = (math.acos(1 - (h/r)))

beta = (1/5*math.sin(alpha))
# print(alpha)git
# print(beta)
P = Lm / math.cos(beta)
T = P* math.sin(alpha+beta)
# print(P, " ", T)
Tcrank = r * T
# print(Tcrank)
px = P * math.sin(beta)
py = P * math.cos(beta)
# print(px)
# reactions y 
sumFy = Fr +Fr - py
sumMy = g6dist* Fr + (g7dist*Fr) - (crankDist*1/2*py) - (crankDist2*1/2*py)
by = np.array([sumFy,sumMy])
ay= np.array([[-1,-1] , [0, -shaftL]])
yr = np.linalg.solve(ay,by)
Mgear6y = g6dist * yr[0]
Mcranky = (crankDist -g6dist) * (yr[0] + Fr) + Mgear6y
Mcranky2 = (crankDist2 - crankDist) * (yr[0] + Fr - 1/2*py)+Mcranky
Mgear7y = (g7dist -crankDist2 ) * (yr[0] + Fr -py) + Mcranky2
Mbearingy = (shaftL -g7dist) * (yr[0] + Fr -py +Fr) +Mgear7y
print(yr)
# print(Mgear6y)
print(Mcranky)
print(Mcranky2)
print(Mgear7y)
print(Mbearingy)
#reactions x
sumFx = F_t + F_t - px
sumMx = g6dist* F_t + (g7dist*F_t) - (crankDist*1/2*px) - (crankDist2*1/2*px)
bx = np.array([sumFx,sumMx])
ax = np.array([[-1,-1] , [0, -shaftL]])
xr = np.linalg.solve(ax,bx)
print(xr)
Mgear6x = g6dist * xr[0]
Mcrankx = (crankDist -g6dist) * (xr[0] + F_t) + Mgear6x
Mcrankx2 = (crankDist2 - crankDist) *(xr[0]+ F_t - 1/2*px) + Mcrankx
Mgear7x = (g7dist -crankDist2 ) * (xr[0] + F_t -px) + Mcrankx2
Mbearingx = (shaftL -g7dist) * (xr[0] + F_t -px +F_t) +Mgear7x
# print(Mgear6x)
print(Mcrankx)
print(Mcrankx2)
print(Mgear7x)
print(Mbearingx)

M = math.sqrt((Mcrankx**2) + (Mcranky**2))
print(M)
Msquared16 = (M**2) * (16 ** 2)
Tsquared16 = (T**2) * (16**2)
d = ((Msquared16 + Tsquared16)/((pi**2) *(maxyield**2))) ** (1/6)
print(d)