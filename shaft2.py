import math
import numpy as np
pi = math.pi
maxyield = 579e6

d4 = 24.4/2/1000 #mm
d5 = 24.4/2/1000 #mm
r6 = (84 + 2/3)/2/1000
# F_t = 183.0765
F_t = 183.0765/2*1000
px = 52.85*1000
py = 24.02*1000
T = px * r6
Fr = 76.94/2*1000
g6dist = round(20/1000,2)
g7dist = round(280/1000,2)
crankDist = round(230/1000,2)
shaftL = round(300/1000,2)



px = 52.85*1000
py = 24.02*1000
print(px)
# reactions y 
sumFy = Fr +Fr - py
sumMy = g6dist* Fr + (g7dist*Fr) - (crankDist*py)
by = np.array([sumFy,sumMy])
ay= np.array([[-1,-1] , [0, -shaftL]])
yr = np.linalg.solve(ay,by)
Mgear6y = g6dist * yr[0]
Mcranky = (crankDist -g6dist) * (yr[0] + Fr) + Mgear6y
Mgear7y = (g7dist -crankDist ) * (yr[0] + Fr -py) + Mcranky
Mbearingy = (shaftL -g7dist) * (yr[0] + Fr -py +Fr) +Mgear7y
print(yr)
print(Mgear6y)
print(Mcranky)
print(Mgear7y)
print(Mbearingy)
#reactions x
sumFx = F_t + F_t - px
sumMx = g6dist* F_t + (g7dist*F_t) - (crankDist*px)
bx = np.array([sumFx,sumMx])
ax = np.array([[-1,-1] , [0, -shaftL]])
xr = np.linalg.solve(ax,bx)
print(xr)
Mgear6x = g6dist * xr[0]
Mcrankx = (crankDist -g6dist) * (xr[0] + F_t) + Mgear6x
Mgear7x = (g7dist -crankDist ) * (xr[0] + F_t -px) + Mcrankx
Mbearingx = (shaftL -g7dist) * (xr[0] + F_t -px +F_t) +Mgear7x
print(Mgear6x)
print(Mcrankx)
print(Mgear7x)
print(Mbearingx)

M = math.sqrt((Mcrankx**2) + (Mcranky**2))
print(M)
Msquared16 = (M**2) * (16 ** 2)
Tsquared16 = (T**2) * (16**2)
d = ((Msquared16 + Tsquared16)/((pi**2) *(maxyield**2))) ** (1/6)
print(d)