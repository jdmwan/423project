import math
import numpy as np
pi = math.pi
maxyield = 515e6/1.1

d4 = 24.4/2/1000 #mm
d5 = 24.4/2/1000 #mm
r6 = (84 + 2/3)/2/1000
# F_t = 183.0765
F_t = 24106.11
px = 11246.54
py = 4093.40
T = px * r6
Fr = 10133.28
g6dist = round(129.01/1000,2)
g7dist = round(599.01/1000,2)
crankDist = round(364.01/1000,2)
shaftL = round(728.02/1000,2)
distances = [round(94.01/1000,2),round(129.01/1000,2), round(164.01/1000,2), round(334.01/1000,2), round(364.01/1000,2),round(394.01/1000,2),round(564.01/1000,2),round(599.01/1000,2),round(634.01/1000,2),round(728.02/1000,2)]


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
Forcesy = [yr[0],0, Fr, 0,0,-py ,0,0, Fr, 0]
print(yr)
print(Mgear6y)
print(Mcranky)
print(Mgear7y)
print(Mbearingy)
Momenty = []
for i in range(len(distances)):
    if i == 0:
        Momenty.append(Forcesy[0] * distances[i])
    else:
        momentold = Momenty[i-1]
        forcessum = 0
        j=0
        while j <= i:
            forcessum += Forcesy[j]
            j +=1
        print(forcessum)
        momenttemp = forcessum * (distances[i] - distances [i-1]) + momentold
        Momenty.append(momenttemp)
print(Momenty)
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
Forcesx = [xr[0],0, F_t, 0,0,-px ,0,0, F_t, 0]
Momentx = []
for i in range(len(distances)):
    if i == 0:
        Momentx.append(Forcesx[0] * distances[i])
    else:
        momentold = Momentx[i-1]
        forcessum = 0
        j=0
        while j <= i:
            
            forcessum += Forcesx[j]
            j+=1
        momenttemp = forcessum * (distances[i] - distances [i-1]) + momentold
        Momentx.append(momenttemp)
print(Momentx)
M = math.sqrt((Momentx[5]**2) + (Momenty[5]**2))
print(M)
Msquared16 = (M**2) * (16 ** 2)
Tsquared16 = (T**2) * (16**2)
d = ((Msquared16 + Tsquared16)/((pi**2) *(maxyield**2))) ** (1/6)
print(d)

#fatigue
altx = 0
meanx = Momentx[5]
alty = 0
meany = Momenty[5]
ma = (meanx**2 + meany**2 )**(1/2)
Mm = (altx**2 + alty**2 )**(1/2)
Ta = T
Tm = 0
UTS = 675e6/1.1
ka = 0.739
kb = 0.5
kc = 1
kd = 1
ke = 0.5
Se = 1/2*UTS * ka *kb * kc * kd *ke
kt = 2.14
kts = 3.0
r = 1.75
n = 5e7
sqrta = (1.24 - 2.25e-3*UTS + 1.60e-6*(UTS**2) -4.11e-10 *(UTS**3))
sqrtashear = (0.958 - 1.83e-3 * UTS + 1.43e-6 * (UTS**2) -4.11E-10*(UTS**3))
kf = 1 + ((kt - 1)/(1+sqrta/(r**(1/2))))
kfs = 1 + ((kts - 1)/(1+sqrtashear/(r**(1/2))))
A = math.sqrt(4 * (kf * ma) ** 2 + 3 * (kfs * Ta) ** 2)
B = math.sqrt(4 * (kf * Mm) ** 2 + 3 * (kfs * Tm) ** 2)
d = ((16 * n / pi) * (A / Se + B / UTS)) ** (1/3)
print(d)