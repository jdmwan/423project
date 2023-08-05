import math
import numpy as np
pi = math.pi
maxyield = 515e6/1.1

d4 = 24.4/2/1000 #mm
d5 = 24.4/2/1000 #mm
r6 = (84 + 2/3)/2/1000
# F_t = 183.0765

px = 1883.01
py = 685.36
T = 559.324
Fr = 100 * 9.81
g6dist = round(293/1000,2)

crankDist = round(723/1000,2)
shaftL = round(786/1000,2)

distances = [round(263/1000,2),round(293/1000,2),round(323/1000,2),round(708/1000,2),round(723/1000,2),round(738/1000,2),round(786/1000,2)]

# px = 52.85*1000
# py = 19.24*1000
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
Forcesy = [yr[0], 0, Fr, 0,0, -py,0]
print(yr)
print(Mgear6y)
print(Mcranky)

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
Forcesx = [xr[0], 0, 0, 0,0, -px,0]
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
M = math.sqrt((Momentx[2]**2) + (Momenty[2]**2))
print(M)
Msquared16 = (M**2) * (16 ** 2)
Tsquared16 = (T**2) * (16**2)
d = ((Msquared16 + Tsquared16)/((pi**2) *(maxyield**2))) ** (1/6)
print(d)

#fatigue
altx = 0
meanx = Momentx[2]
alty = 0
meany = Momenty[2]
ma = (meanx**2 + meany**2 )**(1/2)
Mm = (altx**2 + alty**2 )**(1/2)
Ta = T
Tm = 0
UTS = 675e6/1.1
ka = 1
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