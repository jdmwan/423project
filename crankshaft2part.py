import math
import numpy as np
pi = math.pi
maxyield = 515e6/1.1
h = 4/1000
r = 60/1000
Lm = 437400
d6 = 122/2 #mm
d7 = 122 #mm
# F_t = 183.0765
F_t = 24106.11
Torque = F_t *d6 /1000
print(Torque)
# Lm = 300000 old
Fr = 10133.28
g6dist = round(134/1000,2)
g7dist = round(604/1000,2)
crankDist = round(269/1000,2)
crankDist2 = round(469/1000,2)
shaftL = round(738/1000,2)
distances = [round(99/1000,2),round(134/1000,2), round(169/1000,2), round(269/1000,2), round(469/1000,2), round(569/1000,2), round(604/1000,2), round(639/1000,2),round(738/1000,2)]
alpha = (math.acos(1 - (h/r)))

beta = (1/5*math.sin(alpha))
# print(alpha)git
# print(beta)
P = Lm / math.cos(beta)
T = P* math.sin(alpha+beta)
# print(P, " ", T)
Tcrank = r * T
print(Tcrank)
px = P * math.sin(beta)
py = P * math.cos(beta)
print(px)
# print(px)
#full load
# reactions y 
sumFy = Fr +Fr - py
sumMy = g6dist* Fr + (g7dist*Fr) - (crankDist*1/2*py) - (crankDist2*1/2*py)
by = np.array([sumFy,sumMy])
ay= np.array([[-1,-1] , [0, -shaftL]])
yr = np.linalg.solve(ay,by)
Forcesy = [yr[0],0, Fr, 0,-1/2*py , -1/2*py,0, Fr, 0]
Mgear6y = g6dist * yr[0]
Mcranky = (crankDist -g6dist) * (yr[0] + Fr) + Mgear6y
Mcranky2 = (crankDist2 - crankDist) * (yr[0] + Fr - 1/2*py)+Mcranky
Mgear7y = (g7dist -crankDist2 ) * (yr[0] + Fr -py) + Mcranky2
Mbearingy = (shaftL -g7dist) * (yr[0] + Fr -py +Fr) +Mgear7y
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
sumMx = g6dist* F_t + (g7dist*F_t) - (crankDist*1/2*px) - (crankDist2*1/2*px)
bx = np.array([sumFx,sumMx])
ax = np.array([[-1,-1] , [0, -shaftL]])
xr = np.linalg.solve(ax,bx)
print(xr)
Forcesx = [xr[0],0, F_t,0, -1/2*px , -1/2*px,0, F_t,0] 
Mgear6x = g6dist * xr[0]
Mcrankx = (crankDist -g6dist) * (xr[0] + F_t) + Mgear6x
Mcrankx2 = (crankDist2 - crankDist) *(xr[0]+ F_t - 1/2*px) + Mcrankx
Mgear7x = (g7dist -crankDist2 ) * (xr[0] + F_t -px) + Mcrankx2
Mbearingx = (shaftL -g7dist) * (xr[0] + F_t -px +F_t) +Mgear7x

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
M = math.sqrt((Momentx[3]**2) + (Momenty[3]**2))
print(M)
print(Tcrank)
Msquared16 = (M**2) * (16 ** 2)
Tsquared16 = (Tcrank**2) * (16**2)
d = ((Msquared16 + Tsquared16)/((pi**2) *(maxyield**2))) ** (1/6)
print(d)
Momentymax = Momenty
Momentxmax = Momentx
#no load
# reactions y 
py = 0
sumFy = Fr +Fr - py
sumMy = g6dist* Fr + (g7dist*Fr) - (crankDist*1/2*py) - (crankDist2*1/2*py)
by = np.array([sumFy,sumMy])
ay= np.array([[-1,-1] , [0, -shaftL]])
yr = np.linalg.solve(ay,by)
Forcesy = [yr[0],0, Fr, 0,-1/2*py , -1/2*py,0, Fr, 0]
Mgear6y = g6dist * yr[0]
Mcranky = (crankDist -g6dist) * (yr[0] + Fr) + Mgear6y
Mcranky2 = (crankDist2 - crankDist) * (yr[0] + Fr - 1/2*py)+Mcranky
Mgear7y = (g7dist -crankDist2 ) * (yr[0] + Fr -py) + Mcranky2
Mbearingy = (shaftL -g7dist) * (yr[0] + Fr -py +Fr) +Mgear7y
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
px = 0
sumFx = F_t + F_t - px
sumMx = g6dist* F_t + (g7dist*F_t) - (crankDist*1/2*px) - (crankDist2*1/2*px)
bx = np.array([sumFx,sumMx])
ax = np.array([[-1,-1] , [0, -shaftL]])
xr = np.linalg.solve(ax,bx)
print(xr)
Forcesx = [xr[0],0, F_t,0, -1/2*px , -1/2*px,0, F_t,0] 
Mgear6x = g6dist * xr[0]
Mcrankx = (crankDist -g6dist) * (xr[0] + F_t) + Mgear6x
Mcrankx2 = (crankDist2 - crankDist) *(xr[0]+ F_t - 1/2*px) + Mcrankx
Mgear7x = (g7dist -crankDist2 ) * (xr[0] + F_t -px) + Mcrankx2
Mbearingx = (shaftL -g7dist) * (xr[0] + F_t -px +F_t) +Mgear7x

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
altx = (Momentxmax[3]-Momentx[3])/2
meanx = (Momentxmax[3] + Momentx[3])/2
alty = (Momentymax[3] - Momenty[3])/2
meany = (Momentymax[3] + Momenty[3])/2
ma = (meanx**2 + meany**2 )**(1/2)
Mm = (altx**2 + alty**2 )**(1/2)
Ta = Torque
Tm = 0
UTS = 675e6/1.1
ka = 0.739
kb = 1.51 * ((1000*d)**-.157)
print(kb)
kc = 1
kd = 1
ke = 0.753
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
