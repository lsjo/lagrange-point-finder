import math
import random
from matplotlib import pyplot as plt

random.seed(42)

# defining a distance function. it comes in handy
def dist(x1, y1, x2, y2):
    num = ((x1-x2)**2 + (y1-y2)**2) ** 0.5
    return num

# setting the value of G
G = 6.67 * 10**(-11)

# setting the masses of the bodies. For this simulation, it's the earth and sun
mSun = 1.99 * 10**30 
mEarth =  5.97 * 10**29 # Normally to the power of 26, set to 10**29 to exaggerate the effects of the gravity of the earth

# setting the positions of the bodies
pSun = 0, 0
pEarth = 150000000000, 0

# calculating the point both bodies orbit around
bari = (mSun*pSun[0] + mEarth*pEarth[0]) / (mSun + mEarth), (mSun*pSun[1] + mEarth*pEarth[1]) / (mSun+mEarth)
dSun = dist(pSun[0], pSun[1], bari[0], bari[1])
dEarth = dist(pEarth[0], pEarth[1], bari[0], bari[1])
vEarth = (G/dEarth * ( (mSun*dEarth**2)/(dEarth+dSun) **2) )**0.5
period = (2 * math.pi * dEarth) / (vEarth)

# VARIABLES
# G - the gravitational constant, N m2/kg
# mSun, mEarth - masses of objects, kg
# pSun, pEarth - positions of objects, m
# vEarth - velocity of the earth, m/s
# period - the orbital period of the objects, s

def point_grav(pointx, pointy, param="aPointGravVec"):
    dPoint = dist(pointx, pointy, bari[0], bari[1])
    vPoint = (2 * math.pi * dPoint) / period
    aPointCentrip = vPoint ** 2 / dPoint

    x, y = bari[0]-pointx, bari[1]-pointy
    theta = abs(math.atan(y/x)) # note that all angles are in radians.
    aPointCentripVec = abs(math.cos(theta) * aPointCentrip) * (-1 if 0 >= x else 1), abs(math.sin(theta) * aPointCentrip) * (-1 if 0 >= y else 1)

    a, b = pSun[0]-pointx, pSun[1]-pointy
    theta = abs(math.atan(b/a))
    aPointSun = G * mSun / dist(pointx, pointy, pSun[0], pSun[1])**2
    aPointSunVec = abs(math.cos(theta) * aPointSun) * (-1 if 0 >= a else 1), abs(math.sin(theta) * aPointSun) * (-1 if 0 >= b else 1)

    c, d = pEarth[0]-pointx, pEarth[1]-pointy
    theta = abs(math.atan(d/c))
    aPointEarth = G * mEarth / dist(pointx, pointy, pEarth[0], pEarth[1])**2
    aPointEarthVec = abs(math.cos(theta) * aPointEarth) * (-1 if 0 >= c else 1), abs(math.sin(theta) * aPointEarth) * (-1 if 0 >= d else 1)

    aPointGravVec = aPointSunVec[0] + aPointEarthVec[0], aPointSunVec[1] + aPointEarthVec[1]
    aPointGrav = (aPointGravVec[0]**2 + aPointGravVec[1]**2)**0.5

    paramsdict = {
        "aPointGrav": aPointGrav, "aPointEarth": aPointEarth, "aPointSun": aPointSun,
        "aPointCentrip": aPointCentrip, "aPointEarthVec": aPointEarthVec,
        "aPointSunVec": aPointSunVec, "dPoint": dPoint, "vPoint": vPoint,
        "aPointCentripVec": aPointCentripVec, "aPointGravVec": aPointGravVec
        }

    return paramsdict[param]

# creating a set of random values for the algorithm to move to Lagrange points
rpointsx = [random.randint(int(-2E+11), int(2E+11)) for n in range(100)]
rpointsy = [random.randint(int(-2E+11), int(2E+11)) for n in range(100)]

# making a set of random values colinear with the earth and sun
# this makes sure that the points L1, L2 and L3 will be shown on the final graph
# because most points tend to move to L4 and L5 because they are in a more stable equilibrium
rpointsx.extend([random.randint(int(-3E+11), int(3E+11)) for n in range(50)])
rpointsy.extend([0 for n in range(50)])

finalpointsx = []
finalpointsy = []
finalpointsdiff = []

individual_plots = False # set this variable to True for plots of the paths of individual points

iterations = 10
# this determines how many iterations the algorithm goes through. a higher value means
# longer execution times but higher accuracy. a value of 50000 will have high accuracy
# without taking too long to execute.

for n in range(len(rpointsx)):
    point = rpointsx[n], rpointsy[n]

    pltx = [point[0]]
    plty = [point[1]]

    for n in range(iterations):
        centrip = point_grav(point[0], point[1], "aPointCentripVec")
        gravity = point_grav(point[0], point[1], "aPointGravVec")

        # calculating the difference between the ideal centripetal acceleration required to keep the point in
        # orbit and the actual acceleration due to gravity of the point
        diff = centrip[0] - gravity[0], centrip[1]-gravity[1]

        # adding the (amplified) distance vector to the point to move it closer to a Lagrange point.
        # multiplying the difference by 10000000000 was determined to be the largest amount that still moved
        # them predictably
        point = point[0] + diff[0]*10000000000, point[1] + diff[1]*10000000000

        pltx.append(point[0])
        plty.append(point[1])

    ex = [pEarth[0]]
    ey = [pEarth[1]]

    sx = [pSun[0]]
    sy = [pSun[1]]

    if individual_plots:
        # setting the colours. if individual_plots is True, then it will plot a series of points from blue to red,
        # with green representing the final location of the point. in all plots one body will be green and the other will be red.
        
        colours = [n for n in range(len(pltx))]
        plt.scatter(ex, ey, c="g", s=500)
        plt.scatter(sx, sy, c="r", s=1000)
        plt.scatter(pltx, plty, c=colours, cmap="coolwarm")
        plt.scatter([pltx[-1]], [plty[-1]], c="g")
        plt.axis("equal")
        plt.show()

    # adding the final location of the point to lists so they can be plotted
    finalpointsx.append(pltx[-1])
    finalpointsy.append(plty[-1])
    finalpointsdiff.append((diff[0]**2 + diff[1]**2)**0.5)

# plotting the final location of the points once they have converged towards the Lagrange points.
# if the points have not converged, the points furthest from Lagrange points will appear yellow and
# the points closest to the Lagrange points will be dark blue.

plt.title("Location of Lagrange Points")
plt.tick_params(left = False, right = False, labelleft = False, labelbottom = False, bottom = False)
plt.xlim(int(-3E+11), int(3E+11))
plt.ylim(int(-2E+11), int(2E+11))
plt.xlabel("Iterations: "+str(iterations))
plt.scatter(ex, ey, c="g", s=500)
plt.scatter(sx, sy, c="r", s=1000)
plt.scatter(finalpointsx, finalpointsy, c=finalpointsdiff, cmap="viridis")
plt.axis("equal")
plt.show()
