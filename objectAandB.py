"""
Object A at temp T
Object B at temp < T
we'll assume that these objects are otherwise identical
we'll assume that they are black bodies
we'll assume these are cylinders and ignore the top ends
    radius r
    height h
    separation d
    d >> r
We assume that the temperature of the surroundings is unchanged.
We'll assume that no additional energy is being added.
    That is to say the two objects are heated to their starting temps
    independantly, then the power is switched off and they are
    immediately brought together.
We'll assume that the temperature of each object is uniform.
    So they radiate in a uniform way.
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

#some constants
s = 5.670373e-8 # [W m^-2 K^-4]
Te = 300.0 # [K] Environmental temperature

# set up the time things
t = 0.2 # [s] time to simulate
dt = 0.0001 # [s] time step
t = np.arange(0,t,dt)

# some of the object properties
r = 1.0 # [m]
h = 1.0 # [m]
d = 10.0 # [m]
A = 2.0 * np.pi * r * h # Area of cylinder ignoring tops
m = 1.0 # [kg] mass
c = 100.0 # [J/kg-K] specific heat capacity


def simulation(Taz, Tbz, t, ca, cb):

    Ta = np.array([Taz])
    Tb = np.array([Tbz])

    for i in range(len(t)-1):
        # Power radiated by object A, B from stefan boltzman law
        Pa = A * s * ((Ta[-1])**4)
        Pb = A * s * ((Tb[-1])**4)

        # At B this power is equally distributed in a cylinder of radius d
        # the power/m^2 here is thus
        Pad = Pa / (2*np.pi*d*h)


        # The area of object B that 'sees' this is 2 * r * h, a rectangle.
        # thus power added to b is
        PowerOnB = Pad * 2*r * h


        # Cancelling some terms we can do the same for the power on A
        # due to the radiation from B
        PowerOnA = (A * s * (Tb[-1])**4 * r)/(np.pi * d)


        # Thus the total power loss for A is
        # whatever it radiates
        # minus whatever it gets off of B
        # minus whatever it gets from the surroundings
        # and the energy is just this multiplied by the time step.
        Ela = (Pa - PowerOnA - (A*s*(Te**4)))*dt
        Elb = (Pb - PowerOnB - (A*s*(Te**4)))*dt


        # So the temperatures change
        Ta = np.append(Ta, Ta[-1] - (Ela/(m*ca)))
        Tb = np.append(Tb, Tb[-1] - (Elb/(m*cb)))


    return ((t,Ta,Tb))



def noB(Taz, t):

    Ta = np.array([Taz])

    for i in range(len(t)-1):
        # Power radiated by object Afrom stefan boltzman law
        Pa = A * s * ((Ta[-1])**4)


        # This time these is no energy added from B
        # the energy loss on A is simply what it radiates
        # minus what is supplied by the environment
        Ela = (Pa - (A*s*(Te**4)))*dt

        # So the temperatures change
        Ta = np.append(Ta, Ta[-1] - (Ela/(m*c)))

    return ((t,Ta))


# So let's get some data for if T(B) is pretty cold compared to T(A)
a = simulation(10000.0, 2000.0, t, c, c)

# And let's get some data for if T(B) is almost the same temperature as T(A)
b = simulation(10000.0, 9999.0, t, c, c)

# And what if there was no B and A just radiated to the environment
c = noB(10000.0, t)


# And now plot these things
fig = plt.figure(figsize=(12,7))
plt.plot(a[0], a[1], '-o', label="T(A) - B is cold")
plt.plot(a[0], a[2], '->', label="T(B) - B is cold")

plt.plot(b[0], b[1], '-v', label="T(A) - B is hot")
plt.plot(b[0], b[2], '-*', label="T(B) - B is hot")

plt.plot(c[0], c[1], '-s', label="T(A) - No Object B")
plt.xlabel("Time [s]")
plt.ylabel("Temperature [K]")
plt.legend()
#plt.show()

fig.savefig('fullGraph.png')
plt.xlim(0.1, 0.102)
plt.ylim(950,1000)
fig.savefig('zoom.png')
