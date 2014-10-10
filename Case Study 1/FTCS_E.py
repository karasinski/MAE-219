# 1-D Unsteady Heat Transfer
# J-P Delplanque -- UC Davis (2014)
#
"""
Solution of 1-D non-dimensional unsteady Heat Transfer problem using FTCS explicit method
Created: OCT-2014

@author: J-P Delplanque 
"""

from numpy import *
import matplotlib.pyplot as plt

# Problem Parameters
L         = 1.           # Domain lenghth       [n.d.]
T0        = 0.           # Initial temperature  [n.d.]
T1        = 1.           # Boundary temperature [n.d.]
t_start   = 0.
t_end     = 0.1
s         = 1. / 6.
N         = 21


# Set-up Mesh
#dx        = L / (N - 1)
x         = linspace(0,L,N)
dx        = x[1]-x[0]
 
# Calculate time-step
dt        = s*dx**2.0   
time      = 0.

# Initial Condition
Tnew      = [T0]*N

#Boundary conditions
Tnew[0]   = T1
Tnew[N-1] = T1

Told      = Tnew

plt.axis([0,L,T0,T1])
plt.xlabel('Length [nd]')
plt.ylabel('Temperature [nd]')

while time <= t_end:

    for i in range(1,N-1):
        Tnew[i]= s*Told[i + 1] + (1-2.0*s)*Told[i] + s*Told[i - 1]

    plt.plot(x,Tnew,linewidth=1)
    time = time + dt
    Told = Tnew


plt.show()
print('\n Done.\n')
