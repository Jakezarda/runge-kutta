import numpy as np
import rungekutta   
import matplotlib.pyplot as plt
import sys

k = 0.0
m = 0.0
x_0 = 0.0
h = 0.0
N = 0.0

for i in range(0, len(sys.argv)):
    if (sys.argv[i] == "-k"):
        i += 1
    elif (sys.argv[i] == "-m"):
        m = float(sys.argv[i+1])
        i += 1
    elif (sys.argv[i] == "-x_0"):
        x_0 = float(sys.argv[i+1])
        i += 1
    elif (sys.argv[i] == "-h"):
        h = float(sys.argv[i+1])
        i += 1
    elif (sys.argv[i] == "-N"):
        N = int(sys.argv[i+1])
        i += 1
        
def f(t, y):
    if (len(y) != 2):
        print("Improper y array size in function f(t,y). Can only be 2 elements.")
    F = np.zeros(2)
    F[0] = y[1]
    F[1] = -(k/m)*y[0]
    return F

t = np.zeros(N)
x = np.zeros(N)
v = np.zeros(N)
KE = np.zeros(N)        ##kinetic energy
PE = np.zeros(N)        ##potential energy

y = np.zeros(2)
y[0] = x_0


for i in range(0, N):
    t[i] = i*h + h
    y = rungekutta.rk4(t[i], y, h, f)
    x[i] = y[0]
    v[i] = y[1]
    KE[i] = 0.5*m*y[1]**2       ##used to track kinetic energy for the assignment (may be wrong)
    PE[i] = 0.5*k*y[0]**2       ##used to track potential energy for the assignment (may be wrong)
    
f = open("rk4spring.dat", "w")
for i in range(0, N):
    f.write(str(t[i]) + " " + str(x[i]) + " " + str(v[i]) + str(KE[i]) + " " + str(PE[i]) + "\n")
f.close()

plt.plot(t, x, 'o')

plt.show()




