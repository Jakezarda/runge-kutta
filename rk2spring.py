import rungekutta
import numpy as np
import matplotlib.pyplot as plt
import sys

## Switched to using system arguments. Old user imputs have been commented out.

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

##k = float(input("Enter value of the spring constant: "))
##m = float(input("Enter value of the mass: "))

def f(t,y):
    if (len(y) !=2):
        print("Improper y array size in function f(t,y). Must be 2 elements only.")
    F = np.zeros(2)
    F[0] = y[1]
    F[1] = -(k/m)*y[0]
    return F

##h = float(input("Enter the desired time step: "))
##N = int(input("Enter number of time steps: "))
t = np.zeros(N)
x = np.zeros(N)
v = np.zeros(N)

y = np.zeros(2)
y[0] = x_0
##y[0] = float(input("Enter the starting position: "))

for i in range(0,N):
    t[i] = i*h + h
    y = rungekutta.rk2(t[i],y,h,f)
    x[i] = y[0]
    v[i] = y[1]
    
    
f = open("rk2spring.dat", "w")
for i in range(0,N):
    f.write(str(t[i]) + " " + str(x[i]) + " " + str(v[i]) + "\n")
f.close()


plt.plot(t, x, 'o')

plt.show()
