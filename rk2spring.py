import rungekutta
import numpy as np
import matplotlib.pyplot as plt
import sys

## Switched to using system arguments. Old user imputs have been commented out.

##k = 0.0
m = 0.0
x_0 = 0.0
h = 0.0
N = 0.0

for i in range(0, len(sys.argv)):
    ##if (sys.argv[i] == "-k"):
      ##  k = float(sys.argv[i+1])
        ##i += 1
    if (sys.argv[i] == "-m"):
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
        
k = (4*np.pi**2)*m
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
KE = np.zeros(N)        ##kinetic energy
PE = np.zeros(N)
E_total = np.zeros(N)
x_an = np.zeros(N)
v_an = np.zeros(N)

y = np.zeros(2)
y[0] = x_0
##y[0] = float(input("Enter the starting position: "))

err = 10**-8


def update(existingaggregate, newvalue):
    count, mean, M2 = existingaggregatecount += 1
    delta = newvalue - mean
    mean+= delta / count
    delta2 = newvalue - mean
    M2 += delta*delta2
    return count, mean, M2

def finalize(existingAggregate):
    count, mean, M@ = existingAggregate
    if count < 2:
        return float("nan")
    else:
        mean, variance, samplevariance = mean, M2/count, M2/(count-1))
        return mean, variance, samplevariance

for i in range(0,N):
    t[i] = i*h + h
    y = rungekutta.rk2(t[i],y,h,f)
    x[i] = y[0]
    v[i] = y[1]
    KE[i] = 0.5*m*y[1]**2 
    PE[i] = 0.5*k*y[0]**2
    E_total[i] = KE[i] + PE[i]
    x_an[i] = x_0*np.sin(2*np.pi*t[i]+np.pi/2.0)      ##analytical solution for position
    v_an[i] = x_0*2*np.pi*np.sin(2*np.pi*t[i]+np.pi/2.0)    ##analytical solution for velocity
    ##if((x_an[i] - x[i])/x_an[i] <= err):                  ## this caused an infinite loop
        ##print("The difference between the analytical and numerical answers is too large")

    
    
    
    
f = open("rk2spring.dat", "w")
for i in range(0,N):
     f.write(str(t[i]) + " " + str(x[i]) + " " + str(v[i]) + " " + str(KE[i]) + " " + str(PE[i]) + " " + str(E_total[i]) + "\n")
f.close()


rk2graphposition = plt.plot(t, x, 'o')
plt.ylabel('Position')
plt.xlabel('Time')
plt.axis([0, 2, -6, 6])
plt.savefig('rk2graphposition.pdf')

rk2analytical = plt.plot(t, x_an, 'o')
plt.axis([0, 2, -6, 6])
plt.savefig('rk2graphanalytical.pdf')

##rk2graphvelocity = plt.plot(t, v, 'o')
##plt.ylabel('Velocity')
##plt.xlabel('Time')
##plt.axis([0, 10, -40, 40])
##plt.savefig('rk2graphvelocity.pdf')


plt.show()
