import numpy as np
import rungekutta   
import matplotlib.pyplot as plt
import sys


m = 0.0
x_0 = 0.0
h = 0.0
N = 0

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
E_total = np.zeros(N)
x_an = np.zeros(N)
v_an = np.zeros(N)
err = 10**-8

y = np.zeros(2)
y[0] = x_0


for i in range(0, N):
    t[i] = i*h + h
    y = rungekutta.rk4(t[i], y, h, f)
    x[i] = y[0]
    v[i] = y[1]
    KE[i] = 0.5*m*y[1]**2       ##used to track kinetic energy for the assignment 
    PE[i] = 0.5*k*y[0]**2       ##used to track potential energy for the assignment 
    E_total[i] = KE[i] + PE[i]
    x_an[i] = x_0*np.cos(2*np.pi*t[i])                ## analytical solution for position
    v_an[i] = -x_0*2*np.pi*np.sin(2*np.pi*t[i])        ## analytical solution for velocity
    
    
    
##    if (np.abs(x_an[i] - x[i])/x_an[i] >= err):
##        print("The difference between the analytical and numerical solutions for position is too large")
##        break
##    elif (np.abs(x_an[i] - x[i])/x_an[i] < err):
##        print("This time step produces an accurate result for the position.")
    
##    if (np.abs((v_an[i] - v[i])/v_an[i]) >= err):
##        print("The difference between the analytical and numerical solutions for velocity is too large")
##        break
##    elif (np.abs((v_an[i] - v[i])/v_an[i]) < err):
##        print("This time step produces an accurate result for the velocity.")
    
    
#if (np.abs(np.average(x_an) - np.average(x))/np.average(x_an) >= err):
#    print("The difference between the analytical and numerical solutions for position is too large")

#elif (np.abs(np.average(x_an) - np.average(x))/np.average(x_an) < err):
 #   print("This time step produces an accurate result for the position.")
    
#if (np.abs(np.average(v_an) - np.average(v))/np.average(v_an) >= err):
#    print("The difference between the analytical and numerical solutions for velocity is too large")
    
#elif (np.abs(np.average(v_an) - np.average(v))/np.average(v_an) < err):
#    print("This time step produces an accurate result for the velocity.")
    


average = sum(E_total)/len(E_total)
print(average)
standarddev = np.std(E_total)
print(standarddev)



    
rk4graphposition = plt.plot(t, x, 'o')
plt.ylabel('Position')
plt.xlabel('Time')
plt.axis([0, 10, -6, 6])
plt.savefig('rk4graphposition.pdf')
##rk2analytical = plt.plot(t, x_an, 'o')
##plt.axis([0, 10, -6, 6])
#plt.savefig('rk2positioncompare.pdf')
plt.close()


rk4graphvelocity = plt.plot(t, v, 'o')
plt.ylabel('Velocity')
plt.xlabel('Time')
plt.axis([0, 10, -40, 40])
plt.savefig('rk4graphvelocity.pdf')
##rk2analyticalV = plt.plot(t, v_an,'o')
##plt.axis([0, 10, -40, 40])
##plt.savefig('rk2velocitycompare.pdf')    
    
f = open("rk4spring.dat", "w")
for i in range(0, N):
    f.write(str(t[i]) + " " + str(x[i]) + " " + str(v[i]) + " " + str(KE[i]) + " " + str(PE[i]) + " " + str(E_total[i]) + "\n")
f.close()





