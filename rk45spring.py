import numpy as np
import rungekutta   
import matplotlib.pyplot as plt
import sys


m = 0.0
x_0 = 0.0
h = 0.0
N = 0.0
err_tol = 0.0

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
    elif (sys.argv[i] == "-err_tol"):
        err_tol = float(sys.argv[i+1])
        i += 1
        
k = (4*np.pi**2)*m        

def f(t, y):
    if (len(y) != 2):
        print("Improper y array size in function f(t,y). Can only be 2 elements.")
    F = np.zeros(2)
    F[0] = y[1]
    F[1] = -(k/m)*y[0]
    return F

t_max = N*h
t_now = 0.0
t = np.empty(0)
x = np.empty(0)
v = np.empty(0)
ht = np.empty(0)
h_min = h/64.0
h_max = h*64.0
KE = np.empty(0)
PE = np.empty(0)
E_total = np.empty(0)
x_an = np.empty(0)
v_an = np.empty(0)

err=10**-8

y = np.zeros(2)
y[0] = x_0
i = 0

while (t_now < t_max):
    ht = np.append(ht, h)
    t_now += h
    t = np.append(t, t_now)
    h, y = rungekutta.rk45(t[i-1], y, h, f, err_tol, h_min, h_max)
    x = np.append(x, y[0])
    v = np.append(v, y[1])
    KE = np.append(KE, 0.5*m*y[1]**2)
    PE = np.append(PE, 0.5*k*y[0]**2)
    E_total = np.append(E_total, PE[i] + KE[i])
    x_an = np.append(x_an, x_0*np.cos(2*np.pi*t[i]))
    v_an = np.append(v_an, -x_0*2*np.pi*np.sin(2*np.pi*t[i]))
    i += 1

average = sum(E_total)/len(E_total)
print(average)
standarddev = np.std(E_total)
print(standarddev)


rk45graphposition = plt.plot(t, x, 'o')
plt.ylabel('Position')
plt.xlabel('Time')
plt.axis([0, 10, -6, 6])
plt.savefig('rk45graphposition.pdf')
##rk2analytical = plt.plot(t, x_an, 'o')
##plt.axis([0, 10, -6, 6])
#plt.savefig('rk2positioncompare.pdf')
plt.close()


rk4graphvelocity = plt.plot(t, v, 'o')
plt.ylabel('Velocity')
plt.xlabel('Time')
plt.axis([0, 10, -40, 40])
plt.savefig('rk45graphvelocity.pdf')
##rk2analyticalV = plt.plot(t, v_an,'o')
##plt.axis([0, 10, -40, 40])
##plt.savefig('rk2velocitycompare.pdf')    
    
    
f = open("rk45spring.dat", "w")
for i in range(0, len(t)):
    f.write(str(t[i]) + " " + str(x[i]) + " " + str(v[i]) + " " + str(KE[i]) + " " + str(PE[i]) + " " + str(E_total[i]) + " " + str(ht[i]) + "\n")
f.close()





plt.show()




