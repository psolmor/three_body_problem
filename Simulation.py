import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def rpri(pr):
    return pr

def phipri(pphi,r):
    return pphi/(r**2)

def prpri(pphi,r,phi,w,t):
    return (pphi**2)/(r**3)-G*Mt/(dTL**3)*(1/r**2+ML/(Mt*np.sqrt(1+r**2-2*r*np.cos(phi-w*t))**3)*(r-np.cos(phi-w*t)))

def pphipri(r,phi,w,t):
    return -((G*Mt/(dTL**3)*(ML/Mt)*r)/(np.sqrt(1+r**2-2*r*np.cos(phi-w*t))**3))*np.sin(phi-w*t)

h=10
epsmax=h**5
t=h

#Numero de iteraciones. Introduzco dias
T=14*24*60*60
n=int(T/h)

G=6.67e-11
Mt=5.9736e24
ML=0.07349e24
dTL=3.844e8
w=2.6617e-6
RT=6.378160e6
RL=1.7374

vesc=np.sqrt(2*G*Mt/RT)
v=0.9917*vesc/dTL

#parameters
tetha=np.pi/3
phi=np.pi/4

r=RT/dTL

#Condiciones iniciales r,phi,pr,pphi
y=np.zeros((4,n))

y[0][0]=r
y[1][0]=phi
y[2][0]=v*np.cos(tetha-phi)
y[3][0]=r*v*np.sin(tetha-phi)

posnavex=np.zeros(n)
posnavey=np.zeros(n)
poslunax=np.zeros(n)
poslunay=np.zeros(n)

for j in range(1,n):
    
    k11=h*rpri(y[2][j-1])
    k12=h*phipri(y[3][j-1],y[0][j-1])
    k13=h*prpri(y[3][j-1],y[0][j-1],y[2][j-1],w,t)
    k14=h*pphipri(y[0][j-1],y[2][j-1],w,t)

    k21=h*rpri(y[2][j-1]+k13/2)
    k22=h*phipri(y[3][j-1]+k14/2,y[0][j-1]+k11/2)
    k23=h*prpri(y[3][j-1]+k14/2,y[0][j-1]+k11/2,y[2][j-1]+k12/2,w,t+h/2)
    k24=h*pphipri(y[0][j-1]+k11/2,y[2][j-1]+k12/2,w,t+h/2)

    k31=h*rpri(y[2][j-1]+k23/2)
    k32=h*phipri(y[3][j-1]+k24/2,y[0][j-1]+k21/2)
    k33=h*prpri(y[3][j-1]+k24/2,y[0][j-1]+k21/2,y[2][j-1]+k22/2,w,t+h/2)
    k34=h*pphipri(y[0][j-1]+k21/2,y[2][j-1]+k22/2,w,t+h/2)
        
    k41=h*rpri(y[2][j-1]+k33)
    k42=h*phipri(y[3][j-1]+k34,y[0][j-1]+k31)
    k43=h*prpri(y[3][j-1]+k34,y[0][j-1]+k31,y[2][j-1]+k32,w,t+h)
    k44=h*pphipri(y[0][j-1]+k31,y[2][j-1]+k32,w,t+h)

    y[0][j]=y[0][j-1]+1/6*(k11+2*k21+2*k31+k41)
    y[1][j]=y[1][j-1]+1/6*(k12+2*k22+2*k32+k42)
    y[2][j]=y[2][j-1]+1/6*(k13+2*k23+2*k33+k43)
    y[3][j]=y[3][j-1]+1/6*(k14+2*k24+2*k34+k44)

    posnavex[j-1]=y[0][j-1]*np.cos(y[1][j-1])
    posnavey[j-1]=y[0][j-1]*np.sin(y[1][j-1])
    poslunax[j-1]=1*np.cos(w*t)
    poslunay[j-1]=1*np.sin(w*t)

    t+=h


plt.figure(figsize=(8, 6))
plt.plot(posnavex[:-1], posnavey[:-1], label="nave")
plt.plot(poslunax[:-1],poslunay[:-1], label="Luna" )
plt.xlabel('Posición x')
plt.ylabel('Posición y')
plt.grid(True)
plt.legend()
plt.axis('equal')  # Para que los ejes tengan la misma escala
plt.show()

