# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 14:37:49 2020

@author: Perrotin
"""
import numpy as np
import matplotlib.pyplot as plt
#Premier programme avec un Euler qui prend en compte des fonction de R dans R,
#On généralisera après.
def solve_euler_explicit(f,x0,dt,t0,tf):
    t=[t0]
    x=[x0]
    while t[-1]<tf:
        x.append(x[-1]+dt*f(t[-1],x[-1]))
        t.append(t[-1]+dt)
    return (t,x)

#On va tester le solver avec l'équation différentielle y'+y=0, et y(0)=1
    #sur l'intervalle [0;10], avec un pas de temps de 0.001
#La solution est la fonction définie sur R par f(x)=exp(-x)
  
#Paramètre d'entrée
dt=0.001
t0=0
tf=10
x0=1


def f(t,x):
    return -x


def fonction_reference(dt,t0,tf,x0):
    les_t=[t0]
    les_x=[x0]
    while les_t[-1]<tf:
        les_t.append(les_t[-1]+dt)
        les_x.append(np.exp(-les_t[-1]))
        
    return (les_t,les_x)

les_t2,les_x2=solve_euler_explicit(f,x0,dt,t0,tf)
les_t,les_x=fonction_reference(dt,t0,tf,x0)

plt.plot(les_t,les_x)
plt.plot(les_t2,les_x2,"r--")
plt.title('Courbes superposées de la fonction solution et de la résolution numérique')
plt.ylabel('Solution numérique en rouge')
plt.xlabel('temps t, dt=0.001')
plt.show()

les_x=np.array(les_x)    
les_x2=np.array(les_x2)

plt.plot(les_t,les_x-les_x2)
plt.title("Tracé de l'erreur de troncature en fonction du temps")
plt.xlabel('temps t, dt=0.001')
plt.show()   
    
#On va maintenant chercher à illustrer le fait que notre schéma soit d'ordre 1.
#Pour cela, on va faire un graphe de l'erreur de troncature maximale en fonction
#de dt, dt allant de 0.0001 à 0.1.

les_dt=np.linspace(0.0001,.1,1000)
erreur_tronc=[]
for delta in les_dt:
    les_x2=solve_euler_explicit(f,x0,delta,t0,tf)[1]
    les_x=fonction_reference(delta,t0,tf,x0)[1]
    les_x=np.array(les_x)    
    les_x2=np.array(les_x2)
    erreur_tronc.append(max(abs(les_x-les_x2)))

plt.plot(les_dt,erreur_tronc)
plt.axis('equal')
plt.title('Erreur de troncature maximale en fonction des pas de temps choisi')
plt.xlabel('Pas de temps dt')
plt.show()     

    
