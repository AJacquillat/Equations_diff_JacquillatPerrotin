# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 15:31:02 2020

@author: Perrotin
"""
'''
Là c'est Euler avec des fonctions de R dans R^n
'''
import numpy as np
import matplotlib.pyplot as plt

#Dans un second temps, on écrit une fonction capable de gérer la méthode d'Euler
#avec des fonctions de R dans R^n

def solve_euler_explicit(f,x0,dt,t0,tf):
    t=[t0]
    x=np.array([x0])
    while t[-1]<tf:
        a=x[-1]+dt*f(t[-1],x[-1])
        a=np.array([a])
        x=np.insert(x,len(x),a,axis=0)
        t.append(t[-1]+dt)
    return (t,x)

#On va tester le solver avec l'équation différentielle y''+0.7y'+10y=0, et y(0)=0.5 et y'(0)=0
    #sur l'intervalle [0;10], avec un pas de temps de 0.001
    #Il s'agit d'un oscillateur amortie.
    
#La solution est la fonction définie sur R par f(t)=0.5*exp(-0.35*t)*cos(3.14*t)
  
#Paramètre d'entrée
dt=0.001
t0=0
tf=10
x0=[0.5,0]


def f(t,x):
    return np.array([x[1],-0.7*x[1]-10*x[0]])


def fonction_reference(dt,t0,tf,x0):
    les_t=[t0]
    les_x=[x0]
    while les_t[-1]<tf:
        les_t.append(les_t[-1]+dt)
        les_x.append(0.5*np.exp(-les_t[-1]*0.7/2)*np.cos(6.285698/2*les_t[-1]))
    return (les_t,les_x)


les_t2,les_x2=solve_euler_explicit(f,x0,dt,t0,tf)
les_t,les_x=fonction_reference(dt,t0,tf,x0[0])


plt.plot(les_t2,les_x2[:,0],"r--")
plt.plot(les_t,les_x)
plt.title('Courbes superposées de la fonction solution \net de la résolution numérique')
plt.ylabel('Solution numérique en rouge \nPosition x')
plt.xlabel('temps t, dt=0.001')
plt.show()

les_x=np.array(les_x)    

plt.plot(les_t,les_x-les_x2[:,0])
plt.title("Tracé de l'erreur de troncature en fonction du temps")
plt.xlabel('temps t, dt=0.001')
plt.ylabel('Position x')
plt.show()   

#Une dernière équation différentielle pour la route !
#On va tester le solver avec l'équation différentielle y''+10²y=0, et y(0)=0.5 et y'(0)=0
    #sur l'intervalle [0;10], avec un pas de temps de 0.001
    #Il s'agit d'un oscillateur amortie.
    
#La solution est la fonction définie sur R par f(t)=0.5*cos(10*t)
  
#Paramètre d'entrée
dt=0.001
t0=0
tf=0.6
x0=[0.5,0]


def f(t,x):
    return np.array([x[1],-10**2*x[0]])


def fonction_reference(dt,t0,tf,x0):
    les_t=[t0]
    les_x=np.array([x0])
    while les_t[-1]<tf:
        les_t.append(les_t[-1]+dt)
        les_x=np.insert(les_x, len(les_x), np.array([0.5*np.cos(10*les_t[-1]),-0.5*10*np.sin(10*les_t[-1])]), axis=0)
    return (les_t,les_x)


les_t2,les_x2=solve_euler_explicit(f,x0,dt,t0,tf)
les_t,les_x=fonction_reference(dt,t0,tf,x0)


plt.plot(les_t2,les_x2[:,0],"r--")
plt.plot(les_t,les_x[:,0])
plt.title('Courbes superposées de la fonction solution \net de la résolution numérique')
plt.ylabel('Solution numérique en rouge \n Position x')
plt.xlabel('temps t, dt=0.001')
plt.show()

les_x=np.array(les_x)    

plt.plot(les_t,abs(les_x[:,0]-les_x2[:,0]))
plt.title("Tracé de l'erreur de troncature en fonction du temps \nen valeur absolue")
plt.ylabel("Erreur sur la position")
plt.xlabel('temps t, dt=0.001')
plt.show()   

#On choisit comme norme la norme infinie qui renvoie le maximum des composantes d'un vecteur
    
les_dt=np.linspace(0.0001,0.03,100) #Pas de temps maximal tel qu'on ait au moins 20 itération sur l'intervalle
erreur_tronc=[]
for delta in les_dt:
    les_x2=solve_euler_explicit(f,x0,delta,t0,tf)[1]
    les_x=fonction_reference(delta,t0,tf,x0)[1]
    erreur_tronc.append(np.max(abs(les_x-les_x2)))

plt.plot(les_dt,erreur_tronc)
plt.title('Erreur de troncature maximale en fonction des pas de temps choisi')
plt.xlabel('Pas de temps dt')
plt.show()  

#Ca a pas une allure ouf mais tu peux dire qu'on dégage clairement une droite 
#"moyenne" si l'on veut, donc ça marche ;-)
print('Il se dégage de ce graphe une sorte de droite, \
      donc on peut bien majorer l erreur de troncature par une droite linéaire de dt.')