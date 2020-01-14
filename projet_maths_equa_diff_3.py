# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 17:06:22 2020

@author: Perrotin
"""
#On choisit la méthode de Heun pour le schéma d'ordre 2
import numpy as np
import matplotlib.pyplot as plt

def solve_heun(f,x0,dt,t0,tf):
    t=[t0]
    x=np.array([x0])
    while t[-1]<tf:
        a=x[-1]+dt/2*(f(t[-1],x[-1])+f(t[-1]+dt,x[-1]+dt*f(t[-1],x[-1])))
        a=np.array([a])
        x=np.concatenate((x,a),axis=0)
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

les_t2,les_x2=solve_heun(f,x0,dt,t0,tf)
les_t,les_x=fonction_reference(dt,t0,tf,x0)

plt.plot(les_t,les_x)
plt.plot(les_t2,les_x2,'r--')
plt.title('Courbes superposées de la fonction solution et de la résolution numérique')
plt.ylabel('Solution numérique en rouge')
plt.xlabel('temps t, dt=0.001')
plt.show()


les_x=np.array(les_x)    
les_x2=np.array(les_x2)

les_x=np.array(les_x)    

plt.plot(les_t,les_x-les_x2)
plt.title("Tracé de l'erreur de troncature en fonction du temps")
plt.xlabel('temps t, dt=0.001')
plt.show()   

#On va maintenant chercher à illustrer le fait que notre schéma soit d'ordre 2.
#Pour cela, on va faire un graphe de l'erreur de troncature maximale en fonction
#de dt, dt allant de 0.0001 à 0.1. Si l'on prend une échelle quadratique en dt, on 
#devrait donc avoir une droite.   

les_dt=np.linspace(0.0001,.1,100)
erreur_tronc=[]
for delta in les_dt:
    les_t2,les_x2=solve_heun(f,x0,delta,t0,tf)
    les_t,les_x=fonction_reference(delta,t0,tf,x0)
    les_x=np.array(les_x)    
    erreur_tronc.append(max(abs((les_x-les_x2))))

les_dt2=[k**2 for k in les_dt]
plt.plot(les_dt2,erreur_tronc)
plt.title('Erreur de troncature maximale en fonction des pas de temps choisi')
plt.xlabel('Pas de temps dt en échelle quadratique')
plt.show() 

#On va tester le solver avec l'équation différentielle y''+10²y=0, et y(0)=0.5 et y'(0)=0
    #sur l'intervalle [0;10], avec un pas de temps de 0.001
    #Il s'agit d'un oscillateur amortie.
    
##La solution est la fonction définie sur R par f(t)=0.5*cos(10*t)
#  
##Paramètre d'entrée
#dt=0.001
#t0=0
#tf=0.6
#x0=[0.5,0]
#
#
#def f(t,x):
#    return np.array([x[1],-10**2*x[0]])
#
#
#def fonction_reference(dt,t0,tf,x0):
#    les_t=[t0]
#    les_x=[x0]
#    while les_t[-1]<tf:
#        les_t.append(les_t[-1]+dt)
#        les_x.append(0.5*np.cos(10*les_t[-1]))
#    return (les_t,les_x)
#
#
#les_t2,les_x2=solve_heun(f,x0,dt,t0,tf)
#les_t,les_x=fonction_reference(dt,t0,tf,x0[0])
#
#
#plt.plot(les_t2,les_x2[:,0],"r--")
#plt.plot(les_t,les_x)
#plt.title('Courbes superposées de la fonction solution et de la résolution numérique')
#plt.ylabel('Solution numérique en rouge')
#plt.xlabel('temps t, dt=0.001')
#plt.show()
#
#les_x=np.array(les_x)    
#
#plt.plot(les_t,les_x-les_x2[:,0])
#plt.title("Tracé de l'erreur de troncature en fonction du temps")
#plt.xlabel('temps t, dt=0.001')
#plt.show()   
#    
#les_dt=np.linspace(0.0001,.015,100)
#erreur_tronc=[]
#for delta in les_dt:
#    les_x2=solve_heun(f,x0,delta,t0,tf)[1]
#    les_x=fonction_reference(delta,t0,tf,x0[0])[1]
#    les_x=np.array(les_x)    
#    les_x2=np.array(les_x2[:,0])
#    erreur_tronc.append(max(abs(les_x-les_x2)))
#
#les_dt2=[k**2 for k in les_dt]
#plt.plot(les_dt2,erreur_tronc)
#plt.title('Erreur de troncature maximale en fonction des pas de temps choisi')
#plt.xlabel('Pas de temps dt en échelle quadratique')
#plt.show()  


 
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


les_t2,les_x2=solve_heun(f,x0,dt,t0,tf)
les_t,les_x=fonction_reference(dt,t0,tf,x0)


plt.plot(les_t2,les_x2[:,0],"r--")
plt.plot(les_t,les_x[:,0])
plt.title('Courbes superposées de la fonction solution \net de la résolution numérique')
plt.ylabel('Solution numérique en rouge \n Position x')
plt.xlabel('temps t, dt=0.001')
plt.show()

les_x=np.array(les_x)    

plt.plot(les_t,abs(les_x[:,0]-les_x2[:,0]))
plt.title("Tracé de l'erreur de troncature en fonction du temps en valeur absolue")
plt.ylabel("Erreur sur la position")
plt.xlabel('temps t, dt=0.001')
plt.show()   

#On choisit comme norme la norme infinie qui renvoie le maximum des composantes d'un vecteur


dt=0.001
t0=0
tf=1
x0=[0.5,0]
    
les_dt=np.linspace(0.0001,0.03,100) #Pas de temps maximal tel qu'on ait au moins 20 itération sur l'intervalle
erreur_tronc=[]
for delta in les_dt:
    les_x2=solve_heun(f,x0,delta,t0,tf)[1]
    les_x=fonction_reference(delta,t0,tf,x0)[1]
    erreur_tronc.append(np.max(abs(les_x-les_x2)))

les_dt2=[k**2 for k in les_dt]

plt.plot(les_dt2,erreur_tronc)
plt.title('Erreur de troncature maximale en fonction des pas de temps choisi')
plt.xlabel('Pas de temps dt en échelle quadratique')
plt.show()  

