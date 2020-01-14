# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 20:58:58 2020

@author: Perrotin
"""
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
from scipy.integrate import solve_ivp
import matplotlib; 
from matplotlib.pyplot import *
import seaborn as sns
sns.set()

#Paramètre d'entrée
dt=0.05
t0=0
tf=2
x0=[0.5,0]


def solve_ivp_euler_explicit_variable_step(f, t0, x0, t_f, dtmin = 1e-16, dtmax = 0.05, atol = 1e-6):
    dt = dtmax/10; # initial integration step
    ts, xs = [t0], [x0]  # storage variables
    les_dt=[dt]
    t = t0
    ti = 0  # internal time keeping track of time since latest storage point : must remain below dtmax
    x = x0
    while ts[-1] < t_f:
        while ti < dtmax:
            #1
            t_next, ti_next, x_next = t + dt, ti + dt, x + dt * f(x)
            x_back = x_next - dt * f(x_next)
            ratio_abs_error = atol / (linalg.norm(x_back-x)/2)
            dt = 0.9 * dt * sqrt(ratio_abs_error)
            #2
            if dt < dtmin:
                raise ValueError("Time step below minimum")
            elif dt > dtmax/2:
                dt = dtmax/2
            #3
            t, ti, x = t_next, ti_next, x_next
        #4    
        dt2DT = dtmax - ti # time left to dtmax
        les_dt.append(dt2DT)
        t_next, ti_next, x_next = t + dt2DT, 0, x + dt2DT * f(x)
        ts = vstack([ts,t_next])
        xs = vstack([xs,x_next])
        t, ti, x = t_next, ti_next, x_next
        #5
    return (ts, xs.T,les_dt)
'''
Ce programme va calculer numériquement la solution d'une équation différentielle mais avec un
pas de temps variable. Pour cela, il va reprendre les deux formules trouvées plus haut (evaluation de l'erreur
et stratégie d'adaptation du pas de temps). C'est ce que fait le programme entre #1 et #2.
La deuxième boucle while est là car on souhaite certe améliorer la précision de l'évaluation numérique
avec un pas de temps variable, mais on veut quand même garder un espacement temporel entre deux points successifs
de dtmax au niveau de ts. C'est cela qui va permettre de comparer les deux approche de la méthode d'Euler. 
Pour cela, on introduit :
    - dt, la pas de temps variable, on controle qu'il est bien compris entre dtmin et dtmax/2 (pour faire au moins
    deux évaluation)
    - la variable ti, qui représente la distance temporelle (au niveau de l'axe des abscisses)
entre le dernier point calculé et le dernier point stocké dans (ts,xs) qui a pour abscisse un multiple de dtmax.
    - la variable x, qui stocke la valeur de la fonction à l'abscisse ts[-1]+ti.
On ajoute à ti, après chaque itération de la boucle, le pas de temps variable dt calculé entre #1 et #3. On avance
ainsi dans la résolution de la même façon qu'un Euler classique, à ceci près que le pas est variable et les
valeurs calculées de la fonction ne sont plus toutes stocker, mais seulement celle de la dernière itération.

Une fois que ti>tmax, on sort de la boucle et on se recale ensuite sur un point où l'abscisse est un multiple de dtmax
(en reculant, car dt2DT est négatif). On stocke ensuite ce dernier point dans ts et xs, puis on réitère le procédé 
à partir de ce dernier point. 
'''

'''
On va maintenant tester le programme avec l'équation différentielle  y''+0.7y'+10y=0, et y(0)=0.5 et y'(0)=0
La solution est la fonction définie sur R par f(t)=0.5*exp(-0.35*t)*cos(3.14*t)
'''

def f1(x):
    return np.array([x[1],-0.7*x[1]-10*x[0]])


def fonction_reference(dt,t0,tf,x0):
    les_t=[t0]
    les_x=[x0]
    while les_t[-1]<=tf:
        les_t.append(les_t[-1]+dt)
        les_x.append(0.5*np.exp(-les_t[-1]*0.7/2)*np.cos(6.285698/2*les_t[-1]))
    return (les_t,les_x)


les_t2,les_x2,les_dt=solve_ivp_euler_explicit_variable_step(f1,t0,x0,tf)
les_t,les_x=fonction_reference(dt,t0,tf,x0[0])


plot(les_t2,les_x2[0],"r--")
plot(les_t,les_x)
title('Courbes superposées de la fonction solution \net de la résolution numérique\n Euler à pas de temps variable')
ylabel('Solution numérique en rouge \nPosition x')
xlabel('temps t')
show()

'''
On voit que le programme tourne bien, mais il reste à comparer son erreur de troncature
à chaque point évalué avec celles du premier programme Euler, pour réellement juger de son
efficacité
'''

def f(t,x):
    return np.array([x[1],-0.7*x[1]-10*x[0]])

les_t3,les_x3=solve_euler_explicit2(f,x0,dt,t0,tf)
plot(les_t2,les_x2[0]-les_x, 'r--')
plot(les_t3,les_x3[:,0]-les_x)
title("Courbes superposées de l'erreur de troncature \nsuivant les différentes méthodes\ny''+0.7y'+10y=0")
ylabel('Euler à pas de temps variable en rouge \nErreur sur la position x')
xlabel('temps t, dtmax=0.05')
show()

#_______________autres exemple________________

dt=0.01
t0=0
tf=5
x0=1

def f3(t,x):
    return -x

def f4(x):
    return -x


def fonction_reference2(dt,t0,tf,x0):
    les_t=[t0]
    les_x=[x0]
    while les_t[-1]<tf:
        les_t.append(les_t[-1]+dt)
        les_x.append(np.exp(-les_t[-1]))
        
    return (les_t,les_x)

les_t2,les_x2,les_dt=solve_ivp_euler_explicit_variable_step(f4,t0,x0,tf,dtmax=.01)
les_t,les_x=fonction_reference(dt,t0,tf,x0)

les_t3,les_x3=solve_euler_explicit(f3,x0,dt,t0,tf)
les_t,les_x=fonction_reference2(dt,t0,tf,x0)
plot(les_t2,les_x2[0]-les_x, 'r--')
plot(les_t3,np.array(les_x3)-les_x)
title("Courbes superposées de l'erreur de troncature \nsuivant les différentes méthodes\ny'+y=0")
ylabel('Euler à pas de temps variable en rouge \nErreur sur la position x')
xlabel('temps t, dtmax=0.01')
show()

'''
On voit bien qu'à chaque fois, l'erreur de troncature de la méthode à pas variable est bien plus faible
que le programme initial. Le programme est donc plus efficace. 
Dans le premier graphe, cette tendance s'inverse localement mais cela est du
aux variations de la fonctions cosinus.
'''
