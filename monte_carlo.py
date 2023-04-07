from ising import Observable
from ising import Ising
from ising import ising_aleatoire
import numpy as np
import numba as nb
import matplotlib.pyplot as plt

temperatures = np.arange(1.0,4.1,0.1)
N=32

energie_moy     = open('energie_moyenne.txt','w')
aimantation_moy = open('aimantation_moyenne.txt','w')

energie_err     = open('energie_erreur.txt','w')
aimantation_err = open('aimantation_erreur.txt','w')

energie_tau     = open('energie_tau.txt','w')
aimantation_tau = open('aimantation_tau.txt','w')

n              = 1000000
nombre_niveaux = 16
pas            = 1000

for T in temperatures :
    ising_init = ising_aleatoire(T,N) #on initialise une grille à une température donnée
    print(T)

    energie     = Observable(nombre_niveaux) # doit être constant = 16
    aimantation = Observable(nombre_niveaux) #on initialise nos deux observables avec 16 niveaux
                       
    ising_init.simulation(n) #on fait évoluer la grille aléatoirement n fois pour réchauffer le système
    while not(energie.est_rempli()) :
        ising_init.simulation(1000)
        energie.ajout_mesure(ising_init.calcule_energie()/N**2) 
        aimantation.ajout_mesure(ising_init.aimantation()/N**2) #on ajoute une mesure d'observable à tous les "pas" itérations

    if energie.est_rempli() :
        energie_moy.write(str(energie.moyenne())+',')
        energie_err.write(str(energie.erreur())+',')
        energie_tau.write(str(energie.temps_correlation())+',') #on écrit dans le fichier si le binning est complété
    if aimantation.est_rempli() :
        aimantation_moy.write(str(aimantation.moyenne())+',')
        aimantation_err.write(str(aimantation.erreur())+',')
        aimantation_tau.write(str(aimantation.temps_correlation())+',')

energie_moy     = open('energie_moyenne.txt','r')
aimantation_moy = open('aimantation_moyenne.txt','r')

energie_err     = open('energie_erreur.txt','r')
aimantation_err = open('aimantation_erreur.txt','r')

energie_tau     = open('energie_tau.txt','r')
aimantation_tau = open('aimantation_tau.txt','r')

E_moy = energie_moy.readlines()[0][0:-1].split(',')
E_moy = [float(el) for el in E_moy]

aim_moy = aimantation_moy.readlines()[0][0:-1].split(',')
aim_moy = [float(el) for el in aim_moy]

E_err = energie_err.readlines()[0][0:-1].split(',')
E_err = [float(el) for el in E_err]

aim_err = aimantation_err.readlines()[0][0:-1].split(',')
aim_err = [float(el) for el in aim_err]

E_tau = energie_tau.readlines()[0][0:-1].split(',')
E_tau = [float(el) for el in E_tau]

aim_tau = aimantation_tau.readlines()[0][0:-1].split(',')
aim_tau = [float(el) for el in aim_tau]

plt.figure()
plt.errorbar(temperatures,E_moy,yerr=E_err,fmt='o--',label='Energie')
plt.errorbar(temperatures,aim_moy,yerr=E_err,fmt='s--',label='Aimantation')
plt.xlabel('Temperature')
plt.ylabel('Valeur moyenne')
plt.legend()
plt.savefig('moy.png')
plt.show()

plt.figure()
plt.plot(temperatures,E_tau,'o--',label='Energie')
plt.plot(temperatures,aim_tau,'s--',label='Aimantation')
plt.xlabel('Temperature')
plt.ylabel('Temps autocorrélation')
plt.legend()
plt.savefig('tau.png')
plt.show()