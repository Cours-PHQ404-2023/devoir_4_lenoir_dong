import numpy as np
import numba as nb


@nb.jit(nopython=True)
def ising_aleatoire(temperature, taille):
    """ Génére une grille aléatoire de spins.

    Arguments
    ---------
    temperature : Température du système.
    taille : La grille a une dimension taille x taille.
    """
    # On initialise aléatoirement des spins de valeurs -1 ou +1
    # sur un grille de dimension taille x taille.
    spins = np.random.randint(0, 2, (taille, taille))
    spins = 2 * spins - 1
    return Ising(temperature, spins)

# Numba permet de compiler la classe pour qu'elle
# soit plus rapide. Il faut attention car certaines
# opérations ne sont plus permises.
@nb.experimental.jitclass([
    ("temperature", nb.float64),
    ("spins", nb.int64[:, :]),
    ("taille", nb.uint64),
    ("energie", nb.int64),
])
class Ising:
    """ Modèle de Ising paramagnétique en 2 dimensions.

    Représente une grille de spins classiques avec un couplage J = +1 entre
    les premiers voisins.

    Arguments
    ---------
    temperature : Température du système.
    spins : Tableau carré des valeurs de spins
    """

    def __init__(self, temperature, spins):
        self.temperature = temperature
        self.spins = spins
        self.taille = np.shape(spins)[0]
        self.energie = self.calcule_energie()

    def difference_energie(self, x, y):
        """Retourne la différence d'énergie si le spin à la position (x, y)
        était renversé.
        """

        energie = 0 
        initial = 0 
        final = 0 
        n=self.taille 
        for i in range(-1,2,2): 
            # terme avec le voisin de droite. 
            energie = self.spins[x, y] * self.spins[(x + i) % n, y] 
            initial -= energie 
            final += energie 
            #si le spin (x,y) est flippé il suffit de changer le signe. # terme avec le voisin du bas. 
            energie = self.spins[x, y] * self.spins[x, (y + i) % n] 
            initial -= energie 
            final += energie 
        deltaE = final-initial #flippé - non_flippé. 
        return deltaE
   
   
   
       # n = self.taille

        #energie_init = 0        
        #energie_init -= self.spins[x, y] * self.spins[(x + 1) % n, y]
        #energie_init -= self.spins[x, y] * self.spins[x, (y + 1) % n]

        #self.spins[x, y] = -self.spins[x, y]

        #energie_finale = 0
        #energie_finale -= self.spins[x, y] * self.spins[(x + 1) % n, y]
        #energie_finale -= self.spins[x, y] * self.spins[x, (y + 1) % n] #on ne calcule la différence d'énergie que localement à l'endroit du flip
                                                                        #ailleurs la matrice reste la même => même énergie
        #return energie_finale - energie_init

    def iteration_aleatoire(self):
        """Renverse un spin aléatoire avec probabilité ~ e^(-ΔE / T). #si delta_E est négatif : l'énergie diminue et proba de 1

        Cette fonction met à jour la grille avec la nouvelle valeur de spin
        """
        x = np.random.randint(0,self.taille)
        y = np.random.randint(0,self.taille)
        p = np.random.rand()
        T = self.temperature
        delta_e = self.difference_energie(x,y)
        if delta_e < 0 :
            self.spins[x,y] *= -1
        else :
            #prob = rand.random_sample() #retourne un nombre aleatoire en [0,1)
            prob = np.random.rand()
            expo = np.exp(-delta_e/self.temperature)
            if prob <= expo:
                self.spins[x,y] *= -1

        #if p<np.exp(delta_e/T) :
         #   self.spins[x, y] = - self.spins[x, y]

    def simulation(self, nombre_iterations):
        """Simule le système en effectuant des itérations aléatoires.
        """
        for _ in range(nombre_iterations):
            self.iteration_aleatoire()

    def calcule_energie(self):
        """Retourne l'énergie actuelle de la grille de spins."""
        energie = 0
        n = self.taille
        for x in range(n):
            for y in range(n):
                energie -= self.spins[x, y] * self.spins[(x + 1) % n, y]
                energie -= self.spins[x, y] * self.spins[x, (y + 1) % n]
        return energie

    #@property
    def aimantation(self):
        """Retourne l'aimantation actuelle de la grille de spins."""
        aimantation = 0
        n = self.taille
        for x in range(n) :
            for y in range(n) : 
                aimantation += self.spins[x,y]
        return abs(aimantation)

class Observable:
    """Utilise la méthode du binning pour calculer des statistiques
    pour un observable.

    Arguments
    ---------
    nombre_niveaux : Le nombre de niveaux pour l'algorithme. Le nombre
                     de mesures est exponentiel selon le nombre de niveaux.
    """

    def __init__(self, nombre_niveaux):
        self.nombre_niveaux = nombre_niveaux

        # Les statistiques pour chaque niveau
        self.nombre_valeurs = np.zeros(nombre_niveaux + 1, int)
        self.sommes = np.zeros(nombre_niveaux + 1)
        self.sommes_carres = np.zeros(nombre_niveaux + 1)

        # La dernière valeur ajoutée à chaque niveau.
        self.valeurs_precedentes = np.zeros(nombre_niveaux + 1)

        # Le niveau que nous allons utiliser.
        # La différence de 6 donne de bons résultats.
        # Voir les notes de cours pour plus de détails.
        self.niveau_erreur = self.nombre_niveaux - 6

    def ajout_mesure(self, valeur, niveau=0):
        """Ajoute une mesure.

        Arguments
        ---------
        valeur : Valeur de la mesure.
        niveau : Niveau à lequel ajouter la mesure. Par défaut,
                 le niveau doit toujours être 0. Les autres niveaux
                 sont seulement utilisé pour la récursion.
        """
        self.nombre_valeurs[niveau] += 1
        self.sommes[niveau] += valeur
        self.sommes_carres[niveau] += valeur**2
        # Si un nombre pair de valeurs a été ajouté,
        # on peut faire une simplification.
        if self.nombre_valeurs[niveau] % 2 == 0:
            moyenne = (valeur + self.valeurs_precedentes[niveau]) / 2
            self.ajout_mesure(moyenne, niveau + 1)
        else:
            self.valeurs_precedentes[niveau] = valeur

    def est_rempli(self):
        """Retourne vrai si le binnage est complété."""
        dernier_niv = self.sommes[-1] #peut utiliser le bas de la pyramide
        res = True
        if dernier_niv == 0 :
            res = False #le dernier niveau est vide, le binnage n'est pas complété
        return res

    def erreur(self):
        """Retourne l'erreur sur le mesure moyenne de l'observable.

        Le dernier niveau doit être rempli avant d'utiliser cette fonction.
        """
        erreurs = np.zeros(self.nombre_niveaux + 1)
        for niveau in range(self.niveau_erreur + 1):
            nombre_valeurs = self.nombre_valeurs[niveau]
            erreurs[niveau] = np.sqrt(
                (
                    self.sommes_carres[niveau]
                    - self.sommes[niveau]**2 / self.nombre_valeurs[niveau]
                ) / (
                    np.int64(nombre_valeurs)*(np.int64(nombre_valeurs)-1)#on convertit de 32bit vers 64bit, parce qu'on a des nombres de 2**16*(2**16-1) = 32 bit
                )
            )
        return erreurs[self.niveau_erreur] #retourne delta_l du 10e niveau, environ = delta_infty

    def temps_correlation(self):
        """Retourne le temps de corrélation."""
        erreur_finale   = self.erreur()
        nombre_valeurs = self.nombre_valeurs[0]
        erreur_initiale =  np.sqrt((self.sommes_carres[0] - self.sommes[0]**2/self.nombre_valeurs[0])/(np.int64(nombre_valeurs)*(np.int64(nombre_valeurs)-1)))#on convertit de 32bit vers 64bit, parce qu'on a des nombres de 2**16*(2**16-1) = 32 bit
        temps_cor = 0.5*((erreur_finale/erreur_initiale)**2-1)
        return temps_cor/(32*32)

    def moyenne(self):
        """Retourne la moyenne des mesures."""
        somme = self.sommes[self.niveau_erreur]
        nombre_valeur = self.nombre_valeurs[self.niveau_erreur]
        return somme/nombre_valeur