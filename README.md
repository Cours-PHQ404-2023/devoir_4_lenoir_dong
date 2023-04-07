# devoir_4_lenoir_dong
Les fichiers de ce quatrième devoir sont organisés comme suit. Le rapport principal se nomme PHQ404_Devoir_4.pdf 

Le fichier principal se nomme monte_carlo.py et contient tout le protocole à défiler pour l'analyse de la grille de spins. L'environnement est détaillé dans le fichier pyproject.toml.

Ce fichier principal va utiliser les classes définies dans le second fichier ising.py. Ce fichier contient la classe Ising, décrivant la grille de spins, ainsi que les méthodes utiles pour le calcul des observables. La classe Observable quant à elle, sera utilisée pour la méthode du binning et une analyse plus détaillée de ces observables.

Les fichiers .txt contiennent les données collectées lors de l'algorithme du binning. Il y a deux types de fichier, ceux commençant par "E" pour l'énergie et ceux commençant par "aim" pour l'aimantation. Il y a ensuite trois sous-types de fichiers. Ceux avec "moy" décrivent les moyennes des observables.» Ceux avec "err" décrivent les erreurs. Ceux avec "tau" décrivent les temps d'autocorrélation.

Les deux images "moy.png" et "tau.png" contiennent respectivement les graphes des évolutions des moyennes et des temps d'autocorrélation des deux observables en fonction de la température. 
