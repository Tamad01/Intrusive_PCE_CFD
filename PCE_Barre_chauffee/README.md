# Solveur EF 1D : Équation de Poisson stochastique par Chaos Polynomial Intrusif

Ce script Python illustre l'application de la méthode de Galerkin stochastique (Chaos Polynomial intrusif) sur un cas d'étude 1D académique : une barre chauffée soumise à une incertitude matériau. Il sert de preuve de concept (Toy Model) pour valider l'algorithme de propagation d'incertitude avant l'intégration dans des solveurs CFD complexes.

**Le Problème Physique**

Le code résout l'équation de Poisson 1D $-k(\xi)\frac{\partial^2 u}{\partial x^2} = f$ sur le domaine $x \in [0, 2]$, avec un terme source $f=2$ et des conditions de Dirichlet homogènes aux bords. La conductivité thermique $k$ est incertaine et modélisée par un développement de Karhunen-Loève tronqué à une variable aléatoire gaussienne : $k(\xi) = 1.0 + 0.2\xi$.

**Approche Mathématique**
* **Discrétisation spatiale :** Méthode des Éléments Finis (EF) 1D classique.
* **Modélisation stochastique :** Développement de la solution sur une base de polynômes d'Hermite tronquée à l'ordre $P=1$.
* **Projection de Galerkin Stochastique :** Assemblage d'un système couplé par blocs (Super-Matrice) résolu simultanément pour extraire le mode moyen ($U_0$) et le mode fluctuant ($U_1$).

**Résultats Clés**

L'exécution du script génère un graphique mettant en évidence deux phénomènes mathématiques majeurs capturés par l'approche intrusive :
1. **La localisation spatiale de l'incertitude :** La variance est nulle aux bords (Dirichlet déterministe) et maximale au centre du domaine.
2. **Le décalage de la moyenne (Mean Shift) :** La solution moyenne stochastique diffère de la solution exacte déterministe évaluée à $k=1.0$, illustrant la capture naturelle des effets non linéaires ($u \propto 1/k$) par le solveur couplé.

**Prérequis et Exécution**
Requiert `numpy` et `matplotlib`.

```bash
python PCE_EF_1D_barre_chauffee.py