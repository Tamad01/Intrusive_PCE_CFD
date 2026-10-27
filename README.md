# Analyse de sensibilité paramétrique par polynômes du chaos intrusifs dans un modèle CFD
En simulation numérique de la mécanique des fluides, les données d'entrée (propriétés du fluide, conditions aux limites, paramètres de forçage) sont rarement connues avec certitude. Ces incertitudes se propagent au cours du calcul et affectent directement la fiabilité des résultats.

Ce projet implémente une méthode de développement en polynômes du chaos (PCE) intrusif reposant sur la projection de Galerkin stochastique au sein du solveur spectral Tchebycube. Contrairement à l'approche classique de Monte-Carlo, dont la convergence lente rend le coût prohibitif pour les solveurs fluides, cette méthode intrusive modifie directement les équations gouvernantes pour capturer fidèlement les effets non linéaires et propager l'incertitude de manière hautement efficace.

## Fonctionnalités et apports
- Analyse de sensibilité analytique : Extraction des moments statistiques (espérance, variance) et des indices de Sobol globaux directement à partir des coefficients spectraux calculés, offrant une analyse de sensibilité pour un coût de calcul quasi nul.
- Évaluation multi-paramètres : Modélisation de la convection de Rayleigh-Bénard incluant simultanément plusieurs paramètres incertains indépendants (couplage thermique et température de paroi chaude) afin de hiérarchiser l'influence de chaque source.