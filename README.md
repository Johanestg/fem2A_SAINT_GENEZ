# README.md
Utilisation des programmes de simulation de la méthode des éléments finis: fem2A_SAINT_GENEZ

### Obtention des codes
Pour obtenir les codes, il faut cloner le repository en local puis se placer dans le repository local.


### Tests
Toutes les fonctions implémentées pour la simulation de la méthodes des éléments finis peuvent être testées. Pour se faire, il faut:
1) Ouvrir le fichier main.cpp.
2) Regarder la fonction run_test().
3) Choisir la fonction à tester à partir du nom explicite.
4) Sélectionner la fonction à tester en mettant le booléen correspondant égal à true. Par exemple, pour tester la fonction
   constructeur_elementmapping(), il faut mettre const bool t_constructeur_element_mapping = true.
   #### Afficher tests pour le segment
   Dans le cas des tests:
   - test_constructeur_shapefunction(2,1) (l.57 du main.cpp)
   - test_nbfunction_shapefunction(2,1) (l.58 du main.cpp)
   - test_evaluate_shapefunction(2, 1, 0) (l.59 du main.cpp) 
   - test_evaluategrad_shapefunction(2, 1, 0) (l.60 du main.cpp) 
     Pour afficher le test pour le cas du segment, il faut mettre la valeur 1 au premier paramètre.
    De la même manière, dans le cas des tests:
   - test_assemble_elementary_vector(false) (l.63 du main.cpp) 
   - test_local_to_global_vector(false) (l.64 du main.cpp) 
     Pour afficher le test pour le cas du segment, il faut mettre le booléen true en paramètre.
   #### Attention: lignes à décommenter pour l'affichage du test
   Attention, pour afficher le résultat de test_local_to_global_vector, il faut aller dans le fichier fem.cpp et décommenter les lignes     443, 449 et 457.


### Simulations
Trois types de simulations sont possibles: celle du problème de Dirichlet pur, de Dirichlet à terme source constant et de Dirichlet à terme source sinusoïdal. Pour sélectionner une simulation, il faut:
1) Aller dans le main.cpp
2) Regarder la fonction run_simu()
3) Mettre true au booléen qui caractérise la simulation que l'on souhaite exécuter
#### Changer le maillage 
Pour changer le maillage, il faut changer le string "data/square_fine.mesh" aux lignes 81 ou 82 ou 83 du main.cpp en fonction de la simulation réalisée.

Le fichier solution.bb et le fichier solution.mesh associé au problème sont stockés dans data/output. Le nom des solutions sont modifiables dans le fichier simu.h.

### Compilation et exécution
Pour compiler, entrer la commande make dans le terminal. Pour exécuter lors d'un test, entrer ./build/fem2a -t. Pour exécuter lors d'une simulation, entrer ./build/fem2a -s.

### Affichage des simulations
Pour afficher les simulations, exécuter dans le terminal la commande : path_to_medit data/output/nom_de_la_solution.mesh. Une fois le maillage affiché, appuyer sur la touche m du clavier pour visualiser les valeurs.
