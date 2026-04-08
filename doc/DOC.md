# Documentation - Path Finding

Projet Julia de recherche de chemins sur des cartes au format(fichiers `.map`), avec plusieurs algorithmes, visualisation PDF et un module expérimental de **cross-docking** (plusieurs AMR).

---

## Structure du dépôt

| Élément | Rôle |
|--------|------|
| `dat/` | Cartes `.map` et données d’exemple |
| `res/` | Sorties générées (cartes avec chemin, PDF) |
| `src/` | Code source Julia |

---

## Point d’entrée

- **`Main.jl`** - Inclut les modules : A*, BFS, Dijkstra, Glouton, `MapToPDF`, `solutionChemin`, `CrossDocking`. Décommenter les appels souhaités pour lancer un algorithme.

---

## Format des fichiers `.map`

En-tête typique :

- `type octile`
- `height H` / `width W`
- ligne `map`
- puis **H** lignes de **W** caractères (`@` mur, `.` sol, `T` terre, `S` sable, `W` eau, etc.)

Les cartes peuvent inclure une section **`amr`** (voir `lireDonneeAMR.jl`) avant `map` pour le cross-docking.

---

## `MapEnGraphe.jl`

### Constante

- **`couts`** - Dictionnaire des coûts par caractère : `'@' => 0` (mur), `'W' => 8`, `'S' => 5` ; les autres cases traversables reçoivent le coût **1** par défaut lors de la lecture.

### Fonction

| Fonction | Description |
|----------|-------------|
| **`lireMapPondere(fichier::String)`** | Lit le fichier après la ligne `map` et renvoie une **matrice** de coûts (`Float64`), une ligne de la carte = une ligne de la matrice. Les murs ont le coût 0 (non traversables dans les algos qui testent `!= 0`). |

---

## Algorithmes de plus court chemin

Toutes ces fonctions prennent :

- **`fname`** : nom ou chemin vers un `.map` (souvent `"nom.map"` ; le code utilise `dat/$fname` en interne là où c’est prévu),
- **`vD`** / **`vA`** : tuples `(ligne, colonne)` départ et arrivée.

Elles retournent en général le résultat via **`afficherResultat`** (distance, nombre d’états, chemin, temps, allocations).

| Fichier | Fonction | Idée |
|---------|----------|------|
| **`BFS.jl`** | **`algoBFS(fname, vD, vA)`** | Parcours en largeur : plus court chemin en **nombre de cases** (coût unitaire par pas). |
| **`Astar.jl`** | **`heuristique(v, vA)`** | Distance de Manhattan \( \|i-i_A\| + \|j-j_A\| \). |
| | **`algoAstar(fname, vD, vA)`** | A* avec file de priorité et heuristique. |
| **`Dijkstra.jl`** | **`algoDijkstra(fname, vD, vA)`** | Plus court chemin avec **coûts positifs** sur la grille (matrice pondérée). |
| **`Glouton.jl`** | **`heuristique(v, vA)`** | Même heuristique Manhattan. |
| | **`algoGlouton(fname, vD, vA)`** | Recherche gloutonne vers l’objectif (pas garantie d’optimalité globale). |
| **`Astar2.0.jl`** | **`heuristique(v, vA)`** | Manhattan. |
| | **`algoAstar2_0(fname, vD, vA)`** | Variante A* utilisée notamment par **`CrossDocking.jl`** ; chemin renvoyé sous forme de **`Queue{Tuple{Int,Int}}`**. |

**Dépendance** : `DataStructures` (files, files de priorité) selon les fichiers.

---

## Affichage et export (`Affichage.jl`, `solutionChemin.jl`, `MapToPDF.jl`)

### `Affichage.jl`

| Fonction | Description |
|----------|-------------|
| **`afficherResultat(distance, nb_etats, chemin, temps, allocation, fname, s::String)`** | Affiche les métriques ; si un chemin existe, écrit une carte solution dans `res/` puis appelle **`map_to_pdf`**. Le nom de sortie dépend de **`fname`** et du libellé **`s`** (ex. préfixe d’algorithme). |

### `solutionChemin.jl`

| Fonction | Description |
|----------|-------------|
| **`ecrire_solution_map(input_map, chemin, output_map)`** | Copie le `.map` et remplace par **`'X'`** toutes les cases du **`chemin`** (coordonnées `(ligne, col)`), uniquement après la ligne `map`. |

### `MapToPDF.jl`

| Élément | Description |
|---------|-------------|
| **`_basename_noext`** | Utilitaire : nom de fichier sans extension. |
| **`parse_map(filename)`** | Lit un `.map` et construit une grille d’entiers (0-6) pour l’affichage : murs, sol, types de terrain, **`X`** pour chemin. |
| **`map_to_pdf(input_file, s::String)`** | Génère un **PDF** (CairoMakie) à partir du fichier carte ; applique une **rotation** (`rotr90`) pour l’orientation ; enregistre sous `res/` avec un nom dérivé de **`input_file`** et du paramètre **`s`**. |

---

## Cross-docking (`CrossDocking.jl`, `lireDonneeAMR.jl`)

### Structure

- **`mutable struct AMR`** (`lireDonneeAMR.jl`) - Champs : position actuelle, file **`chemin`**, arrivée, compteur de temps **`t`**.

### Fonctions

| Fonction | Description |
|----------|-------------|
| **`lireDonneeAMR(fichier)`** | Lit la section après la ligne **`amr`** jusqu’à **`map`** ; remplit un **`Dict{Int, AMR}`** (identifiants d’AMR → état). |
| **`crossDocking(fname)`** | Charge la carte et les AMR depuis **`dat/$fname`**, calcule des chemins avec **`algoAstar2_0`**, simule le temps pas à pas avec gestion de **collisions** et recalcul de chemins. |

