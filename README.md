# 🚀 Path Finding (PF)

> Le but de ce projet est de déterminer le plus court chemin entre deux points connu.
> Pour cela nous allons utilise plusieur algorithmes de Path Finding.

---

## 📋 Table des matières

- [Aperçu](#-aperçu)    
- [Prérequis](#-prérequis)
- [Installation](#-installation)
- [Utilisation](#-utilisation)
- [Structure du projet](#-structure-du-projet)
---

## 🎯 Aperçu
 Des applications du Path Finding (PF)

   ### • Les jeux vidéos
        ◦ comment éviter des obstacles habilement, comment trouver le chemin le plus efficace
        sur différents terrains...

   ### • Les véhicules autonomes
        ◦ comment se déplacer dans un espace contenant des obstacles connus et inconnus en
        minimisant les collisions et le temps moyen pour accomplir le chemin...

   ### • La logistique d’entrepôt (MAPF : variante du PF dite Multi-agent path finding)
        ◦ comment déplacer un ensemble d’agents depuis leur point de départ individuel a leur
        point d’arrivée sans collision

   ### • La circulation de trains (MTPF : variante du MAPF dite Multi-Train Path Finding)
        ◦ comment établir un plan de circulation de trains sans collision.



---

## 🛠 Prérequis

Avant de commencer, assure-toi d'avoir installé :

- [Julia](https://julialang.org/downloads/)
- [Git](https://git-scm.com/)

---

## 📦 Installation

### 1. Cloner le dépôt
```bash
git clone https://github.com/TimeoBossard44/Path_Finding.git
cd Path_Finding
```

### 2. Installer les dépendances Julia
```julia
# Lance Julia dans le terminal
julia

# Dans le REPL Julia
Appuier sur ']'

# Puis ajouter le package :
add DataStructures

```

---

## 🚀 Utilisation

### Lancement basique dans le REPL
```julia
include("Main.jl")
```

### Lancer les tests
```julia
algoBFS(fname, D, A)
algoDijkstra(fname, D, A)
algoGlouton(fname, D, A)
algoAstar(fname, D, A)

avec comme paramètres :

• fname | type : String | exemple : "Aurora.map"
• D | type : Tuple{Int64, Int64} | exemple : (170, 10)
• A | type : Tuple{Int64, Int64} | exemple : (200, 120)
```

---

## 📁 Structure du projet
```
Projet_d_informatique_scientifique/
├── src/
│   ├── Main.jl          # Point d'entrée principal
│   ├── BFS.jl          
│   ├── Dikjstra.jl
│   ├── Astar.jl
│   ├── Glouton.jl
│   ├── Affichage.jl
│   └── MapEnGraphe.jl  
├── dat/
│   ├── bootybay.map       
│   ├── Aurora.map
│   └── Berlin_0_256.map
├── res/
│   └── resultat.txt  
├── docs/                # Documentation
└── README.md
```
