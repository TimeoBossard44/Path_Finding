include("MapEnGraphe.jl")
include("Affichage.jl")
using DataStructures

function algoBFS(fname, vD, vA)

    # Matrice pondérée du fichier fname
    M = lireMapPondere("../dat/$fname")

    temps = @elapsed begin
        allocation = @allocated begin

            rows = size(M, 1)   # nombre de lignes
            cols = size(M, 2)   # nombre de colonnes

            # Les 4 voisins possibles : haut, droite, bas, gauche
            directions = [(-1,0), (0,1), (1,0), (0,-1)]

            dist = fill(Inf, rows, cols)
            dist[vD[1], vD[2]] = 0

            visites      = falses(rows, cols)          # Matrice de booléens
            predecesseur = fill((0,0), rows, cols)   # pour reconstruire le chemin

            file = Queue{typeof(vD)}() # Création d'une file du même type que vD
            enqueue!(file,vD)

            nb_etats = 0 # nombre d'etats visitée
            cheminValide = false

            while !isempty(file)
                pos = dequeue!(file)   # position actuelle (ligne, colonne)
                i, j = pos

                if visites[i, j] # si la case est deja visité" on repart de la boucle while
                    continue
                end

                visites[i, j] = true
                nb_etats += 1

                if pos == vA
                    cheminValide = true
                    break
                end

                for (dl, dc) in directions
                    ni = i + dl   # nouvelle ligne
                    nj = j + dc   # nouvelle colonne

                    # Vérifications :
                    if 1 <= ni <= rows && # on reste dans la Map
                        1 <= nj <= cols && 
                        M[ni, nj] != 0 &&# la case n'est pas un mur
                        !visites[ni, nj] 

                        newdist = dist[i,j] + 1 # (+1) car c'est le coût de 1 pour chaque case
                        dist[ni,nj] = newdist
                        predecesseur[ni, nj] = pos  # Ajoute de la positon actuelle sur la nouvelle position
                        enqueue!(file, (ni, nj))
                    end
                end
            end
        end
    end

    # Reconstruction du chemin
    chemin = []
    if cheminValide
        etape = vA
        while etape !== (0,0)
            push!(chemin, etape)
            etape = predecesseur[etape[1], etape[2]]
        end 
        chemin = reverse(chemin)
    end

    distance = cheminValide ? dist[vA[1], vA[2]] : 0

    return afficherResultat(distance, nb_etats, chemin, temps, allocation)
end 

