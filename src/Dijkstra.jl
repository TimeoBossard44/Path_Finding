include("MapEnGraphe.jl")
include("Affichage.jl")
using DataStructures


function algoDijkstra(fname, vD, vA)

    # Matrice pondérée du fichier fname
    M = lireMapPondere("dat/$fname") 

    temps = @elapsed begin
        allocation = @allocated begin

            rows = size(M, 1)
            cols = size(M, 2)

            # 4 voisins : haut, droite, bas, gauche
            directions = [(-1,0), (0,1), (1,0), (0,-1)]

            # Distance minimale connue jusqu'à chaque case
            dist = fill(Inf, rows, cols) # Création d'une matrice avec que des infinis
            dist[vD[1], vD[2]] = 0

            visites      = falses(rows, cols)
            predecesseur = fill((0,0), rows, cols)

            # File de priorité : clé = (i,j), priorité = distance
            pq = PriorityQueue{Tuple{Int,Int}, Int64}() # File de priorité pour avoir accées à la plus court distance en premier
            enqueue!(pq, vD => 0) 

            nb_etats = 0
            cheminValide = false

            while !isempty(pq)
                pos = dequeue!(pq)   # récupère la case de plus petite distance
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
                    ni = i + dl
                    nj = j + dc

                    # Reste dans la Map et pas un mur
                    if 1 <= ni <= rows && 1 <= nj <= cols && M[ni, nj] != 0
                        newdist = dist[i, j] + M[ni, nj]
                        if newdist < dist[ni, nj] # Nouvelle distance plus court que l'ancienne on la remplace
                            dist[ni, nj] = newdist
                            predecesseur[ni, nj] = pos # Ajoute de la positon actuelle sur la nouvelle position
                            if haskey(pq, (ni, nj))
                                pq[(ni, nj)] = newdist
                            else
                                enqueue!(pq, (ni, nj) => newdist) 
                            end
                        end
                    end
                end
            end
        end
    end

    # Reconstruction du chemin
    chemin = Tuple{Int,Int}[]
    if cheminValide
        etape = vA
        while etape !== (0,0)
            push!(chemin, etape)
            etape = predecesseur[etape[1], etape[2]]
        end
        chemin = reverse(chemin)
    end

    distance = cheminValide ? dist[vA[1], vA[2]] : 0

    return afficherResultat(distance, nb_etats, chemin, temps, allocation, fname,"solutionDijkstra")
end
