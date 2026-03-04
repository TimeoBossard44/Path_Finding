include("map_en_graphe.jl")
include("affichage.jl")
using DataStructures

function algoDijkstra(fname, vD, vA)

    temps = @elapsed begin
        allocation = @allocated begin
            # Matrice pondérée (utilise lireMapPondere de map_en_graphe.jl)
            M = lireMapPondere("../dat/$fname")

            rows = size(M, 1)
            cols = size(M, 2)

            # 4 voisins : haut, droite, bas, gauche
            directions = [(-1,0), (0,1), (1,0), (0,-1)]

            # Distance minimale connue jusqu'à chaque case
            dist = fill(Inf, rows, cols)
            dist[vD[1], vD[2]] = 0

            visites      = falses(rows, cols)
            predecesseur = fill((0,0), rows, cols)

            # File de priorité : clé = (i,j), priorité = distance
            pq = PriorityQueue{Tuple{Int,Int}, Float64}()
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

                    # reste dans la Map et pas un mur
                    if 1 <= ni <= rows && 1 <= nj <= cols && M[ni, nj] != 0
                        newdist = dist[i, j] + M[ni, nj]
                        if newdist < dist[ni, nj]
                            dist[ni, nj] = newdist
                            predecesseur[ni, nj] = pos
                            enqueue!(pq, (ni, nj) => newdist)
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
        end
    end

    return afficherResultat(distance, nb_etats, chemin, temps, allocation)
end
