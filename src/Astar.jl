include("MapEnGraphe.jl")
include("Affichage.jl")
using DataStructures

function heuristique(v, vA, c_min)
    (i, j), (iA, jA) = v, vA
    (abs(i - iA) + abs(j - jA)) * c_min
end

function algoAstar(fname, vD, vA)

    temps = @elapsed begin
        allocation = @allocated begin
            M = lireMapPondere("../dat/$fname")

            rows = size(M, 1)
            cols = size(M, 2)

            # Coût minimum sur les cases franchissables (pour une heuristique admissible)
            couts_positifs = M[M .> 0]
            c_min = isempty(couts_positifs) ? 1.0 : Float64(minimum(couts_positifs))

            directions = [(-1,0), (0,1), (1,0), (0,-1)]

            # g(v) : coût réel du chemin optimal de vD à v dans la partie explorée
            g = fill(Inf, rows, cols)
            g[vD[1], vD[2]] = 0

            visites      = falses(rows, cols)
            predecesseur = fill((0,0), rows, cols)

            # Tas min : (f(v), pos) pour accepter plusieurs fois la même case (priorités différentes)
            heap = BinaryMinHeap{Tuple{Float64, Tuple{Int,Int}}}()
            hD = heuristique(vD, vA, c_min)
            push!(heap, (0.0 + hD, vD))

            nb_etats = 0
            cheminValide = false

            while !isempty(heap)
                (f_val, pos) = pop!(heap)
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

                    if 1 <= ni <= rows && 1 <= nj <= cols && M[ni, nj] != 0
                        new_g = g[i, j] + M[ni, nj]
                        if new_g < g[ni, nj]
                            g[ni, nj] = new_g
                            predecesseur[ni, nj] = pos
                            h_voisin = heuristique((ni, nj), vA, c_min)
                            push!(heap, (new_g + h_voisin, (ni, nj)))
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

            distance = cheminValide ? g[vA[1], vA[2]] : 0
        end
    end

    return afficherResultat(distance, nb_etats, chemin, temps, allocation)
end
