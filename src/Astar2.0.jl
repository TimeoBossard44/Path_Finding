include("MapEnGraphe.jl")
using DataStructures

function heuristique(v, vA)
    return (abs(v[1] - vA[1]) + abs(v[2] - vA[2]))
end

function algoAstar2_0(fname, vD, vA)

    # Matrice pondérée du fichier fname
    M = lireMapPondere("dat/$fname")


    rows = size(M, 1)
    cols = size(M, 2)

    directions = [(-1,0), (0,1), (1,0), (0,-1)]

    dist = fill(Inf, rows, cols)
    dist[vD[1], vD[2]] = 0

    visites      = falses(rows, cols)
    predecesseur = fill((0,0), rows, cols)

    pq = PriorityQueue{Tuple{Int,Int}, Int64}() # File de priorité pour avoir accées à la plus court distance en premier
    enqueue!(pq, vD => 0) 

    cheminValide = false

    while !isempty(pq)
        pos = dequeue!(pq)
        i, j = pos

        if visites[i, j]
            continue
        end

        visites[i, j] = true

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
                newh = heuristique((ni,nj),vA)
                newvalue = newdist + newh
                if newvalue < dist[ni, nj] + newh
                    dist[ni, nj] = newdist
                    predecesseur[ni, nj] = pos
                    # PriorityQueue: une même position (clé) ne peut apparaître qu'une fois.
                    # Si elle est déjà dedans, on met à jour la priorité si elle s'améliore.
                    if haskey(pq, (ni, nj))
                        if newvalue < pq[(ni, nj)]
                            pq[(ni, nj)] = newvalue
                        end
                    else
                        enqueue!(pq, (ni, nj) => newvalue)
                    end
                end
            end
        end

    end

    # Reconstruction du chemin
    temp_chemin = Vector{Tuple{Int,Int}}()
    if cheminValide
        etape = vA
        while etape !== (0,0)
            push!(temp_chemin, etape)
            etape = predecesseur[etape[1], etape[2]]
        end
    end

    chemin = Queue{Tuple{Int,Int}}()
    for i in length(temp_chemin):-1:1
        enqueue!(chemin, temp_chemin[i])
    end
    dequeue!(chemin) #on retire la position de départ

    return chemin
end
