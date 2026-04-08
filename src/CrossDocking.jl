
include("lireDonneeAMR.jl")
include("MapEnGraphe.jl")
include("Astar2.0.jl")

function crossDocking(fname)
    AMR_dict = lireDonneeAMR("dat/$fname")
    M = lireMapPondere("dat/$fname")
    println(M)
    println(AMR_dict)

    directionsGaucheDroite = [(-1,0), (1,0), (0,1), (0,-1)]
    directionsHautBas = [(0,1), (0,-1), (-1,0), (1,0)]

    rows = size(M, 1)
    cols = size(M, 2)
    
    s = "" # string pour les résultats

    for id in keys(AMR_dict)
        chemin = algoAstar2_0(fname, AMR_dict[id].actuelle_pos, AMR_dict[id].arrivee)
        AMR_dict[id].chemin = chemin
        M[AMR_dict[id].actuelle_pos[1], AMR_dict[id].actuelle_pos[2]] = 0
        println("chemin de $id : $chemin")
    end

    println(AMR_dict)
    t = 0


    while !isempty(AMR_dict)

        #sleep(1)
        println("--------------------------------")
        println("Temps : $t")

        for id in sort(collect(keys(AMR_dict))) # pour tous les AMR
            if AMR_dict[id].t == t # on regarde si l'AMR est au bon temps

                println("id : $id")
                chemin_suivant = dequeue!(AMR_dict[id].chemin)
                #println(chemin_suivant)

                for id2 in sort(collect(keys(AMR_dict))) # pour tous les autres AMR
                    if id2 != id && AMR_dict[id2].actuelle_pos == chemin_suivant # on regarde le chemin est valide (pas de collision avec un autre AMR)
                    println("Collision entre AMR $id et AMR $id2")


                        pq = PriorityQueue{Tuple{Int,Int}, Int64}()

                        x = AMR_dict[id].actuelle_pos[1] - chemin_suivant[1]

                        if x == 1 || x == -1
                            directions = directionsHautBas # pour regarder les cases au dessus et en dessous avant de regarder les cases à gauche et à droite
                        else
                            directions = directionsGaucheDroite # pour regarder les cases à gauche et à droite avant de regarder les cases au dessus et en dessous
                        end

                        for (dl, dc) in directions # deplacer l'AMR qui a causé la collision d'une case
                            cout = M[AMR_dict[id2].actuelle_pos[1], AMR_dict[id2].actuelle_pos[2]]

                            ni = AMR_dict[id].actuelle_pos[1] + dl
                            nj = AMR_dict[id].actuelle_pos[2] + dc

                            M[AMR_dict[id2].actuelle_pos[1], AMR_dict[id2].actuelle_pos[2]] = 0
                            # Reste dans la Map et pas un mur
                            if 1 <= ni <= rows && 1 <= nj <= cols && M[ni, nj] != 0 
                                enqueue!(pq, (ni, nj) => M[ni, nj]) 
                            end

                            M[AMR_dict[id2].actuelle_pos[1], AMR_dict[id2].actuelle_pos[2]] = cout
                        end

                        println(pq)

                        if isempty(pq)
                            println("Aucun déplacement possible pour AMR $id")
                            # aucun déplacement possible (bloqué) on passe au prochain AMR
                            chemin_suivant = AMR_dict[id].actuelle_pos
                            AMR_dict[id].chemin = algoAstar2_0(fname, chemin_suivant, AMR_dict[id].arrivee)
                            continue
                        end
                    
                        chemin_suivant = dequeue!(pq)  # renvoie la position (Tuple), pas un Pair

                        AMR_dict[id].chemin = algoAstar2_0(fname, chemin_suivant, AMR_dict[id].arrivee)

                        break
                    end
                end
                M[AMR_dict[id].actuelle_pos[1], AMR_dict[id].actuelle_pos[2]] = 1
                AMR_dict[id].actuelle_pos = chemin_suivant
                println("chemin de $id : $(AMR_dict[id].chemin)")
                M[AMR_dict[id].actuelle_pos[1], AMR_dict[id].actuelle_pos[2]] = 0

                println("position de $id : $(AMR_dict[id].actuelle_pos)")
                AMR_dict[id].t += 1

            end 
        end

        t += 1

        for id in sort(collect(keys(AMR_dict)))
            if AMR_dict[id].actuelle_pos == AMR_dict[id].arrivee
                println("AMR $id a atteint son objectif")
                s = s * ("AMR $id a atteint son objectif apres $(AMR_dict[id].t) iterations") * "\n"
                delete!(AMR_dict, id)
            end
        end

    end
    println("--------------------------------")
    println(s)
end