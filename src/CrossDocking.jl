#= creer un tuple (depart, arrivee, temps) pour chaque AMR
  calculer avec Astar le plus court chemin de chaque AMR
  ajouter un compteur de pas pour chaque AMR et regarder si a un pas k le chemin est toujours 
  valide (pas de collision avec un autre AMR) si il n est plus valide alors 
  le recalcluer a partir de la position de l'AMR qui a causé la collision et je le deplace
  d une case pour ne plus avoir de collision quand le deuxieme AMR devra se déplacer.

  un AMR regarder si son chemin est valide et si plus valide alors on recalcule son chemin 
  avec la condition de rajouter un mur temporaire a la position de l'AMR qui a causé la collision
  sur le fichier.map et apres enlever le mur temporaire.

=#
# structure des AMR pour les utiliser dans un dictionnaire 

# a chaque t+1 on fait un popfirst! sur le chemin
# on creer un dictionnaire pour les AMR avec comme clé des int de 1 a n (n nombre d'AMR) et comme valeur l'AMR
# pour regarder si il y a  une collision on va regarder si pour un amr la nouvelle position 
# est la meme que pour un autre amr et si c est le cas alors on recalcule le chemin de l'AMR qui a causé la collision


# il faut cas chaque iteration +1 tous les t dans le dictionnaire on le meme t a la fin de l'iteration

include("lireDonneeAMR.jl")
include("MapEnGraphe.jl")
include("Astar2.0.jl")

function crossDocking(fname)
    AMR_dict = lireDonneeAMR("dat/$fname")
    M = lireMapPondere("dat/$fname")
    println(M)
    println(AMR_dict)

    directions = [(-1,0), (0,1), (1,0), (0,-1)]

    rows = size(M, 1)
    cols = size(M, 2)
    

    for id in keys(AMR_dict)
        chemin = algoAstar2_0(fname, AMR_dict[id].actuelle_pos, AMR_dict[id].arrivee)
        AMR_dict[id].chemin = chemin
        println("chemin de $id : $chemin")
    end

    println(AMR_dict)
    t = 0


    while !isempty(AMR_dict)

        #sleep(1)
        #println("Temps : $t")

        for id in sort(collect(keys(AMR_dict))) # pour tous les AMR
            if AMR_dict[id].t == t # on regarde si l'AMR est au bon temps

                #println(id)
                chemin_suivant = dequeue!(AMR_dict[id].chemin)
                #println(chemin_suivant)

                for id2 in sort(collect(keys(AMR_dict))) # pour tous les autres AMR
                    if id2 != id && AMR_dict[id2].actuelle_pos == chemin_suivant # on regarde le chemin est valide (pas de collision avec un autre AMR)
                    #println("Collision entre AMR $id et AMR $id2")


                        pq = PriorityQueue{Tuple{Int,Int}, Int64}()

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

                        #println(pq)

                        if isempty(pq)
                            #println("Aucun déplacement possible pour AMR $id")
                            # aucun déplacement possible (bloqué) → on passe au prochain AMR
                            chemin_suivant = AMR_dict[id].actuelle_pos
                            AMR_dict[id].chemin = algoAstar2_0(fname, chemin_suivant, AMR_dict[id].arrivee)
                            continue
                        end
                    
                        chemin_suivant = dequeue!(pq)  # renvoie la position (Tuple), pas un Pair

                        AMR_dict[id].chemin = algoAstar2_0(fname, chemin_suivant, AMR_dict[id].arrivee)

                        if isempty(AMR_dict[id].chemin) 
                            #println("Aucun chemin trouver pour AMR $id")
                            continue
                        end
                        break
                    end
                end

                AMR_dict[id].actuelle_pos = chemin_suivant
                #println("chemin de $id : $(AMR_dict[id].chemin)")

                #println("position de $id : $(AMR_dict[id].actuelle_pos)")
                AMR_dict[id].t += 1

            end 
        end

        t += 1

        for id in sort(collect(keys(AMR_dict)))
            if AMR_dict[id].actuelle_pos == AMR_dict[id].arrivee
                println("AMR $id a atteint son objectif apres $(AMR_dict[id].t) iterations")
                delete!(AMR_dict, id)
            end
        end

    end
end