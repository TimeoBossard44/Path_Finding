mutable struct AMR
    actuelle_pos::Tuple{Int,Int}
    chemin::Queue{Tuple{Int,Int}}
    arrivee::Tuple{Int,Int}
    t::Int
end

function lireDonneeAMR(fichier)
    AMR_dict = Dict{Int, AMR}()
    lecture_active = false    # on commence sans lire la map

    open(fichier, "r") do f # lire le fichier et le refermer 
        for ligne in eachline(f)

            # Dès qu'on voit "amr", on active la lecture
            if strip(ligne) == "amr"
                lecture_active = true
                continue       
            end

            # On ne traite les lignes que si on a vu "amr" (jusqu'à la section « map »)
            if lecture_active && !isempty(strip(ligne))
                s = strip(ligne)
                s == "map" && break
                donneeAMR = split(s)
                length(donneeAMR) == 4 || error("ligne AMR invalide (attendu: id x,y x,y t) : $s")
                id, actuelle_pos, arrivee, t = donneeAMR

                id = parse(Int, id)

                a,b = split(actuelle_pos, ',')
                actuelle_pos = (parse(Int, a), parse(Int, b))

                a,b = split(arrivee, ',')
                arrivee = (parse(Int, a), parse(Int, b))

                t = parse(Int, t)

                AMR_dict[id] = AMR(actuelle_pos, Queue{Tuple{Int,Int}}(), arrivee, t)
            end
        end
    end

    return AMR_dict
end