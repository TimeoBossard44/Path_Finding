
function afficherResultat( distance, nb_etats, chemin, temps, allocation)

    println("Distance D -> A : ", distance)
    println("Number of states evaluated : ", nb_etats)

    if isempty(chemin)
        println("Path D -> A : no path found")
    else
        println("Path D -> A : ", join(chemin, "→"))
    end

    println("Temps d'execution : ", temps*1000 , " ms")
    println("Nombre d'allocation : ", allocation )
    
end
