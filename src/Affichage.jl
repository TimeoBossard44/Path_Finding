function afficherResultat( distance, nb_etats, chemin, temps, allocation, fname,s::String)

    println("Distance D -> A : ", distance)
    println("Number of states evaluated : ", nb_etats)

    if isempty(chemin)
        println("Path D -> A : no path found")
    else
        println("Path D -> A : ", join(chemin, "→"))
        outmap = "res/" * splitext(fname)[1] * "_$(s).map"
        ecrire_solution_map("dat/$fname", chemin, outmap)
        map_to_pdf(outmap,s)
    end

    println("Temps d'execution : ", temps*1000 , " ms")
    println("Nombre d'allocation : ", allocation )
    
end
