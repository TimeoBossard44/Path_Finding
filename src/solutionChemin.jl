"""
Écrit un fichier `.map` modifié où toutes les positions du `chemin` sont remplacées par 'X'.

- `input_map` : chemin du fichier `.map` source
- `chemin`    : vecteur de positions (ligne, colonne) = `Vector{Tuple{Int,Int}}`
- `output_map`: chemin du fichier `.map` destination

La modification ne touche que la section après la ligne `map`.
"""
function ecrire_solution_map(input_map::String, chemin::AbstractVector{<:Tuple{Int,Int}}, output_map::String)
    lines = readlines(input_map)

    i_map = findfirst(l -> strip(l) == "map", lines)
    i_map === nothing && error("ligne 'map' introuvable dans $input_map")

    map_start = i_map + 1

    # Pour accès O(1) : ensemble de positions à marquer
    to_mark = Set{Tuple{Int,Int}}(chemin)

    for (idx, rawline) in enumerate(lines[map_start:end])
        row = idx
        line = chomp(rawline)
        isempty(line) && continue

        chars = collect(line)
        for col in 1:length(chars)
            if (row, col) in to_mark
                chars[col] = 'X'
            end
        end

        lines[map_start + idx - 1] = String(chars)
    end

    open(output_map, "w") do io
        for l in lines
            println(io, l)
        end
    end

    return output_map
end

