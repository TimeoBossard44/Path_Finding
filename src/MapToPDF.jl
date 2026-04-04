using CairoMakie

function parse_map(filename::String)
    lines = readlines(filename)
    
    # Lire l'en-tête
    height = parse(Int, split(lines[2])[2])
    width  = parse(Int, split(lines[3])[2])
    
    # Lire la carte (après "map")
    map_lines = lines[5:end]
    
    grid = zeros(Int, height, width)
    
    for (i, line) in enumerate(map_lines)
        for (j, char) in enumerate(line)
            grid[i, j] = if char == '@'
                1   # mur
            elseif char == '.'
                2   # sol
            elseif char == 'T'
                3   # terre
            elseif char == 'S'
                4   # sable
            elseif char == 'W'
                5   # eau
            else
                0   # inconnu
            end
        end
    end
    
    return grid
end

function map_to_pdf(input_file::String, output_file::String)
    grid = parse_map(input_file)

    # Couleurs associées à chaque valeur
    # 0=inconnu, 1=mur, 2=sol, 3=terre, 4=sable, 5=eau
    colors = [
        RGBf(0.2, 0.2, 0.2),   # 0 - inconnu   → gris foncé
        RGBf(0.1, 0.1, 0.1),   # 1 - mur   '@' → noir
        RGBf(0.9, 0.9, 0.9),   # 2 - sol   '.' → blanc cassé
        RGBf(0.6, 0.4, 0.2),   # 3 - terre 'T' → marron
        RGBf(0.9, 0.8, 0.5),   # 4 - sable 'S' → jaune sable
        RGBf(0.2, 0.5, 0.9),   # 5 - eau   'W' → bleu
    ]
    
    cmap = cgrad(colors, 6, categorical = true)

    fig = Figure(size = (900, 900))
    ax = Axis(fig[1, 1],
        title = input_file,
        aspect = DataAspect()
    )

    heatmap!(ax, permutedims(grid), colormap = cmap, colorrange = (0, 5))


    save(output_file, fig)
    println("PDF généré : $output_file")
end