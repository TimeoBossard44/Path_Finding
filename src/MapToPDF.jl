using CairoMakie

_basename_noext(p::AbstractString) = splitext(basename(p))[1]

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
            elseif char == 'X'
                6   # chemin
            else
                0   # inconnu
            end
        end
    end
    
    return grid
end

function map_to_pdf(input_file::String, s::String)
    # `input_file` peut être "theglaive.map", "dat/theglaive.map", "res/theglaive.map", etc.
    # On ne préfixe pas automatiquement par "dat/" pour éviter les chemins du type "dat/res/...".
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
        RGBf(0.9, 0.1, 0.1),   # 6 - chemin 'X' → rouge
    ]
    
    cmap = cgrad(colors, 7, categorical = true)

    fig = Figure(size = (900, 900))
    ax = Axis(fig[1, 1],
        title = input_file,
        aspect = DataAspect()
    )

    # Orientation : certains .map sont "couchés" avec permutedims seulement.
    # On applique une rotation 90° pour remettre la carte dans le bon sens.
    heatmap!(ax, rotr90(grid), colormap = cmap, colorrange = (0, 6))


    # `s` peut être un préfixe de chemin (ex: \"solutionBFS_\") ou un dossier (ex: \"solutionBFS_res\")
    # On construit un chemin de sortie valide et on crée le dossier si nécessaire.
    base = _basename_noext(input_file) * ".pdf"
    out = joinpath(s, base)

    mkpath(dirname("res/"*out))

    save("res/"*out, fig)
    println("PDF généré : res/"*out)
end