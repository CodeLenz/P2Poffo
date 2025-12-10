#
# Rotina que cria a malha e salva em um .yaml
#


"""
nome - nome que vai ser salvo o arquivo yaml
mesh - arquivo.portic com o dominio da malha do portico (existe função que cria "Exporta_1d")
secao - arquivo.geo ou arquivo.msh com os dados do Gmsh
mat - arquivo.mat com os dados do material (existe função que cria "Exporta_mat")
ap - arquivo.ap com as coordenadas dos nós que queremos restrição(existe função que cria "Exporta_apoios")
"""


function criayaml(nome::String, mesh::AbstractString, secao::AbstractString, mat::AbstractString, ap::AbstractString)

    #   Ler o arquivo do portico
    Portic_mesh = readlines(mesh)

    # Converter linhas para números
    Lx      =  parse(Float64, Portic_mesh[2])
    nx      =  parse(Int,     Portic_mesh[3])
    Ly      =  parse(Float64, Portic_mesh[4])
    ny      =  parse(Int,     Portic_mesh[5])
    Lz      =  parse(Float64, Portic_mesh[6])
    nz      =  parse(Int,     Portic_mesh[7])
    origin  = eval(Meta.parse(Portic_mesh[8]))


    # Criação da malha usando o BMesh
    bmesh = Bmesh_truss_3D(Lx,nx,Ly,ny,Lz,nz,origin=origin) 
    Lgmsh_export_init(nome*".pos",bmesh)

    # Recupera os dados da Malha
    ne      = bmesh.ne
    nn      = bmesh.nn
    coord   = bmesh.coord
    connect = bmesh.connect

    # Seção 
    geo = Exporta_sec(secao)
    geo_sec = geo * ".sec"

    #   Ler o arquivo da seção
    linha = readlines(geo_sec)

    # Converter linhas para números
    A    = parse(Float64, linha[2])
    Iz   = parse(Float64, linha[3])
    Iy   = parse(Float64, linha[4])
    J    = parse(Float64, linha[5])
    α    = parse(Float64, linha[6])

    # Recupera os dados do material (todas as linhas)
    material = readlines(mat)

    # Converter linhas para valores
    nome_mat = strip(material[1])
    Ex       = parse(Float64, strip(material[2]))
    G        = parse(Float64, strip(material[3]))
    S_esc    = parse(Float64, strip(material[4]))

    
    # Ler todas as linhas
    linhas = readlines(ap)

    coord_apoios = []
    gdl = Int[]
    valor = Int[]

    for lin in linhas
        # Divide a linha por espaços
        cols = split(lin)

        # Converte cada coluna
        x = parse(Float64, cols[1])
        y = parse(Float64, cols[2])
        z = parse(Float64, cols[3])
        g = parse(Int, cols[4])
        v = parse(Int, cols[5])

        push!(coord_apoios, (x, y, z))
        push!(gdl, g)
        push!(valor, v)
    end 

    
    nos = Int[]
    # loop pelas quantidade de coordenadas 
    for i in axes(coord_apoios, 1)

        # coordenadas do arquivo
        x, y, z = coord_apoios[i]

        #encontra o nó mais perto
        no = Find_node(bmesh, x, y, z)

        push!(nos, no)
    end

    # Nome que será salvo o yaml
    arquivo = string(nome, ".yaml")

    open(arquivo, "w") do io

        # Materiais 
        write(io, "materiais:\n")
        write(io, """  - nome: \"$(nome_mat)\"\n""")
        write(io,   "    G: $(G)\n")
        write(io,   "    Ex: $(Ex)\n")
        write(io,   "    S_esc: $(S_esc)\n\n")

        # Carregamentos (vazio por enquanto)
        write(io, "floads:\n")
        write(io, "  # i FX FY FZ MX MY MZ\n\n")

        write(io, "loads:\n")
        write(io, " \n")

        write(io, "geometrias:\n")
        write(io, "  - Iz: $(Iz)\n")
        write(io, "    A: $(A)\n")
        write(io, "    Iy: $(Iy)\n")
        write(io, "    nome: $(geo)\n")
        write(io, "    α: $α\n")
        write(io, "    J0: $(J)\n\n")

        # Coordenadas da malha
        write(io, "coordenadas:\n")
        for i in 1:nn
            x, y, z = coord[i, :]
            write(io, "  $(x) $(y) $(z)\n")
        end
        write(io, "\n")

        # Conectividades
        write(io, "conectividades:\n")
        for e in 1:ne
            n1, n2 = connect[e, :]
            write(io, "  $(n1) $(n2)\n")
        end
        write(io, "\n")

        # Apoios 
        write(io, "apoios:\n")

        for i in eachindex(nos)
            write(io, "  $(nos[i]) $(gdl[i]) $(valor[i])\n")
        end
        

        # Dados dos elementos
        write(io, "dados_elementos:\n")
        write(io, "  $nome_mat $(geo) \n\n")
    

    end

    println("Arquivo salvo: $arquivo")

    return bmesh
end