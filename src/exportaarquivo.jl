#
# Exporta um arquivo em .sec com a ordem experada do LFrame
#

# Ordem que precisa ser salvo
#*.sec
# nome
# A
# Iz 
# Iy
# J0
# α
function Exporta_sec(arquivo, gera_pos=false)

    # muda a terminação do arquivo para .sec
    if occursin(".msh",arquivo) 
         
        nome_sec = replace(arquivo,".msh"=>".sec")

    else 
        nome_sec = replace(arquivo,".geo"=>".sec")
        
    end

    # Nome da geometria
    nome_ = basename(nome_sec) 
    geo = nome_[1:end-4]
    
    # Obtém os valores da análise
    centroide, area, Izl, Iyl, Jeq, α, _  = Pre_processamento(arquivo, gera_pos)

    # Abre o arquivo para escrita e grava as informações uma por linha
    open(nome_sec, "w") do file
        println(file, geo)
        println(file, area)
        println(file, Izl)
        println(file, Iyl)
        println(file, Jeq)
        println(file, α)
    end

    # Fim da rotina
    return geo

end


# nome      = Nome que vai ser salvo o arquivo.yaml
# Lx        = Tamanho em x
# nx        = numero de divisões em x
# Ly        = Tamanho em y
# ny        = numero de divisões em y
# Lz        = Tamanho em z
# nz        = numero de divisões em z
# origin    = coordenadas do nó de origem (nó 1)

function Exporta_1d(nome::String,Lx::Float64,nx::Int64,Ly::Float64,ny::Int64,Lz::Float64,nz::Int64; origin=(0.0,0.0,0.0))

    # Nome do arquivo com a malha do portico
    nome_port = nome*".portic"

    # Abre o arquivo para escrita e grava as informações uma por linha
    open(nome_port, "w") do file
        println(file, nome)
        println(file, Lx)
        println(file, nx)
        println(file, Ly)
        println(file, ny)
        println(file, Lz)
        println(file, nz)
        println(file, origin)
    end

    # Fim da rotina
    return nome

end

function Exporta_mat(nome::String,Ex::Float64,G::Float64,S_esc::Float64)

    # Nome do arquivo com a malha do portico
    nome_mat = nome*".mat"

    # Abre o arquivo para escrita e grava as informações uma por linha
    open(nome_mat, "w") do file
        println(file, nome)
        println(file, Ex)
        println(file, G)
        println(file, S_esc)
    end

    # Fim da rotina
    return nome

end