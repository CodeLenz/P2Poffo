#
# Exporta um arquivo em .dat com a ordem experada do LFrame
#

# Ordem que precisa ser salvo
#*.dat
# nome
# A
# Iz 
# Iy
# J0
# α
function ExportaDat(arquivo)

    # Recupera o nome do arquivo sem caminho
    nome_msh = basename(arquivo)

    # tira a extensão .msh
    nome = nome_msh[1:end-4]   

    # Obtém os valores da análise
    centroide, area, Izl, Iyl, Jeq, α  = AnaliseTorcao(arquivo)

    # Define o nome do arquivo .dat
    nome_dat = string(nome, ".dat")

    # Abre o arquivo para escrita e grava as informações uma por linha
    open(nome_dat, "w") do file
        println(file, nome)
        println(file, area)
        println(file, Izl)
        println(file, Iyl)
        println(file, Jeq)
        println(file, α)
    end
end