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
function ExportaDat(arquivo, gera_pos=false)

    # Recupera o nome do arquivo sem caminho
    caminho = pathof(P2Poffo)[1:end-14]*"\\geometria"

    # Nome do arquivo .msh
    nome_msh = basename(arquivo)

    # repassa o nome do arquivo para .sec 
    nome_sec = replace(nome_msh,".msh"=>".sec")

    # Cria o arquivo completo do .pos com o nome do yaml
    nome = joinpath(caminho, nome_sec)
 
    geo = nome_msh[1:end-4]
    # Obtém os valores da análise
    centroide, area, Izl, Iyl, Jeq, α, _  = Pre_processamento(arquivo, gera_pos)

    # Abre o arquivo para escrita e grava as informações uma por linha
    open(nome, "w") do file
        println(file, geo)
        println(file, area)
        println(file, Izl)
        println(file, Iyl)
        println(file, Jeq)
        println(file, α)
    end
end