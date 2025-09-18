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
function Exporta(arquivo, gera_pos=false)

    # muda a terminação do arquivo para .sec
    if occursin(".msh",arquivo) 
         
        nome_sec = replace(arquivo,".msh"=>".sec")

    else 
        nome_sec = replace(arquivo,".geo"=>".sec")
        
    end

    # Nome da geometria 
    geo = nome_[1:end-4]
    
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

    # Fim da rotina
    return nothing

end