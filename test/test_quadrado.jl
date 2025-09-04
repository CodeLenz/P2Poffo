#
# Propriedades de um quadrado de lado 1cm
#
@testset "Quadrado" begin

    # Valor de referência 
    area, Iz, Jeq = P2Poffo.Quadrado(1E-2)

    # Vamos utilizar o .msh de referência
    caminho = pathof(P2Poffo)[1:end-14]*"\\test"
    arquivo = joinpath([caminho,"data","quadrado.msh"])

    # Calcula com o programa
    centroide, area, Izl, Iyl, Jeq, α = Pre_processamento(arquivo)


end
