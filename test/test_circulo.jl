#
# Propriedades de um círculo de raio 1cm
#
@testset "Circulo" begin

    # Valor de referência 
    area, Iz, Jeq = P2Poffo.Circulo(1E-2)

    caminho = pathof(P2Poffo)[1:end-14]*"\\test"
    arquivo = joinpath([caminho,"data","quadrado.msh"])

    # Calcula com o programa
    centroide, area, Izl, Iyl, Jeq, α = Pre_processamento(arquivo)


end
