#
# Propriedades de um círculo de raio 1cm
#
@testset "Circulo" begin

    # Valor de referência 
    area, Iz, Jeq = P2Poffo.Circulo(1E-2)

    # Vamos utilizar o .msh de referência
    arquivo = "circular.msh"

    # Calcula com o programa
    centroide, area, Izl, Iyl, Jeq, α = Pre_processamento(arquivo)


end
