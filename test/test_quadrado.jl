#
# Propriedades de um quadrado de lado 1cm
#
@testset "Quadrado" begin

    # Valor de referência 
    area, Iz, Jeq = Torcao.Quadrado(1E-2)

    # Vamos utilizar o .msh de referência
    arquivo = "qudrado.msh"

    # Calcula com o programa
    centroide, area, Izl, Iyl, Jeq, α = AnaliseTorcao(arquivo)


end
