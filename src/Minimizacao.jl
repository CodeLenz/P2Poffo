#
# Rotinhas de integração
#

function Criterio_Otimo(arquivo::AbstractString,posfile=true; verbose=false,ρ0=[])

    # Calcula os deslocamentos e propriedades da malha inicial
    U, malha = LFrame.Analise3D(arquivo,posfile,verbose=verbose,ρ0=ρ0)

    # Calcula o volume inicial da malha    
    V = volumes(malha)

    # Volume total ta estrutura
    V0 = sum(V)

end

function volumes(malha)

    # Recupera os dados da malha
    L = malha.L 
    dados_ele = malha.dados_elementos
    dicionario_mat = malha.dicionario_materiais
    dicionario_geo = malha.dicionario_geometrias


    V = zeros(malha.ne)
    # Loop pelos elementos
    for ele = 1:malha.ne

        Iz, Iy, J0, A, α, E, G, geo = LFrame.Dados_fundamentais(ele, dados_ele, dicionario_mat, dicionario_geo)

        V[ele] = A * L[ele]
    end

    return V
end