#
# Rotinhas de integração
#

# V(ρ) <= V_sup = vf * V0
# 

function Criterio_Otimo(arquivo::AbstractString,posfile=true; verbose=false,ρ0=[],vf = 0.5, μ1=0.0, μ2=1E5, δ=0.1, tol=1E-6,ρ_min = 1E-3)

    # Calcula os deslocamentos e propriedades da malha inicial
    U, malha = LFrame.Analise3D(arquivo,posfile,verbose=verbose,ρ0=ρ0)

    # Calcula o volume inicial da malha    
    V = volumes(malha)

    # Volume total ta estrutura
    V0 = sum(V)

    # Volume máximo
    V_sup = vf * V0

    # Copia do ρ0
    ρ_estimado = copy(ρ0)

    # iniciando as derivadas
    dC = zeros(malha.ne)
    dV = zeros(malha.ne)


    # Loop da bisseção
    for k=1:1000  
       
       # Ponto médio do Intervalo 
       μ = (μ1 + μ2)/2
          
        # Loop pelos elementos
        for i in malha.ne     

            # Derivada do compliance em relação ao ρ
            dC[i] = dCompliance(malha, U, ρ_estimado)[i]

            # Derivada do volume em relação ao ρ
            dV[i] = V[i]

            # fator beta(aqui tinha um min)
            β = -dC[i]/(μ*dV[i])
        
            # Novo valor ρ estimado 
            ρ_e = ρ0[i] * (β^0.5)
            
            # Limites
            ρ_dir = min(ρ0[i] + δ, 1.0) 
            ρ_esq = max(ρ0[i] - δ, ρ_min) 
           
            # 
            #-------0[----|--x--|----]1---------  
            #           x-δ    x+δ 
            #
            
            # Verificando os Limites
            ρ_estimado[i] = max( min( ρ_e, ρ_dir), ρ_esq)
        end 

       # Volume novo
       V_atual = sum(ρ_estimado .* V)
 
       # Testando se atingiu a tolerancia
       if abs(μ2 - μ1)<=tol
            break
       end
        
       # Verificando se o volume atingiu o máximo
       if (V_atual > V_sup)
            μ1 =  μ
       else
            μ2 =  μ
       end
 
        
    end # k
 
    # Return o ρ_estimado
    return ρ_estimado
    
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

    # Retorna o volume
    return V
end


function dCompliance(malha,U::Vector,ρ::Vector)    

    # iniciando a derivada
    D = zeros(malha.ne)

    # Loop pelos elementos
    for ele in malha.ne
        dados_ele = malha.dados_elementos
        dicionario_mat = malha.dicionario_materiais
        dicionario_geo = malha.dicionario_geometrias

        Iz, Iy, J0, A, α, E, G, geo = LFrame.Dados_fundamentais(ele, dados_ele, dicionario_mat, dicionario_geo)

        # gdls - 6 graus por nó
        dofs = malha.nnos * 6   

        # deslocamentos no global
        ug = U[dofs] 

        # Recupera dados da malha
        ij = malha.conect
        coordenadas = malha.coord

        # matriz de tranformação de coord(ver se é só rotacao ou translação)
        T = LFrame.Rotacao3d(ele,ij,coordenadas,α)

        # Para o local
        ul = T'ug
        
        # Rigidez
        Ke = LFrame.Ke_portico3d(E,Iz,Iy,G,J0,malha.L[ele],A)

        # Vetor das derivada
        D[ele]  =  -dot(ul,Ke,ul)

    end
   
    return D

end
