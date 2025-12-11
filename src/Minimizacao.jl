#
# Rotinas de integração
#
#
# V(ρ) <= V_sup = vf * V0
# 

# Parametrização SIMP 
simp(ρ,p=3) = ρ^p

# Derivada da parametrização 
dsimp(ρ,p=3) = p*ρ^(p-1)

#
# Rotina principal
#
function Main_Otim(arquivo::AbstractString, fkparam::Function, fdkparam::Function, posfile=true; verbose=false,vf = 0.5, niter=100)

    # Chama o Analise3D com o nome do arquivo, para receber a estrutura de malha
    U0, malha = Analise3D(arquivo,posfile,verbose=verbose)

    # Numero de elementos 
    ne = malha.ne

    # Estimativa inicial das densidades relativas
    ρ0 = vf*ones(ne)

    # Calcula o volume de cada elemento sem considerar a 
    # parametrização 
    V = Volumes(malha)

    # Volume total ta estrutura
    V0 = sum(V)

    # Volume limite
    V_sup = vf * V0
    
    # Derivada do volume é fixa 
    dV = V
    
    # Loop externo de otimização 
    for iter=1:niter

        # Calcula os deslocamentos
        U, _ = Analise3D(malha,posfile,ρ0=ρ0,kparam=[fkparam])
       
        @show iter, sum(ρ0.*dV), U[8]

        # Deriva da compliance
        dC = dCompliance(malha,U,ρ0,fdkparam)    

        # Atualiza as densidades relativas utilizando o OC
        ρ0 .= OC(ρ0,dC,dV,V_sup,ne)
        
    end # loop externo
    
    # Adiciona uma vista com as densidades relativas ao arquivo de saída
    if posfile

        # Nome do arquivo com .pos 
        arquivo_pos =  replace(arquivo,".yaml"=>".pos")
        
        # Exporta a vista escalar 
        Lgmsh.Lgmsh_export_element_scalar(arquivo_pos,ρ0,"ρ")

    end

   # Retorna as densidades relativas 
   return ρ0

end

#
# Retorna um vetor com o volume de cada elemento da malha
#
function Volumes(malha::LFrame.Malha)

    # Recupera os dados da malha
    L = malha.L 
    dados_ele = malha.dados_elementos
    dicionario_mat = malha.dicionario_materiais
    dicionario_geo = malha.dicionario_geometrias


    # Inicializa o vetor de volumes
    V = zeros(malha.ne)

    # Loop pelos elementos
    for ele = 1:malha.ne

        Iz, Iy, J0, A, α, E, G, geo = LFrame.Dados_fundamentais(ele, dados_ele, dicionario_mat, dicionario_geo)

        V[ele] = A * L[ele]
    end

    # Retorna o volume
    return V

end


#
# Critério de ótimo (loop local)
#
function OC(ρ0::Vector, dC::Vector, dV::Vector, V_sup::Float64, ne::Int64, μ1=0.0, μ2=1E5, δ=0.1, tol=1E-6,ρ_min = 1E-3)

    # Copia do ρ0
    ρ_estimado = copy(ρ0)

    # Aqui, V e dV são iguais
    V = dV

    # Loop da bisseção
    for k=1:1000  
       
       # Ponto médio do Intervalo 
       μ = (μ1 + μ2)/2
          
        # Loop pelos elementos
        for i=1:ne     
            
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

function dCompliance(malha::LFrame.Malha,U::Vector,ρ::Vector,fdkparam::Function)    

    # Iniciando a derivada
    D = zeros(malha.ne)

    # Alias dos dados da estrutura de malha
    dados_ele = malha.dados_elementos
    dicionario_mat = malha.dicionario_materiais
    dicionario_geo = malha.dicionario_geometrias
    conect = malha.conect
    coordenadas = malha.coord
    L = malha.L

    # Loop pelos elementos
    for ele in 1:malha.ne
        
        # Recupera dados do elemento 
        Iz, Iy, J0, A, α, E, G, geo = LFrame.Dados_fundamentais(ele, dados_ele, dicionario_mat, dicionario_geo)

        # Graus de liberdade do elemento 
        dofs = LFrame.Gls(ele,conect)

        # deslocamentos no global
        ug = U[dofs] 
       
        # matriz de tranformação de coord(ver se é só rotacao ou translação)
        T = LFrame.Rotacao3d(ele,conect,coordenadas,α)

        # Para o local
        ul = T'*ug
        
        # Rigidez local do elemento 
        Ke = LFrame.Ke_portico3d(E,Iz,Iy,G,J0,L[ele],A)

        # Deriada da compliance C em relação ao elemento
        D[ele]  =  -dot(ul,Ke,ul)*fdkparam(ρ[ele])

    end
   
    return D

end
