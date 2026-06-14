#
# Rotinas que calcula as derivadas para otimização
#

#####################################################################################
                        # Derivada da Compliance
#####################################################################################
function dCompliance(malha::LFrame.Malha,U::Vector,x::Vector,fdkparam::Function)    

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
        D[ele]  =  -dot(ul,Ke,ul)*fdkparam(x[ele])

    end
   
    return D

end

#####################################################################################
                        # Derivada da frequencia
#####################################################################################
function dω(ωn::Float64,U0::Vector,malha::LFrame.Malha,x::Vector,fdkparam::Function,fdmparam::Function,fmparam::Function)    

    # Iniciando a derivada
    D = zeros(malha.ne)

    # Dados da estrutura de malha
    dados_ele = malha.dados_elementos
    dicionario_mat = malha.dicionario_materiais
    dicionario_geo = malha.dicionario_geometrias
    conect = malha.conect
    coordenadas = malha.coord
    apoios = malha.apoios
    L = malha.L
    nos = malha.nnos

    # Inicializa o vetor com os modos globais
    ϕg = zeros(nos*6)

    # Descobre quais gdls são livres
    dofs_l = LFrame.dofs_livres(nos,apoios)

    # Associa o gdl livre ao vetor de modos
    ϕg[dofs_l] .= U0

    Mg = LFrame.Monta_Mg(malha,x,fmparam)

    pm =ϕg'*Mg*ϕg

    # Loop pelos elementos
    for ele in LinearIndices(x)
        
        # Recupera dados do elemento 
        Iz, Iy, J0, A, α, E, G, geo, ρ = LFrame.Dados_fundamentais(ele, dados_ele, dicionario_mat, dicionario_geo)

        # Rigidez local do elemento 
        Ke = LFrame.Ke_portico3d(E,Iz,Iy,G,J0,L[ele],A)

        # Massa local do elemento
        Me = LFrame.Me_portico3d(A,ρ,J0,L[ele])

        # Derivada da rigidez em relação a variavel
        dK = Ke*fdkparam(x[ele])

        # Derivada da massa em relação a variavel
        dM = Me*fdmparam(x[ele])
         
        # Graus de liberdades globais do elemento
        dofs = LFrame.Gls(ele,conect)
        
        # Modos no global
        ϕ0g = ϕg[dofs]
       
        # Matriz de tranformação de coord(ver se é só rotacao ou translação)
        T = LFrame.Rotacao3d(ele,conect,coordenadas,α)
        
        # Para o local
        ϕ0l = T*ϕ0g
    
        # variavel auxiliar para facilitar
        aux = dK - (ωn^2)*dM

        # Derivada do elemento
        D[ele] = dot(ϕ0l,aux,ϕ0l)/(2*ωn)

    end

    # Deriadas da frequencia ωn em relação a todos os elemento
    d = D ./ pm

    # Retorna a derivada
    return d
end

# ω é um valor de referencia 
function norma_dω(ωn::Vector,U0::Matrix,malha::LFrame.Malha,x::Vector,fdkparam::Function,fdmparam::Function,fmparam::Function, P::Float64)

    # frequência de referência
    ω = 1.0 #mean(ωn)
    
    # Inicializa o somatório
    sum1 = 0.0

    # loop pelas frequências
    for i in eachindex(ωn)

        # ∑ (ωi/ω)^(-p)
        sum1 += (ωn[i] / ω)^(-P)

    end

    # fator global
    S = sum1^((-1/P) - 1)

    # vetor derivada
    D = zeros(length(x))

    for i in eachindex(ωn)

        # dωi/dx
        dωi = dω(ωn[i],U0[:, i],malha,x,fdkparam,fdmparam,fmparam)

        # derivada normalizada
        D .+= S * (ωn[i] / ω)^(-P - 1) * (dωi ./ ω)

    end

    return D
end


# Rotina que volta uma matriz de derivadas da tensão equivalente em relação as variáveis de projeto
function norma_dσ(Λ_tio::Vector{Vector{Float64}},σe_tio::Vector{Matrix{Float64}},S::Vector{Matrix{Float64}},Pi::Vector{Vector{Matrix{Float64}}},malha::LFrame.Malha,U::Vector{Float64},x::Vector,fkparam::Function,fdkparam::Function,fσparam::Function,fdσparam::Function,P::Float64,s::Float64,σesc::Float64)


    # Dados da estrutura de malha
    dados_ele = malha.dados_elementos
    dicionario_mat = malha.dicionario_materiais
    dicionario_geo = malha.dicionario_geometrias
    conect = malha.conect
    coordenadas    = malha.coord
    ndof           = 6 * malha.nnos
    ne = malha.ne
    L = malha.L
    nos = malha.nnos
    apoios = malha.apoios

    # Descobre quais gdls são livres
    dofs_l = LFrame.dofs_livres(nos,apoios)
    
    # Monta a matriz K global e pega o numero de graus de liberdade
    KG = LFrame.Monta_Kg(malha,x,fkparam)

    Kg = KG[dofs_l, dofs_l]

    V = [1 0; 0 3]

    #  Inicializa T0 e Q
    Q = zeros(ndof, 2*ne)
    T0 = zeros(2*ne)

    # Inicializa o problema adjunto
    Γ = zeros(ndof, 2*ne)
    
    Λ  = similar(Λ_tio)
    σe = similar(σe_tio)
    # loop pelos elementos
    for ele in 1:ne

        # dados do elemento  nao estao no central principal
        Iz, Iy, J0, A, α, E, G, geo, ρ = LFrame.Dados_fundamentais(ele, dados_ele, dicionario_mat, dicionario_geo)

        # Descobre o gdls do elemento
        gls = LFrame.Gls(ele,conect)

        # transformação de coordenadas
        T = LFrame.Rotacao3d(ele,conect,coordenadas,α)

        # Rigidez do elemento 
        Ke = LFrame.Ke_portico3d(E,Iz,Iy,G,J0,L[ele],A)
        
        ## loop pelos nos
        for no in 1:2

            ## indice do vetor de tensões equivalente
            idx = 2*(ele-1) + no

            # vetor local TEMPORÁRIO
            T1 = zeros(ndof)

            # parametrizada
            Λ[idx] = fσparam(x[ele]) .* Λ_tio[idx]

            ## termo T0 de cada indice
            T0[idx] = (sum(Λ[idx].^P))^((1/P) - 1)


            σe[idx] = fσparam(x[ele]).*σe_tio[idx]
            
            # olhar se vai denovo o sigmaparam
            for ino in eachindex(Pi[idx])
                T1[gls] .+= (Λ[idx][ino]^(P-2)) * vec(σe[idx][ino,:]' * (V * fσparam(x[ele]) * Pi[idx][ino] * S[idx] *fkparam(x[ele]) * Ke * T))
            end


            # monta coluna da matriz adjunta
            Q[:,idx] .= (s/σesc) .* T0[idx] .* T1

        end

    end

    # resolvendo o problema adjunto 
    Γ[dofs_l,:] = -(Kg \ Q[dofs_l,:])

    ## derivada da tensão equivalente em relação a cada elemento
    dσ = zeros(2*ne, ne)

    # Loop pelos elementos - calcula o termo direto (dKe/dxm é zero para m ≠ ele no SIMP)
    for ele in 1:ne

        # dados do elemento não estão no central principal
        Iz, Iy, J0, A, α, E, G, geo, ρ = LFrame.Dados_fundamentais(ele, dados_ele, dicionario_mat, dicionario_geo)

        # Descobre os gdls do elemento
        gls = LFrame.Gls(ele, conect)

        # transformação de coordenadas
        T = LFrame.Rotacao3d(ele, conect, coordenadas, α)

        # Rigidez do elemento no sistema local 
        Ke = LFrame.Ke_portico3d(E, Iz, Iy, G, J0, L[ele], A)

        # derivada da rigidez no sistema local 
        dKe = Ke * fdkparam(x[ele])

        # Converte o deslocamento para o sistema local 
        Ue = T * U[gls]

        ## loop pelos nós do elemento de pórtico
        for no in 1:2

            # indice do vetor de tensões equivalente
            idx = 2*(ele-1) + no

            # parametrizada
            Λ[idx] = fσparam(x[ele]) .* Λ_tio[idx]
            σe[idx] = fσparam(x[ele]).*σe_tio[idx]


            # Inicializa o termo direto para esse indice (Temporario)
            termo_direto = 0.0

            # Loop pelos nós da seção transversal - termo direto
            # dKe/dxm é não nulo apenas para m == ele (SIMP)
            for ino in eachindex(Pi[idx]) 

                # derivada da param * tensao sem param
                termo1 = fdσparam(x[ele]) .* σe_tio[idx][ino,:]

                # param * derivada da tensao
                termo2 = fσparam(x[ele]) .* (Pi[idx][ino] * S[idx] * dKe * Ue)

                # termo todo
                termo_direto += (s / σesc) * T0[idx] * ((Λ[idx][ino])^(P - 2)) * ((σe[idx][ino,:])' * V * (termo1 .+ termo2))[1]
            end
            dσ[idx, ele] += termo_direto

        end
    end

    # termo indireto: γe_m'*dKm*Um para TODOS os pares (idx, m)
    # separado para evitar recalcular dados de m dentro do loop de ele
    for m in 1:ne

        # dados do elemento m
        Iz_m, Iy_m, J0_m, A_m, α_m, E_m, G_m, geo_m, ρ_m = LFrame.Dados_fundamentais(m, dados_ele, dicionario_mat, dicionario_geo)

        # gdls do elemento m
        gls_m = LFrame.Gls(m, conect)

        # transformação de coordenadas do elemento m
        T_m = LFrame.Rotacao3d(m, conect, coordenadas, α_m)

        # Rigidez do elemento m no sistema local
        Ke_m = LFrame.Ke_portico3d(E_m, Iz_m, Iy_m, G_m, J0_m, L[m], A_m)

        # derivada da rigidez do elemento m no sistema local
        dKe_m = Ke_m * fdkparam(x[m])

        # deslocamento do elemento m no sistema local
        Ue_m = T_m * U[gls_m]

        # termo indireto para TODOS os idx - vetorizado
        dKUe_m = dKe_m * Ue_m                  
        Γe_m   = T_m * Γ[gls_m, :]         
        

        # loop por todas as restrições idx
        for idx in 1:(2*ne)
            dσ[idx, m] += (Γe_m[:, idx]' * dKUe_m)[1]
        end
    end

    # retorna a matriz de derivadas da tensão equivalente em relação as variáveis de projeto
    return  dσ

end
