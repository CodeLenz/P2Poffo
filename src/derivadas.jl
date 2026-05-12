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
function dω(ωn::Float64,U0::Vector,malha::LFrame.Malha,x::Vector,fdkparam::Function,fdmparam::Function)    

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

    # Iniciando o denominador
    pm = 0.0

    # Inicializa o vetor com os modos globais
    ϕg = zeros(nos*6)

    # Descobre quais gdls são livres
    dofs_l = LFrame.dofs_livres(nos,apoios)

    # Associa o gdl livre ao vetor de modos
    ϕg[dofs_l] .= U0

    # Loop pelos elementos
    for ele in 1:malha.ne
        
        # Recupera dados do elemento 
        Iz, Iy, J0, A, α, E, G, geo, ρ = LFrame.Dados_fundamentais(ele, dados_ele, dicionario_mat, dicionario_geo)

        # Rigidez local do elemento 
        Ke = LFrame.Ke_portico3d(E,Iz,Iy,G,J0,L[ele],A)

        # Massa local do elemento
        Me = LFrame.Me_portico3d(A,ρ,L[ele])

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
        ϕ0l = T'*ϕ0g
    
        # variavel auxialiar para facilitar
        aux = dK - (ωn^2)*dM

        # Denominador
        pm += dot(ϕ0l,Me,ϕ0l)

        # Derivada do elemento
        D[ele] = dot(ϕ0l,aux,ϕ0l)/(2*ωn)
    end

    # Deriadas da frequencia ωn em relação a todos os elemento
    d = D / pm

    # Retorna a derivada
    return d
end

# ω é um valor de referencia 
function norma_dω(ωn::Vector,U0::Matrix,malha::LFrame.Malha,x::Vector,fdkparam::Function,fdmparam::Function, P::Float64)

    # frequência de referência
    ω = mean(ωn)
    
    # Inicializa o somatorio
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
        dωi = dω(ωn[i],U0[:, i],malha,x,fdkparam,fdmparam)

        # derivada normalizada
        D += S * (ωn[i] / ω)^(-P - 1) * (dωi / ω)
    end

    return D
end


### Tensoes
function norma_dσ(σeq::Vector,σe::Vector{Vector{Float64}},S::Vector{Matrix{Float64}},Pi::Vector{Matrix{Float64}},malha::LFrame.Malha,U::Vector,x::Vector,fkparam::Function,fdkparam::Function,P::Float64)

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
    
    ## termo constante
    T0 = (sum(σeq.^P))^((1/P) - 1)

    # Monta a matriz K global e pega o numero de graus de liberdade
    KG = LFrame.Monta_Kg(malha,x,fkparam)

    Kg = KG[dofs_l, dofs_l]

    # inicializando T1
    T1 = zeros(ndof)

    V = [1 0;
         0 3]

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
            
            # Parte T1 para o problema adjunto
            T1[gls] .+= (σeq[idx]^(P - 2)) *vec(σe[idx]' * (V * Pi[idx] * S[idx] * Ke * T))
        end
    end

    # variavel auxiliar para a solução do problema adjunto
    aux = (T0 * T1)[dofs_l]

    # resolvendo o problema adjunto
    λ = Kg\aux

    # vetor de multiplicadores de lagrange no global
    λg = zeros(ndof)
    # associando os multiplicadores de lagrange aos gdls livres
    λg[dofs_l] = λ

    ## derivada da tensão equivalente em relação a cada elemento
    dσ = zeros(ne)

    
    for ele in 1:ne

        # dados do elemento  nao estao no central principal
        Iz, Iy, J0, A, α, E, G, geo, ρ = LFrame.Dados_fundamentais(ele, dados_ele, dicionario_mat, dicionario_geo)

        # Descobre o gdls do elemento
        gls = LFrame.Gls(ele,conect)

        # transformação de coordenadas
        T = LFrame.Rotacao3d(ele,conect,coordenadas,α)

        # Rigidez do elemento 
        Ke = LFrame.Ke_portico3d(E,Iz,Iy,G,J0,L[ele],A)

        # derivada da rigidez 
        dKe = Ke*fdkparam(x[ele])

        Ue = T' * U[gls]

        ## loop pelos nos
        for no in 1:2

            ## indice do vetor de tensões equivalente
            idx = 2*(ele-1) + no
            
            # termo direto da derivada da tensão equivalente em relação a cada elemento
            termo_direto = (σeq[idx]^(P - 2)) * (σe[idx]' * V * Pi[idx] * S[idx] * dKe * Ue)[1]

            # termo indireto da derivada da tensão equivalente em relação a cada elemento
            termo_indireto = (λg[gls]' * T' * dKe * U[gls])[1]

            dσ[ele] += termo_direto + termo_indireto
        end
    end

    return  dσ
end

