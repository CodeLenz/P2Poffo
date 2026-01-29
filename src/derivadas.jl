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
    dofs_l = dofs_livres(nos,apoios)

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
