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
function dω(ωn::Vector,U0::Matrix,malha::LFrame.Malha,x::Vector,fdkparam::Function,gdmparam::Function)    

    # numero de frequencia 
    nω = length(ωn)

    # Iniciando a derivada
    D = zeros(nω,malha.ne)

    # Alias dos dados da estrutura de malha
    dados_ele = malha.dados_elementos
    dicionario_mat = malha.dicionario_materiais
    dicionario_geo = malha.dicionario_geometrias
    conect = malha.conect
    coordenadas = malha.coord
    L = malha.L

    for i in eachindex(ωn)
        
        # Recupera todas as frequencia
        ω = ωn[i]

        ϕ = U0[:,i]

        # Loop pelos elementos
        for ele in 1:malha.ne
        
            # Recupera dados do elemento 
            Iz, Iy, J0, A, α, E, G, geo,ρ = LFrame.Dados_fundamentais(ele, dados_ele, dicionario_mat, dicionario_geo)

            # modos no global
            u0g = ϕ[nω] 
       
            # matriz de tranformação de coord(ver se é só rotacao ou translação)
            T = LFrame.Rotacao3d(ele,conect,coordenadas,α)

            # Para o local
            u0l = T'*u0g
        
            # Rigidez local do elemento 
            Ke = LFrame.Ke_portico3d(E,Iz,Iy,G,J0,L[ele],A)

            # derivada da rigidez em relação a variavel
            dK = T'*Ke*T*fdkparam(x[ele])

            # Massa local do elemento
            Me = LFrame.Me_portico3d(A,ρ,L[ele])

            # derivada da massa em relação a variavel
            dM = T'*Me*T*gdmparam(x[ele])

            # variavel auxialiar para facilitar
            aux = dK - (ω^2)*dM

            # denominador
            pm = dot(u0l,Me,u0l)

            # Deriada da frequencia ωn em relação ao elemento
            D[i,ele]  =  (1/(2*ω)*pm)*dot(u0l,aux,u0l)

        end
    end
    return D
end
