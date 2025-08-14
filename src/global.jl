#
# Monta a matriz global de rigidez
#
function RigidezGlobal(nn,ne,XY,IJ)

    # Inicializa a matriz de rigidez
    K = spzeros(nn,nn)

    # Loop pelos elementos finitos
    for e = 1:ne
        
        # Coordenadas dos nós do elemento 
        X,Y = MontaXY(e,IJ,XY) 

        # Matriz de rigidez local do elemento
        Ke = MontaKe(X,Y)

        # Determina os gls globais do elemento
        # que neste caso são os próprios nós
        ge = IJ[e,:]

        # Sobreposição da matriz de rigidez
        for i=1:4
            for j=1:4
                K[ge[i],ge[j]] = K[ge[i],ge[j]] + Ke[i,j]
            end 
        end 

    end #e

    # Retorna a matriz de rigidez global 
    return K

end


#
# Acumula os valores de Inércia de cada elemento da malha
#
function Propriedades_secao(Φ, ne, IJ, XY)

    # Inicializa os somatórios globais 
    Jeq = 0.0
    Qx  = 0.0
    Qy  = 0.0
    Ix  = 0.0
    Iy  = 0.0
    Ixy = 0.0
    
    # Vamos devolver um vetor com as áreas de cada elemento
    A = zeros(ne)

    # Loop pelos elementos 
    for ele = 1:ne

        # Recupera os nós do elemento
        nos = IJ[ele,:]

        # Recupera as coordenadas do elemento 
        X,Y = MontaXY(ele,IJ,XY) 

        # Recupera os valores de Φ nos nós do elemento 
        ϕ_e = Φ[nos]

        # Contribuição do elemento 
        areae, Qxe, Qye, Ixe, Iye, Ixye, Je = Inercias_elemento(ϕ_e,X,Y)

        # Soma a contribuição do elemento para as inércias globais
        Qx  += Qxe
        Qy  += Qye
        Ix  += Ixe
        Iy  += Iye
        Ixy += Ixye
        Jeq += Je
        
        # Armazena a área de cada elemento
        A[ele] = areae

    end

    # Agora podemos retornar os valores, lembrando que são 
    # em relação à origem do sistema de referência
    return A, Qx, Qy, Ix, Iy, Ixy, Jeq

end