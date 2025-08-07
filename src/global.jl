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
