#
# Rotina que mudas as coordenadas para o centroide e rotaciona para o estado principal
#

function Muda_coordenada(x,y,α,xc,yc)

    # Inicia o vetor coordenadas finais 2x1
    coordl = zeros(2)

    # Translação da coordenada para a referência do centroide
    xl = xc - x
    yl = y - yc

    # Armazena em um vetor 2x1
    coord = [xl;yl]

    # Recupera a matriz de rotação 4x4
    R = Matriz_rotacao(α)

    # Rotaciona as coordenadas 
    coordl = R * coord

    return coordl
end

function Matriz_rotacao(α)

    # Calcula a matriz de rotação
    R = [cosd(α) -sind(α);
         sind(α)  cosd(α)]
    
    # Retorna a matriz     
    return R
end
