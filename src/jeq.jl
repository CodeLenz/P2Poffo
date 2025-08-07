#
# Calcula o J_eq da seção, pelo somatório das integrais 
# de Φ em cada elemento finito da malha
#
function Jequivalente(Φ, ne, IJ, XY)

    # Inicializa o somatório para o Jeq
    jeq = 0.0

    # Inicializa o somatório para a área
    area = 0.0

    # Loop pelos elementos 
    for ele = 1:ne

        # Recupera os nós do elemento
        nos = IJ[ele,:]

        # Recupera as coordenadas do elemento 
        X,Y = MontaXY(ele,IJ,XY) 

        # Recupera os valores de Φ nos nós do elemento 
        ϕ_e = Φ[nos]

        # Contribuição do elemento 
        je,ae = Jelemento(ϕ_e,X,Y)

        # Soma a contribuição do elemento para o Jeq
        jeq = jeq + je

        # Soma a área
        area = area + ae
       
    end

    # J equivalente do elemento e área
    return 2*jeq, area

end