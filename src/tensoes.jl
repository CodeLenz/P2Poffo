#
# Calcula as tensões tangenciais no centro de cada 
# elemento da malha
#
function Tensoes(Φ,ne,IJ,XY)

    # Cada linha é um elemento. Coluna 1 é τ_zx e coluna 2 é τ_zy
    τ = zeros(ne,2)

    # Loop por cada elemento da malha
    for ele=1:ne
        
        # Recupera os nós do elemento
        nos = IJ[ele,:]

        # Recupera as coordenadas do elemento 
        X,Y = MontaXY(ele,IJ,XY) 

        # Recupera os valores de Φ nos nós do elemento 
        ϕ_e = Φ[nos]

        # Monta a matriz B do elemento no ponto central 
        # r=s=0
        # Derivadas das funções de interpolação 
        # em relação a rs neste ponto
        dNrs = MontadN(0,0)

        # Calcula a matriz J no ponto
        J = MontaJ(dNrs,X,Y)

        # Mapeia as derivadas para xy
        dNxy = CorrigedN(dNrs,J)
  
        # Monta a matriz B no ponto
        B = MontaB(dNxy)

        # Tensões no centro do elemento 
        τ[ele,:] .= B*ϕ_e

    end

    # Retorna as tensões 
    return τ
    
end