#
# Derivada das funções de interpolação 
# em relação a r e s
#
function MontadN(r,s)

    # Matriz com as derivadas das funções de interpo
    # primeira linha em relação a r
    # segunda linha em relação a s
    dNrs = (1/4)*[s-1   1-s   1+s -(1+s) ;
                  r-1 -(1+r)  1+r   1-r ]

    # Retorna a matriz
    return dNrs

end


#
# Monta a matriz J em um ponto (r,s)
#
function MontaJ(dNrs,X,Y)

    # Aloca a matriz J
    J = zeros(2,2)

    # Loop do somatório
    for i=1:4

        # dx/dr
        J[1,1] = J[1,1] + dNrs[1,i]*X[i]

        # dy/dr
        J[1,2] = J[1,2] + dNrs[1,i]*Y[i]

        # dx/ds
        J[2,1] = J[2,1] + dNrs[2,i]*X[i]

        # dy/ds
        J[2,2] = J[2,2] + dNrs[2,i]*Y[i]

    end #i

    # Retorna a matriz J no ponto (r,s)
    return J

end

#
# Corrige as derivadas dNrs para dNxy em um ponto 
# (r,s) do elemento
#
function CorrigedN(dNrs,J)

    # Inverte a matriz
    invJ = inv(J)

    # Calcula a correção 
    dNxy = invJ*dNrs

    # Retorna a correção 
    return dNxy

end

#
# Monta a matriz B no ponto (r,s)
#
function MontaB(dNxy)

    # Aloca a matriz B
    B = zeros(2,4)

    # Loop pelas funções de interpolação 
    for i=1:4
 
       # Derivada em relação a x
       B[1,i] = dNxy[1,i]

       # Derivada em relação a y
       B[2,i] = dNxy[2,i]
 
    end

    # Retorna a matriz B no ponto (r,s)
    return B

end



#
# Monta XY
#
function MontaXY(e,IJ,XY)

    # Aloca os vetores de saída
    X = zeros(4)
    Y = zeros(4)

    # Loop pelos nós deste elemento
    for n=1:4

        # Descobre quem é o nó i do elemento 
        no = IJ[e,n]

        # Recupera as coordeandas deste nó 
        # na matriz XY
        X[n] = XY[no,1]
        Y[n] = XY[no,2]

    end

    # Retorna os vetores 
    return X,Y

end


#
# Monta a matriz de rigidez local do elemento
#
function MontaKe(X,Y)

    # Aloca a matriz
    Ke = zeros(4,4)

    # Define os pontos de Gauss-Legendre
    pg = (1/sqrt(3))*[-1  1 1 -1 ;
                      -1 -1 1  1  ]

    # Laço pelos pontos de Gauss
    for i=1:4

       # Recupera as coordenadas do PG
       r = pg[1,i]
       s = pg[2,i]

       # Derivadas das funções de interpolação 
       # em relação a rs neste ponto
       dNrs = MontadN(r,s)

       # Calcula a matriz J no ponto
       J = MontaJ(dNrs,X,Y)

       # Calcula o determinante no ponto 
       dJ = det(J)

       # Mapeia as derivadas para xy
       dNxy = CorrigedN(dNrs,J)
  
       # Monta a matriz B no ponto
       B = MontaB(dNxy)

       # Acumula a matriz (somatório da integração)
       Ke = Ke + transpose(B)*B*dJ

    end #i

    # Retorna a matriz de rigidez do elemento
    return Ke
    
end


#
# Integra o campo Φ dentro do elemento a partir dos 
# valores nodais ϕ_e
#
function Jelemento(ϕ_e,X,Y)

    # Matriz dos pontos de Gauss (r, s)
    pg = (1/sqrt(3)) * [-1  1  1 -1;
                        -1 -1  1  1]

    # Inicializa o somátorio para o je
    je = 0.0  

    # Incializa o somatório para a área 
    area = 0.0

    # Loop pelos 4 pontos de Gauss
    for i=1:4

        # Recupera o ponto 
        r = pg[1, i]
        s = pg[2, i]
        
        # Funções de interpolação neste ponto
        N = MontaN(r, s) 
        
        # Derivadas das funções de interpolação 
        # em relação a rs neste ponto
        dNrs = MontadN(r,s)

        # Calcula a matriz J no ponto
        J = MontaJ(dNrs,X,Y)

        # Calcula o determinante no ponto 
        dJ = det(J)

        # Valor interpolado no ponto (r,s)
        ϕ = dot(N,ϕ_e)

        # Acumula o valor da integral do je
        je = je + ϕ*dJ

        # Acumula a área
        area = area + dJ

    end
    
    # Retorna o valor da integral de Φ em Ωe
    # e a área
    return je, area 

end
