
#
# Monta a matriz N em um ponto (r,s)
#
function MontaN(r,s)

   # Calcula as funções no ponto 
   N1 = (1/4)*(1-r)*(1-s)
   N2 = (1/4)*(1+r)*(1-s)
   N3 = (1/4)*(1+r)*(1+s)
   N4 = (1/4)*(1-r)*(1+s)

   # Devolve a matriz
   N = [N1  N2  N3  N4 ]

   return N
end

#
# Monta o vetor local de forças de corpo devido a ao 
# 2
#
function MontaFBe(X,Y)

    # Aloca o vetor
    Fb = zeros(4)

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
  
       # Monta a matriz N no ponto
       N = MontaN(r,s)

       # Acumula o vetor 
       Fb .= Fb + transpose(N)*2*dJ

    end #i

    # Retorna o vetor
    return Fb
    
end


#
# Monta o vetor global de forças devido as forças de corpo
#
function ForcabGlobal(nn,ne,XY,IJ)

    # Inicializa o vetor global 
    F = zeros(nn)

    # Laço pelas linhas de FB
    for ele=1:ne

        # Coordenadas dos nós do elemento 
        X,Y = MontaXY(ele,IJ,XY) 

        # Calcula a integral no domínio do elemento, retornando 
        # um vetor 8 × 1
        Fb = MontaFBe(X,Y)

        # Determina os graus de liberdade globais do elemento 
        ge = IJ[ele,:]

        # Posiciona o vetor local 8 × 1 no vetor global 
        for j=1:4
            F[ge[j]] = F[ge[j]] + Fb[j]
        end
        
    end #i

    # Retorna o vetor global de forças de corpo
    return F

end

