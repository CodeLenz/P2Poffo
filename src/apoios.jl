#
# Aplica as condições de contorno essenciais
# Homogêneas
#
function AplicaCCH(nn,na,AP,K,F)

    # Laço sobre as linhas de AP
    for l = 1:na

        # Nó
        no = Int(AP[l,1])

        # Loop para zerar as linhas e colunas
        for i=1:nn
            
            # Zera a linha glg
            K[no,i] = 0

            # Zera a coluna glg
            K[i,no] = 0

        end

        # Coloca 1 na diagonal da matriz
        K[no,no] = 1

        # Zera a linha em F
        F[no] = 0

    end #l

    return K,F

end

