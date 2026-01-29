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

function dofs_livres(nos,apoios)

    # numero total de GDLs
    ngdl = nos * 6
    
    # Vetor de GDLs fixos
    fixos = Int[]

    # Loop pelos apoios
    for i in axes(apoios,1) 

        # Recupera o nó
        no    = Int(apoios[i,1])  
        
        # Recupera o GDL
        gdl   = Int(apoios[i,2])   

        # Recupera o valor
        valor = apoios[i,3]        

        # Aloca o gdl local no global 
        glg = 6*(no-1) + gdl

        # Verifica se o deslocamento prescrito é zero 
        if valor == 0
            push!(fixos, glg)
            #print("O nó $no está com o GDL $gdl fixo \n")
        else
            error("O nó $no está com o GDL $gdl livre")
        end
    end

    # GDLs livres
    livres = Int[]
    
    # Loop por todos os GDLs globais
    for i in 1:ngdl

        # Verifica se i não é um gdl fixo
        if !(i in fixos)

            # Armazena no vetor livres
            push!(livres, i)
        end
    end
    
    return livres
end