# rotina para o calculo da tensao equivalente 
function tensoes(arquivoEsf,ne,P)

    tensao = zeros(2 * ne)

    contador = 1
    # loops pelos elementos
    for ele in 1:ne

        ## tensao do elemento no nó 1
        tensao[contador] = Pos_processamento(arquivoEsf,ele,1,P,true)
        tensao[contador+1] = Pos_processamento(arquivoEsf,ele,2,P,true)
        
        contador += 2
    end
    return tensao
end
