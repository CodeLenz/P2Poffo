# rotina para o calculo da tensao equivalente 
function tensoes(arquivoEsf,ne,P,iter,posfile)

    tensao = zeros(2 * ne)
    SS = Vector{Matrix{Float64}}(undef, 2 * ne)
    Pi = Vector{Matrix{Float64}}(undef, 2 * ne)

    contador = 1
    # loops pelos elementos
    for ele in 1:ne

        ## tensao do elemento no nó 1
        tensao[contador],SS[contador],Pi[contador] = Pos_processamento(arquivoEsf,ele,1,P,iter,posfile)
        tensao[contador+1],SS[contador+1],Pi[contador+1] = Pos_processamento(arquivoEsf,ele,2,P,iter,posfile)
        
        contador += 2
    end
    # Tensao_ele1_no1,Tensao_ele1_no2...
    return tensao,SS,Pi
end
