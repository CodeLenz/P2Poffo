function tensoes(arquivo::AbstractString)

    # Recupera os dados da malha 
    U,malha = LFrame.Analise3D(arquivo, true; verbose=false)

    # Recupera o nome do arquivo que foi salvo os dados de esforços
    arquivoEsf = malha.nome_arquivo*".esf"
    
    tensao = zeros(2 * malha.ne)

    contador = 1
    # loop pelos elementos
    for ele in 1:malha.ne

        ## tensao do elemento no nó 1
        tensao[contador] = Pos_processamento(arquivoEsf,ele,1)
        tensao[contador+1] = Pos_processamento(arquivoEsf,ele,2)
        
        contador += 2
    end
    return tensao
end
