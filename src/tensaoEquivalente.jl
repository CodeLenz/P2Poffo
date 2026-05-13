# rotina para o calculo da tensao equivalente 
function tensoes(arquivoEsf,malha,P,iter,posfile)

    # numero de elementos 
    ne = malha.ne 

    # tensao eqv maxima de cada no 
    tensao = zeros(2 * ne)

    # matriz S associada a cada no
    SS = Vector{Matrix{Float64}}(undef, 2 * ne)

    # matriz Pi associada a cada no
    Pi = Vector{Matrix{Float64}}(undef, 2 * ne)
    
    # vetor da tensao local critica [σxx σxy]
    σe = Vector{Vector{Float64}}(undef, 2 * ne)


    # Cria um dicionário para armazenar os dados das seções transversais, evitando ler o mesmo arquivo várias vezes
    cache_secoes = Dict()
    linhas = readlines(arquivoEsf)
    path_base = dirname(arquivoEsf)

    contador = 1
    # loops pelos elementos
    for ele in 1:ne

        ## tensao do elemento no nó 1
        tensao[contador],SS[contador],Pi[contador],σe[contador] = tensao_vonMises(linhas,path_base,ele,1,P,iter,cache_secoes,posfile)
        tensao[contador+1],SS[contador+1],Pi[contador+1],σe[contador+1] = tensao_vonMises(linhas,path_base,ele,2,P,iter,cache_secoes,posfile)

        
        contador += 2
    end
    # Tensao_ele1_no1,Tensao_ele1_no2...
    return tensao,SS,Pi,σe
end


#
# Calcula as tensões em todos os pontos (nós) de um elemento <ele> e nó <1/2>
# devolvendo a tesão equivalente máxima na seção transversal, a matriz S e a matriz Pi associada a esse ponto
#
# arquivo_esforcos é gerado pelo LFrame
#
function tensao_vonMises(linhas, path_base, ele, no, P, iter, cache_secoes, posfile=false)

    # Testa se nó é válido
    no in [1;2] || error("Pos_processamento:: nó inválido $no")
    
    # Testa se elemento é maior ou igual a 1
    ele >=1 || error("Pos_processamento:: elemento deve ser >=1")

    # verifica cabecalho
    occursin("Esforcos", linhas[1]) ||error("Pos_processamento:: arquivo de esforços inválido")

    # pega linha do elemento (primeira linha é o cabecalho, entao elemento 1 começa na linha 2)
    linha = linhas[ele + 1]

    # Separa a linha por tokens
    dados = split(linha)

    # Testa para ver se temos 13 informações na linha
    length(dados)==13 || error("Pos_processamento:: Dados para o elemento $ele não tem a dimensão correta")


    # Agora que lemos os esforços do elemento, podemos processar a malha de elementos finitos 
    # associada a essa seção transversal

    # O nome da seção transversal é a primeira informação da linha 
    nome_secao = dados[1]

    # O nome da seção será
    arquivo_secao = joinpath(path_base,nome_secao*".geo")

    # O nome do arquivo da malha será
    arquivo_malha = joinpath(path_base,nome_secao*".msh")

    # Verifica se os dados da seção já estão no cache
    if haskey(cache_secoes, nome_secao)

        # pega do dicionário se já tiver sido lido antes
        dados_secao = cache_secoes[nome_secao]

    else

        # utiliza a função de pré-processamento para ler os dados da seção transversal
        centroide, area, Izl, Iyl, Jeq, α, ∇Φ = Pre_processamento(arquivo_secao, false)

        # Vamos precisar dos dados da malha
        nn,XY,ne2d,IJ,MAT,na,AP,etypes,centroides = ConversorFEM(arquivo_malha)

        # Cria um vetor de vetores com os elementos que 
        # são vizinhos de cada nó da malha
        vizinhos = Vizinhos_no(nn,IJ)

        # armazena os dados da seção no dicionário para evitar ler o mesmo arquivo novamente
        dados_secao = (centroide = centroide,area = area,Izl = Izl,Iyl = Iyl,Jeq = Jeq,α = α,∇Φ = ∇Φ,nn = nn,XY = XY,ne2d = ne2d,IJ = IJ,etypes = etypes,vizinhos = vizinhos)

        # Armazena os dados da seção no cache
        cache_secoes[nome_secao] = dados_secao

    end

    # Dados da seção transversal do dict
    centroide = dados_secao.centroide
    area      = dados_secao.area
    Izl       = dados_secao.Izl
    Iyl       = dados_secao.Iyl
    Jeq       = dados_secao.Jeq
    α         = dados_secao.α
    ∇Φ        = dados_secao.∇Φ

    nn        = dados_secao.nn
    XY        = dados_secao.XY
    IJ        = dados_secao.IJ
    etypes    = dados_secao.etypes
    vizinhos  = dados_secao.vizinhos
    ne2d      = dados_secao.ne2d

    # Dependendo do nó, pegamos os esforços internos e a matriz S de acordo com a teoria
    if no==1
        
        N  = -parse(Float64, dados[2])
        T  = -parse(Float64, dados[5])
        My = -parse(Float64, dados[6])
        Mz = -parse(Float64, dados[7])


        S = -1*[1 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0 0 0 0 0;
            0 0 0 0 0 1 0 0 0 0 0 0;
            0 0 0 0 1 0 0 0 0 0 0 0]

    else

        N  = parse(Float64, dados[8])
        T  = parse(Float64, dados[11])
        My = parse(Float64, dados[12])
        Mz = parse(Float64, dados[13])

        S = [0 0 0 0 0 1 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 0 0 0 0 0 0 1;
            0 0 0 0 0 0 0 0 0 0 1 0]
    end

   

    # Agora podemos alocar a matriz de saída
    # (σN σMy σMz σxy σxz)
    σ = zeros(nn,5)

    # A tensão normal devido a teoria de barra é cte
    σN = N/area

    # Podemos calcular as constantes de proporcionalidade
    cteMy =  My/Iyl
    cteMz = -Mz/Izl

    # Calcula a cte αμ
    αμ = T/Jeq

    ∇ = zeros(nn)
    Pi = Vector{Matrix{Float64}}(undef, nn)
    # Loop pelos nós da malha da seção
    for ino = 1:nn

        # Coordenadas do nó 
        x,y = XY[ino,1:2]

        # Converte para (z',y')
        zl,yl = Muda_coordenada(x,y,α,centroide[1],centroide[2])
        
        # Agora podemos fazer as médias nodais SIMPLES das tensões 
        # tangenciais para cada nó (tentar colocar fora do loop, gerar um dicionario para reaproveitar)
        vizi = vizinhos[ino]

        # Calcula a média nodal SIMPLES do gradiente
        ∇xΦ = mean(∇Φ[vizi,1])
        ∇yΦ = mean(∇Φ[vizi,2])

        ∇[ino] = sqrt(∇xΦ^2 + ∇yΦ^2)

        # Guarda nas colunas 
        σ[ino,:] = [σN cteMy*zl cteMz*yl αμ*∇xΦ αμ*∇yΦ]

        # Calcula a matriz Pi associada a esse nó 2d
        Pi[ino] = [1/area   0   -yl/Izl   zl/Iyl;
                    0   (1/Jeq)*∇[ino]   0   0
                ]
    end
    
    # matriz V para o produto quadratico 
    V = [1 0;
         0 3]

    #tensao em cada elemento
    σxx  = σ[:,1] + σ[:,2] + σ[:,3]
    σxy  =  sqrt.(σ[:,4].^2 + σ[:,5].^2)

    # matriz de tensao
    σe = [σxx σxy]

    nn = size(σ, 1)
    σeq = zeros(nn)

    # loop por todos os nos e calcula a tensao eqv
    for i in 1:nn
        σeq[i] = sqrt(σe[i,:]' * V * σe[i,:])
    end

    # tira o maximo
    σeq_max = norm(σeq,P)
    @show σeq_max 

    # Recupera a matriz Pi associada a pior tensao, pois esse P * S * F = tensao equivalente máxima no 2d
    idx_max = argmax(σeq)
    Pi_critico = Pi[idx_max]
    
    # tensão equivalente crítica 
    σe_critico = [σxx[idx_max], σxy[idx_max]]

    if posfile 
        # Caminho para a pasta POS
        pos_file_node = joinpath(path_base,nome_secao*"_iter$(iter)_ele$(ele)_No$(no).pos")

        # Inicializa o arquivo de saída
        Lgmsh_export_init(pos_file_node,nn,ne2d,XY,etypes,IJ) 

        # Exporta o campo σ (falta "somar as tensões")
        Lgmsh_export_nodal_scalar(pos_file_node, σ[:,1],"σxxN")
        Lgmsh_export_nodal_scalar(pos_file_node, σ[:,2],"σxxMY")
        Lgmsh_export_nodal_scalar(pos_file_node, σ[:,3],"σxxMz")
        Lgmsh_export_nodal_scalar(pos_file_node, σ[:,4],"σzyT")
        Lgmsh_export_nodal_scalar(pos_file_node, σ[:,5],"σzxT")
        Lgmsh_export_nodal_scalar(pos_file_node, σeq,"σvon-Mises")

    end


    # Retorna o valor maximo na seção e elemento
    return σeq_max,S,Pi_critico,σe_critico
end


#
# Retorna um vetor de vetores com os elementos vizinhos a cada nó
#
function Vizinhos_no(nn,IJ)

    # Aloca um vetor de vetores 
    vizinhos = Vector{Vector{Int64}}(undef,nn)

    # Loop pelos nós da malha, chamando Vizinhos_no(no,IJ)
    for no=1:nn

        vizinhos[no] = _Vizinhos_no(no,IJ)

    end

    # Retorna os vizinhos 
    return vizinhos

end


#
# Dado um nó, devolve uma lista com todos os elementos que o contem 
#
function _Vizinhos_no(no,IJ)

    # Aloca um vetor para armazenar os vizinhos 
    vizinhos = Int64[]

    # Loop por todos os elementos da malha 
    ele = 1
    for nos in eachrow(IJ)

        # Ve se a conectividade do elemento contém o nó
        if no in nos
           push!(vizinhos,ele) 
        end

        # Incrementa o contador de elementos 
        ele += 1

    end

    # Retorna os vizinhos 
    return vizinhos

end