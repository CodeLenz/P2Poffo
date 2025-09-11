#
# Calcula as tensões em todos os pontos (nós) de um elemento <ele> e nó <1/2>
# devolvendo uma matriz com nnos × 5 valores de tensão (σN σMy σMz σxy σxz)
#
# arquivo_esforcos é gerado pelo LFrame
#
function Pos_processamento(arquivo_esforcos, ele, no)

    # Testa se nó é válido
    no in [1;2] || error("Pos_processamento:: nó inválido $no")
    
    # Testa se elemento é maior ou igual a 1
    ele >=1 || error("Pos_processamento:: elemento deve ser >=1")

    # Abre o arquivo de esforços e lê a linha relativa ao elemento 
    fd = open(arquivo_esforcos,"r")

    # Verifica se o arquivo é válido
    linha = readline(fd)
    occursin("Esforcos",linha) || error("Pos_processamento:: arquivo de esforços inválido")

    # Le a linha associada ao elemento ele
    contador_linha = 1
    while !eof(fd) && contador_linha<=ele

        # Lê a linha 
        linha = readline(fd)

        # Incrementa o contador 
        contador_linha += 1

    end

    # Se contador_linha-1 for menor do que ele, então não temos o elemento no arquivo
    if contador_linha-1 < ele 
       error("Pos_processamento:: elemento $ele não existe no arquivo de esforços")
    end

    # Separa a linha por tokens
    dados = split(linha)

    # Testa para ver se temos 13 informações na linha
    length(dados)==13 || error("Pos_processamento:: Dados para o elemento $ele não tem a dimensão correta")

    # O nome da seção transversal é a primeira informação da linha 
    nome_secao = dados[1]

    # O nome da seção será 
    arquivo_secao = nome_secao*".geo"

    # O nome do arquivo da malha será
    arquivo_malha = nome_secao*".msh"

    # Roda o pré-processamento e obtem todos os dados da seção transversal
    centroide, area, Izl, Iyl, Jeq, α, ∇Φ = Pre_processamento(arquivo_secao, false)

    # Dependendo do nó, pegamos os esforços internos
    if no==1
       
        N  = -1*parse(Float64,dados[2])
        My = -1*parse(Float64,dados[6])
        Mz = -1*parse(Float64,dados[7])
        T  = -1*parse(Float64,dados[5])

    else

        N  = parse(Float64,dados[8])
        My = parse(Float64,dados[12])
        Mz = parse(Float64,dados[13])
        T  = parse(Float64,dados[11])

    end

    # Vamos precisar dos dados da malha
    nn,XY,ne,IJ,MAT,na,AP,etypes,centroides = ConversorFEM(arquivo_malha)

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

    # Loop pelos nós 
    for no = 1:nn

        # Coordenadas do nó 
        x,y = XY[no,1:2]

        # Converte para (z',y')
        zl,yl = Muda_coordenada(x,y,α,centroide[1],centroide[2])
        
        # Agora podemos fazer as médias nodais SIMPLES das tensões 
        # tangenciais para cada nó 
        vizinhos =  Vizinhos_no(no,IJ)

        # Calcula a média nodal SIMPLES do gradiente
        ∇xΦ = mean(∇Φ[vizinhos,1])
        ∇yΦ = mean(∇Φ[vizinhos,2])

        # Guarda nas colunas 
        σ[no,:] = [σN cteMy*zl cteMz*yl αμ*∇xΦ αμ*∇yΦ]


    end
    # Retira os caminhos do nome do arquivo
    mshfile2 = basename(arquivo_malha)
      
    # Gera o nome do arquivo .pos 
    posfile = replace(mshfile2,".msh"=>".pos")

    # Inicializa o arquivo de saída
    Lgmsh_export_init(posfile,nn,ne,XY,etypes,IJ) 

    # Exporta o campo σ 
    ### OBS para a arrumar: aqui só está mostrando o sigma do Mz 
    Lgmsh_export_nodal_scalar(posfile, σ[:,3],"σ")
    
    # Retorna a matriz com as tensões
    return σ

end

#
# Dado um nó, devolve uma lista com todos os elementos que o contem 
#
function Vizinhos_no(no,IJ)

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