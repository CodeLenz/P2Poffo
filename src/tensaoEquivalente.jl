# rotina para o calculo da tensao equivalente 
function tensoes(arquivoEsf,malha,iter,posfile,cache_secoes = Dict())

    # numero de elementos 
    ne = malha.ne 

    # tensao eqv maxima de cada no
    Λ = Vector{Vector{Float64}}(undef, 2 * ne)

    # matriz S associada a cada no
    S = Vector{Matrix{Float64}}(undef, 2 * ne)

    # matriz Pi associada a cada no
    Pi = Vector{Vector{Matrix{Float64}}}(undef, 2 * ne)

    # vetor da tensao local critica [σxx σxy]
    σe = Vector{Matrix{Float64}}(undef, 2 * ne) 

    linhas = readlines(arquivoEsf)
    path_base = dirname(arquivoEsf)

    contador = 1
    # loops pelos elementos
    for ele in 1:ne

        ## tensao do elemento no nó 1 e 2
        Λ[contador],S[contador],Pi[contador],σe[contador] = tensao_vonMises(linhas,path_base,ele,1,iter,cache_secoes,posfile)
        Λ[contador+1],S[contador+1],Pi[contador+1],σe[contador+1] = tensao_vonMises(linhas,path_base,ele,2,iter,cache_secoes,posfile)

        
        contador += 2
    end
    # Tensao_ele1_no1,Tensao_ele1_no2...
    return Λ,S,Pi,σe

end


#
# Calcula as tensões em todos os pontos (nós) de um elemento <ele> e nó <1/2>
# devolvendo a tesão equivalente máxima na seção transversal, a matriz S e a matriz Pi associada a esse ponto
#
# arquivo_esforcos é gerado pelo LFrame
#
function tensao_vonMises(linhas, path_base, ele, no, iter, cache_secoes, posfile=false, ϵ=1E-6)

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

        zl_yl, grad_xy, grad_n, Pi_geo = precomputa_geometria(
            (nn=nn, XY=XY, α=α, centroide=centroide,
            ∇Φ=∇Φ, vizinhos=vizinhos,
            area=area, Izl=Izl, Iyl=Iyl, Jeq=Jeq)
        )

        dados_secao = (
            centroide=centroide, area=area, Izl=Izl, Iyl=Iyl,
            Jeq=Jeq, α=α, ∇Φ=∇Φ, nn=nn, XY=XY, ne2d=ne2d,
            IJ=IJ, etypes=etypes, vizinhos=vizinhos,
            # geometria pré-computada:
            zl_yl=zl_yl, grad_xy=grad_xy,
            grad_n=grad_n, Pi_geo=Pi_geo
        )
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

        S = [0 0 0 0 0 0 1 0 0 0 0 0;
             0 0 0 0 0 0 0 0 0 1 0 0;
             0 0 0 0 0 0 0 0 0 0 0 1;
             0 0 0 0 0 0 0 0 0 0 1 0]
    end

    # A tensão normal devido a teoria de barra é cte
    σN = N/area

    # Podemos calcular as constantes de proporcionalidade
    cteMy =  My/Iyl
    cteMz = -Mz/Izl

    # Extrai do cache
    zl_yl   = dados_secao.zl_yl
    grad_xy = dados_secao.grad_xy
    grad_n  = dados_secao.grad_n
    Pi_geo  = dados_secao.Pi_geo
    nn      = dados_secao.nn

    αμ = T / Jeq

    # Monta Pi com o sign(T) e αμ já aplicados — vetorizado, sem loop
    Pi = Vector{Matrix{Float64}}(undef, nn)
    for i in 1:nn
        p = copy(Pi_geo[i])
        p[2, 2] = sign(T) * grad_n[i] / Jeq
        Pi[i] = p
    end

    # Tensões — totalmente vetorizadas
    σxx = σN .+ cteMy .* zl_yl[:, 1] .+ cteMz .* zl_yl[:, 2]
    σxy = αμ .* sqrt.(grad_xy[:, 1].^2 .+ grad_xy[:, 2].^2)  # = αμ * grad_n

    σe = [σxx σxy]

    
    # matriz V de Von Mises para calcular a tensão equivalente
    V = [1 0;0 3]

    Λ = sqrt.(σe[:, 1].^2 .* V[1,1] .+ σe[:, 2].^2 .* V[2,2] .+ ϵ^2)

   if posfile
        σ = zeros(nn, 5)
        σ[:, 1] .= σN
        σ[:, 2] .= cteMy .* zl_yl[:, 1]
        σ[:, 3] .= cteMz .* zl_yl[:, 2]
        σ[:, 4] .= αμ .* grad_xy[:, 1]
        σ[:, 5] .= αμ .* grad_xy[:, 2]

        pos_file_node = joinpath(path_base, nome_secao*"_iter$(iter)_ele$(ele)_No$(no).pos")
        Lgmsh_export_init(pos_file_node, nn, ne2d, XY, etypes, IJ)
        Lgmsh_export_nodal_scalar(pos_file_node, σ[:,1], "σxxN")
        Lgmsh_export_nodal_scalar(pos_file_node, σ[:,2], "σxxMY")
        Lgmsh_export_nodal_scalar(pos_file_node, σ[:,3], "σxxMz")
        Lgmsh_export_nodal_scalar(pos_file_node, σ[:,4], "σzyT")
        Lgmsh_export_nodal_scalar(pos_file_node, σ[:,5], "σzxT")
        Lgmsh_export_nodal_scalar(pos_file_node, Λ,    "σvon-Mises")
    end

    # retorna a tensao equivalente de todos os nos da seção, S do nó do portico, 
    # Pi de todos os nós da seção e a matriz com os estados de tensao de cada nó da seção (σxx σxy)
    return Λ,S,Pi,σe
end
