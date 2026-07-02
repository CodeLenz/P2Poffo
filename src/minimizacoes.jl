#
# Rotinas de integração
#
#
# V(x) <= V_sup = vf * V0
# 

# Parametrização SIMP 
Kparam(x,p=3) = x^p

# Derivada da parametrização 
dKparam(x,p=3) = p*x^(p-1)

# Parametrização SIMP 
Mparam(x,p=1) = x^p

# Derivada da parametrização 
dMparam(x,p=1) = p*x^(p-1)

σparam(x, p=3.0, q=2.0) = x^(p-q)

dσparam(x, p=3.0, q=2.0) = (p-q)*(x^(p-q-1))

#
# Rotina principal
#
function Main_Otim_OC(arquivo::AbstractString, fkparam::Function, fdkparam::Function, posfile=true; verbose=false,vf = 0.5, niter=3)

    # Chama o Analise3D com o nome do arquivo, para receber a estrutura de malha
    U0, malha = Analise3D(arquivo,posfile,verbose=verbose)

    # Numero de elementos 
    ne = malha.ne

    # Estimativa inicial das densidades relativas
    x0 = vf*ones(ne)

    # Calcula o volume de cada elemento sem considerar a 
    # parametrização 
    V = x0 .* Volumes(malha)

    # Volume total ta estrutura
    V0 = sum(V)

    # Volume limite
    V_sup = vf * V0
    
    # Derivada do volume é fixa 
    dV = V
    
    # Loop externo de otimização 
    for iter=1:niter

        # Calcula os deslocamentos
        U, _ = Analise3D(malha,posfile,x0=x0,kparam=[fkparam])
       
        # @show iter, sum(x0.*dV), U[8]

        # Deriva da compliance
        dC = dCompliance(malha,U,x0,fdkparam)    

        # Atualiza as densidades relativas utilizando o OC
        x0 .= OC(x0,dC,dV,V_sup,ne)
        
    end # loop externo
    
    # Adiciona uma vista com as densidades relativas ao arquivo de saída
    if posfile

        # Nome do arquivo com .pos 
        arquivo_pos =  replace(arquivo,".yaml"=>".pos")
        
        # Exporta a vista escalar 
        Lgmsh.Lgmsh_export_element_scalar(arquivo_pos,x0,"x")

    end

   # Retorna as densidades relativas 
   return x0

end


"""
    Main_Otim_Modal("arquivo.yaml",P2Poffo.Kparam,P2Poffo.dKparam,P2Poffo.Mparam,P2Poffo.dMparam,P2Poffo.fσparam,P2Poffo.fdσparam,n_modos=3,vf=0.5,P=8.0,niter=100)

Rotina para calcular a maior frequência natural de uma estrutura,
minimizando o volume de uma malha.

# Descrição
Esta função realiza um processo de otimização topológica com o objetivo
de maximizar a frequência natural da estrutura, sujeito a uma restrição
de volume. A malha e os parâmetros do problema são fornecidos por meio
de um arquivo `.yaml`.

# Argumentos
- `arquivo_yaml`: Arquivo contendo os dados da malha e propriedades do problema.

# Parâmetros Nomeados
- `vf`: Fração de volume utilizada para definir o volume máximo permitido.
- `P`: Parâmetro da norma `-P` usada na agregação das frequências.
       Valores maiores tornam a aproximação mais precisa, porém aumentam o custo computacional.
- `niter`: Número máximo de iterações caso o critério de convergência não seja atingido.

# Modelos de Material (SIMP)
- Rigidez:
    - Penalização com `p = 3`
    - Função: `P2Poffo.Kparam`
    - Derivada: `P2Poffo.dKparam`

- Massa:
    - Penalização com `p = 1`
    - Função: `P2Poffo.Mparam`
    - Derivada: `P2Poffo.dMparam`

# Observações
- As funções SIMP não estão exportadas diretamente, sendo necessário acessá-las
  via o módulo `P2Poffo`.
- O desempenho da otimização depende da escolha de `P`, podendo impactar
  significativamente o tempo de execução.

# Retorno
Retorna os parâmetros otimizados da malha e as frequências naturais associadas.
"""
function Main_Otim_Modal(arquivo::AbstractString,x=0.5, fkparam=P2Poffo.Kparam, fdkparam=P2Poffo.dKparam,fmparam=P2Poffo.Mparam, fdmparam=P2Poffo.dMparam,fσparam=P2Poffo.σparam, fdσparam=P2Poffo.dσparam, posfile=true; verbose=false,n_modos=6, vf = 0.5 ,P=8.0, niter=300,tol_f=1E-4,tol_g=1E-4,s = 2.0)

    # analise modal somente para a malha
    ωn,U0,malha = Modal3D(arquivo,posfile,n_modos=n_modos)

    # recupera a tensao de escoamento do material (assumindo que é a mesma para todos os elementos)
    σesc = malha.dicionario_materiais[malha.dados_elementos[1,1]]["S_esc"]

    # Numero de elementos 
    ne = malha.ne

    # Estimativa inicial das densidades relativas,
    x0 = x*ones(ne)

    # Calcula o volume de cada elemento sem considerar a 
    # parametrização 
    V = Volumes(malha)

    # Volume total ta estrutura
    V0 = sum(V)

    # Volume limite
    V_sup = vf * V0
    
    # Derivada do volume é fixa = volume inicial
    dV = V

    # Inicializa os vetores de densidades relativas das iterações
    x1 = copy(x0)
    x2 = copy(x1)
    xn = copy(x0)    

    # Limites máximos e mínimos para os deltas
    δ_max = 0.1
    δ_min = 0.01

    # inicializa o veltor de deltas
    δ = δ_min*ones(ne)

    # inicializa os Deltas
    Δ1 = x0 .- x1
    Δ2 = x1 .- x2

    # inicializa os valores inferiores e superiores de densidades relativas
    x_inf = 1E-3*ones(length(x0))
    x_sup = ones(length(x0))        

    # inicializa o vetor 
    ωx0 = ωn[1]
    ω1 = 0

    # nome do arquivo de esforços
    nomeEsf = malha.nome_arquivo

    # inicializa o vetor de derivadas da tensão equivalente em relação as variáveis de projeto
    dσ = zeros(length(x0))

    # inicializa o vetor de tensões equivalentes 
    Λ_tio = zeros(length(x0))

    # dicionário para armazenar os dados das seções transversais, evitando ler o mesmo arquivo várias vezes
    cache_secoes = Dict()
    hist = nothing

    # Chao Le
    ChaoLe = ones(2*malha.ne)

    # Incializa uma penalização inicial para usar no SLP
    penal = 1.0

    # Loop externo de otimização 
    for iter=1:niter

        # Calcula frequências e modos
        ωn,U0,_ = Modal3D(malha,posfile,x0=x0,kparam=[fkparam],mparam=[fmparam],n_modos=n_modos,iter=iter)
        
        # Deslocamentos da malha 1
        U, = Analise3D(malha, posfile, x0=x0, kparam=[fkparam],iter=iter)

        ## arquivo de esforços da iteracao
        arquivoEsf = nomeEsf * "_iter$(iter).esf"

        # tensão equivalente dos nos e elementos 
        # Λ_tio - vetor de tensões equivalentes de cada nó da seção
        # S - matriz s da seção
        # Pi - matriz Pi de cada nó da seção
        # σe_tio estado de tensão (σxx, σxy) para cada nó da seção
        Λ_tio,S,Pi,σe_tio = tensoes(arquivoEsf,malha,iter,false,cache_secoes)

        # so para armanezar a primeira frequencia da primeira iteração
        if iter == 1
            ω1 = ωn[1]
        end

        # Deriva da norma da frequencia - valor mínimo
        dω = norma_dω(ωn,U0,malha,x0,fdkparam,fdmparam,fmparam,P) 
    
        # derivada da tensao em relacao as variaveis de projeto
        dσ = ChaoLe .* norma_dσ(Λ_tio,σe_tio,S,Pi,malha,U,x0,fkparam,fdkparam,fσparam,fdσparam,P,s,σesc)

        # Determina os limites móveis, baseados nas variações das
        # variáveis de projeto. Isso só faz sentido para iter > 2
        atualiza_δ!(iter,δ,Δ1,Δ2,δ_min,δ_max)

        # Lineariza o problema 
        c,A,b,xi,xs,n,m = Lineariza(x0, δ, dω, dV,V_sup,x_inf,x_sup,dσ,Λ_tio,σesc,P,s,ChaoLe)

        #@show (A*x0 .- b)

        # Chama a solução interna do problema
        xn, violacoes ,penal = LP(c,A,b,xi,xs,n,m,penal)

        # Calcula o valor da função objetivo pela norma
        ωxn = norm(ωn,-P)
    
        # Teste de convergência do problema 
        #if iter > 2
        #    if convergencia(xn,ωx0,ωxn,A,b,tol_f,tol_g)
        #        break
        #    end
        #end

        # Copia o valor da funcao 
        ωx0 = copy(ωxn)

        # Roda e apita
        x2 .= x1
        x1 .= x0
        x0 .= xn

        # usando os valores das 3 ultimas iterações
        Δ2 = x1 .- x2
        Δ1 = x0 .- x1

        # Escreve no arquivo de saída
        println("iter $iter")
        println("frequencia ", ωn[1])
        println("x   ",xn)
        println("δ   ",δ)
        nσ = length(Λ_tio)   
        for i in 1:nσ
            Fs = σesc /maximum(Λ_tio[i])
            println("FS:$i ", Fs)
        end
        println("------------------------------")

        if hist === nothing
            hist = HistoricoOtim(length(ωn), ne, nσ)
        end

        push!(hist.iters, iter)
        hist.freq       = vcat(hist.freq,      ωn')
        hist.volume     = vcat(hist.volume,    [xn' * V])
        hist.densidades = vcat(hist.densidades, xn')

        for i in 1:nσ
            push!(hist.FS[i], σesc / maximum(Λ_tio[i]))
        end


    end # loop externo

    # Volume da estrutura final
    Vfinal = xn'*V
 
    # Recalcula frequências e modos para a solução final
    ωn, U0, _ = Modal3D(malha, posfile, x0=xn, kparam=[fkparam], mparam=[fmparam], n_modos=n_modos, iter=0)
 
    println()
    println("╔══════════════════════════════════════════════════╗")
    println("║           RESULTADOS DA OTIMIZAÇÃO               ║")
    println("╠══════════════════════════════════════════════════╣")
    println("║  Volume inicial  : ", lpad(round(V0,     digits=6), 12), " m³         ║")
    println("║  Volume final    : ", lpad(round(Vfinal, digits=6), 12), " m³         ║")
    println("║  Fração de volume: ", lpad(round(Vfinal/V0*100, digits=2), 11), " %          ║")
    println("╠══════════════════════════════════════════════════╣")
    println("║  1ª frequência inicial  : ", lpad(round(ω1,    digits=4), 10), " rad/s    ║")
    println("║  1ª frequência otimizada: ", lpad(round(ωn[1], digits=4), 10), " rad/s    ║")
    println("╠══════════════════════════════════════════════════╣")
    println("║  Densidades relativas dos elementos:             ║")
    for i in eachindex(xn)
        println("║    Elemento $i: ", lpad(round(xn[i], digits=6), 10), "                    ║")
    end
    println("╠══════════════════════════════════════════════════╣")
    println("║  Fatores de segurança (tensão):                  ║")
    nσ = size(dσ, 1)
    for i in 1:nσ
        Fs = σesc / maximum(Λ_tio[i])
        println("║    FS[$i]: ", lpad(round(Fs, digits=4), 12), "                       ║")
    end
    println("╚══════════════════════════════════════════════════╝")
    println()
 
    return ωn, U0,hist
end







