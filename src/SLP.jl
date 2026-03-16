###########################################################################
#            PROCEDIMENTO DE SOLUÇÂO PROBLEMA-AGNÓSTICO 
###########################################################################

#
#   Função que lineariza o problema
#
function Lineariza(x::Vector, δ::Vector, x_inf::Vector,x_sup::Vector,
                   Objetivo::Function, Restricoes::Function,
                   dObjetivo::Function, dRestricoes::Function)


    # Calcula as restrições no ponto e já descobre o número de restrições
    g_real = Restricoes(x)

    # número de restrições
    m = length(g_real)
    
    # número de variáveis de projeto
    n = length(x)

    # Aloca os arrays do LP
    c = zeros(n)
    A = zeros(m,n)
    b = zeros(m)
    xi = zeros(n)
    xs = zeros(n)

    # Coeficientes da função objetivo linearizada
    c = dObjetivo(x)

    # Gradiente das restrições linearizadas
    A = dRestricoes(x)

    # rhs das restrições linearizadas
    for i=1:m
        b[i] = dot(A[i,:],x) - g_real[i]
    end

    # Limites móveis do problema
    for i = 1:n

        if x[i]>0

            xi[i] = max((1 - δ[i])*x[i], x_inf[i])
            xs[i] = min((1 + δ[i])*x[i], x_sup[i])

        else

            xi[i] = max(x[i] - δ[i], x_inf[i])
            xs[i] = min(x[i] + δ[i], x_sup[i])

        end # if

    end # for i

    # Retorna os Coeficientes do problema linearizado
    return c,A,b,xi,xs,n,m

end

#
#   Soluciona o LP (interno)
#
function LP(c,A,b,xi,xs,n,m)

    
    # Aloca vetor de soluções atuais
    xn = zeros(n)

    # Cria o modelo linear atual
    modelo = Model(HiGHS.Optimizer)

    # cale-se cale-se cale-se
    set_silent(modelo)

    # Define as variáveis de projeto
    @variable(modelo, xi[i] <= x[i = 1:n] <= xs[i])

    # Define as variáveis de folga
    @variable(modelo, s[1:m] >= 0)

    # Define a penalização
    r = 1E3

    # Define a função objetivo
    @objective(modelo, Min, c'*x + r*sum(s))

    # Agora com a restrição nl linearizada com variáveis de folga
    @constraint(modelo, g1, A*x + s .<= b)

    # Soluciona o sub-problema
    optimize!(modelo)

    for i = 1:n

        # Recupera a solução atual
        xn[i] = value(x[i])

    end # for i

    return xn,(A*xn .- b)
    
end




###########################################################################
#            PROGRAMA PRINCIPAL - RECEBE AS DEFINIÇÕES DO PROBLEMA 
#                  E EXECUTA TODOS OS PASSOS DE SOLUÇÃO
###########################################################################

#
#   Otimização sequencial
#
function OtimizaSLP(x0::Vector, x_inf::Vector, x_sup::Vector)

    # Inicializa os vetores de iterações anteriores
    x1 = copy(x0)
    x2 = copy(x1)

    # Só inicializa aqui...
    Δ1 = x0 .- x1
    Δ2 = x1 .- x2

    # Vamos inicializar o veltor de deltas
    δ = 0.1*ones(2)

    # Limites máximos e mínimos para os deltas
    δ_max = 0.2
    δ_min = 0.01

    # Abre arquivo para gravar o histórico
    fd = open("convergencia_LP.txt","w")

    # Loop de otimização sequencial
    for iter = 1:100

        # Determina os limites móveis, baseados nas variações das
        # variáveis de projeto. Isso só faz sentido para iter > 2

        if iter > 2

            # Loop pelas variáveis de projeto
            for i = 1:length(x0)

                # identifica o zig-zag
                if Δ1[i]*Δ2[i] < 0 

                    # Diminui o multiplicador do limite móvel
                    δ[i] = max(δ[i]*0.7, δ_min)

                else

                    # aumenta o multiplicador do limite móvel
                    δ[i] = min(δ[i]*1.2, δ_max)

                end # if
            end # for i
        end # if iter > 2

        # Lineariza o problema 
        c,A,b,xi,xs,n,m = Lineariza(x0, δ, x_inf, x_sup, Objetivo, Restricoes,dObjetivo, dRestricoes)
 
        # Chama a solução interna do problema
        xn, gs_lin = LP(c,A,b,xi,xs,n,m)

        # Define as tolerâncias
        tol_f = 1E-4 # se tol_f menor que 1E-3 para por limite de iterações
        tol_g = 1E-4

        # Calcula o valor da função objetivo em x0 e xn
        f_x0 = Objetivo(x0)
        f_xn = Objetivo(xn)

        # Calcula a diferença relativa das variáveis de projeto
        # Calcula a diferença relativa da função objetivo
        dif_x = (norm(xn - x0))/(norm(x0) + 1E-6)
        dif_fx = (abs(f_xn - f_x0))/(abs(f_x0) + 1E-6)

        # Calcula g_real e o termo violation
        g_real = Restricoes(xn)
        violation = max(0.0, minimum(g_real))

        # Verifica as tolerâncias
        if dif_fx < tol_f && violation < tol_g

            # Avisa convergencia
            println(fd,"Convergiu")
            println(fd,"iter $iter")
            println(fd,"x   ",xn)
            println(fd,"δ   ",δ)
            println(fd,"gs_lin   ",gs_lin)
            println(fd,"dif_x   ",dif_x)
            println(fd,"dif_fx   ",dif_fx)
            println(fd,"g_real   ",g_real)
            println(fd,"violação   ",violation)
            println(fd,"------------------------------")

            close(fd)

            x0 .= xn

            break

        end

        # Roda e apita
        x2 .= x1
        x1 .= xn
        x0 .= xn

        # Escreve no arquivo de saída
        println(fd,"iter $iter")
        println(fd,"x   ",xn)
        println(fd,"δ   ",δ)
        println(fd,"gs_lin   ",gs_lin)
        println(fd,"dif_x   ",dif_x)
        println(fd,"dif_fx   ",dif_fx)
        println(fd,"g_real   ",g_real)
        println(fd,"violação   ",violation)
        println(fd,"------------------------------")

    end # loop for iter
    
    # Fecha o arquivo de saída
    close(fd)

    return x0
end