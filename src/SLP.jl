###########################################################################
#            PROCEDIMENTO DE SOLUÇÂO PROBLEMA-AGNÓSTICO 
###########################################################################

#
#   Função que lineariza o problema
#
function Lineariza(x0, δ, dω,dV,V_sup,x_inf,x_sup,dσ,Λ,σesc,P,s,ChaoLe)

    # número de restrições uma de volume e nσ de tensão
    nσ = size(dσ,1)
    m = 1 + nσ
    
    # número de variáveis de projeto
    n = length(x0)

    # Aloca os arrays do LP
    c = zeros(n)
    A = zeros(m,n)
    b = zeros(m)

    # Coeficientes da função objetivo linearizada
    # Lembrando que queremos MAXIMIZAR a menor 
    # frequência (obtida pela norma)
    c = dω

    # Gradiente das restrições linearizadas
    A[1,:] = dV'/V_sup

    # a derivada da tensao ja leva em consideracao o fator s/σesc
    A[2:end,:] = dσ

    # aq eu mudei estava sum(x0.*dV)/V_sup - 1.0 e nao estava indo.
    # rhs das restrições linearizadas - a restrição de volume é linear e já tem o rhs definido como V_sup
    b[1] = 1.0

    for i in 1:nσ
        
        σi = ChaoLe[i]*norm(Λ[i], P) / (σesc/s)
        gi = σi - 1.0

        b[i+1] = -gi + dot(dσ[i,:], x0)

        ChaoLe[i] = maximum(Λ[i])/norm(Λ[i],P)
    end

    # inicializa xs e xi
    xi = zeros(n)
    xs = zeros(n)

    # Limites móveis do problema
    for i = 1:n

        # sempre será maior que 0 , pois x : [0,1]
        xi[i] = max((1 - δ[i])*x0[i], x_inf[i])
        xs[i] = min((1 + δ[i])*x0[i], x_sup[i])

    end

    # Retorna os Coeficientes do problema linearizado
    return c,A,b,xi,xs,n,m

end

#
#   Soluciona o LP (interno)
#
function LP(c,A,b,xi,xs,n,m,penal)

    # Checa antes de passar pro solver
    @assert all(isfinite, c)  "c tem valor não-finito: $c"
    @assert all(isfinite, A)  "A tem valor não-finito"
    @assert all(isfinite, b)  "b tem valor não-finito: $b"
    @assert all(isfinite, xi) "xi tem valor não-finito: $xi"
    @assert all(isfinite, xs) "xs tem valor não-finito: $xs"

    # Checa bounds invertidos
    for i in 1:n
        if xi[i] > xs[i]
            error("Bounds invertidos na variável $i: xi=$(xi[i]) > xs=$(xs[i]), x0 provavelmente é zero")
        end
    end

    # Aloca vetor de soluções atuais
    xn = zeros(n)

    # Cria o modelo linear atual
    modelo = Model(HiGHS.Optimizer)

    # cale-se cale-se cale-se
    set_silent(modelo)

    #set_attribute(modelo, "log_file", "meleca")

    # Define as variáveis de projeto
    @variable(modelo, xi[i] <= x[i = 1:n] <= xs[i])

    # Define as variáveis de folga
    @variable(modelo, s[1:m] >= 0)

    # Define a penalização mínima
    r = 10.0

    # Define a função objetivo
    @objective(modelo, Min, -c'*x/norm(c) + penal*sum(s))

    # Agora com a restrição nl linearizada com variáveis de folga
    # vou guardar como uma variável para poder pegar os duais depois 
    restricoes = @constraint(modelo, [i=1:m], dot(A[i,:], x) - s[i] <= b[i])

    # Soluciona o sub-problema
    optimize!(modelo)
    status = termination_status(modelo)

    if status != MOI.OPTIMAL
        error("LP não convergiu. Status = $status")
    end

    for i = 1:n

        # Recupera a solução atual
        xn[i] = value(x[i])
        
    end # for i

    for i in 1:m
        @show value(s[i])
    end

    # recupera os multiplicadores de lagrange
    λ = dual.(restricoes)

    # a penalização vai ser o máximo entre a nossa penalização "fixa" e 
    # o maior λ
    penal = max(r,maximum(abs.(λ)))

    @show λ, penal

    return xn,(A*xn .- b), penal
    
end


#
# Calcula o ajuste de limites móveis para a iteração
#
function atualiza_δ!(iter::Int,δ::Vector{T},Δ1::Vector{T},Δ2::Vector{T},δ_min::T,δ_max::T) where T

    # Determina os limites móveis, baseados nas variações das
    # variáveis de projeto. Isso só faz sentido para iter > 2

    if iter > 2

        # Loop pelas variáveis de projeto
        for i in eachindex(δ)

            # identifica o zig-zag
            if Δ1[i]*Δ2[i] < 0 

                # Diminui o multiplicador do limite móvel
                δ[i] = max(δ[i]*0.7, δ_min)

            else

                # aumenta o multiplicador do limite móvel
                δ[i] = min(δ[i]*1.2, δ_max)

            end 
        end 
        
    end 


end

#
# Usar a norma P da frequência - objetivo
#
function convergencia(xn,ωx0,ωxn,A,b,tol_f,tol_g)

    # Calcula a diferença relativa da função objetivo
    dif_fx = norm(ωxn - ωx0)/(norm(ωx0))

    violation = maximum(A*xn .- b)
     
    # Verifica as tolerâncias
    if dif_fx < tol_f && violation < tol_g
        return true
    else 
        return false
    end

end