###########################################################################
#            PROCEDIMENTO DE SOLUÇÂO PROBLEMA-AGNÓSTICO 
###########################################################################

#
#   Função que lineariza o problema
#
function Lineariza(x0, δ, dω,dV,V_sup,x_inf,x_sup)

    # número de restrições por enquanto 1 só dv
    m = 1
    
    # número de variáveis de projeto
    n = length(x0)

    # Aloca os arrays do LP
    c = zeros(n)
    A = zeros(m,n)
    b = zeros(m)

    # Coeficientes da função objetivo linearizada
    # Lembrando que queremos MAXIMIZAR a menor 
    # frequência (obtida pela norma)
    c = -dω

    # Gradiente das restrições linearizadas
    A = dV'

    # rhs das restrições linearizadas
    #for i=1:m
        #b[i] = V_sup - dV[i] + dot(A[i,:],x0) 
    #end

    # como tem só uma restrição e ela já é linear (adaptar para quando não for usando o for acima)
    b[1] = V_sup

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
    @constraint(modelo, [i=1:m], dot(A[i,:], x) + s[i] <= b[i])

    # Soluciona o sub-problema
    optimize!(modelo)

    for i = 1:n

        # Recupera a solução atual
        xn[i] = value(x[i])

    end # for i

    return xn,(A*xn .- b)
    
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
function convergencia(x0,xn,ωx0,ωxn,dV,V_sup,tol_f,tol_g)

    # Calcula a diferença relativa da função objetivo
    dif_fx = norm(ωxn - ωx0)/(norm(ωx0))

    # Calcula o termo violation
    violation = max(0.0, dot(dV, xn) - V_sup)
     
    # Verifica as tolerâncias
    if dif_fx < tol_f && violation < tol_g
        return true
    else 
        return false
    end

    

end