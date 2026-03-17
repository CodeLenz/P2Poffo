#
# Rotinas de integração
#
#
# V(x) <= V_sup = vf * V0
# 

# Parametrização SIMP 
simp(x,p=3) = x^p

# Derivada da parametrização 
dsimp(x,p=3) = p*x^(p-1)

# Parametrização SIMP 
gimp(x,p=1) = x^p

# Derivada da parametrização 
dgimp(x,p=1) = p*x^(p-1)

#
# Rotina principal
#
function Main_Otim_OC(arquivo::AbstractString, fkparam::Function, fdkparam::Function, posfile=true; verbose=false,vf = 0.5, niter=100)

    # Chama o Analise3D com o nome do arquivo, para receber a estrutura de malha
    U0, malha = Analise3D(arquivo,posfile,verbose=verbose)

    # Numero de elementos 
    ne = malha.ne

    # Estimativa inicial das densidades relativas
    x0 = vf*ones(ne)

    # Calcula o volume de cada elemento sem considerar a 
    # parametrização 
    V = Volumes(malha)

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


#
# Isso aqui vai ter que juntar com o OtimizaSLP(x0,dω,dV,V_sup,ne) <- recebe as coisas gerais, 
# pois tu vais ter que calcular as derivadas e funções a cada iteração externa do SLP
#
function Main_Otim_Modal(arquivo::AbstractString, fkparam::Function, fdkparam::Function, posfile=true; verbose=false, vf = 0.5 , niter=100)

    # analise modal somente para a malha
    ωn,U0,malha = Modal3D(arquivo,posfile,verbose=verbose)

    # Numero de elementos 
    ne = malha.ne

    # Estimativa inicial das densidades relativas
    x0 = vf*ones(ne)

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

    # Vamos inicializar o veltor de deltas
    δ = 0.1*ones(ne)

    # Limites máximos e mínimos para os deltas
    δ_max = 0.2
    δ_min = 0.01

    # calcula a derivada no ponto x0 para copiar depois
    dω = norma_dω(ωn,U0,malha,x0,fdkparam,fdkparam)   

    # inicializando derivada normalizada no ponto xn
    dω_xn = copy(dω)

    # inicializa os deltas
    Δ1 = x0 .- x1
    Δ2 = x1 .- x2

    # inicializa os valores inferiores e superiores de densidades relativas
    x_inf = zeros(x0)
    x_sup = ones(x0)

    # Loop externo de otimização 
    for iter=1:niter

        # Calcula frequências e modos
        ωn,U0,_ = Modal3D(malha,posfile,x0=x0,kparam=[fkparam])
       
        # Deriva da norma da frequencia - valor mínimo
        dω = norma_dω(ωn,U0,malha,x0,fdkparam,fdkparam)   
         
        # Determina os limites móveis, baseados nas variações das
        # variáveis de projeto. Isso só faz sentido para iter > 2
        atualiza_δ!(iter,δ,Δ1,Δ2,δ_min,δ_max)

        # Lineariza o problema 
        c,A,b,xi,xs,n,m = Lineariza(x0, δ, dω,dV,V_sup,x_inf,x_sup)

        # Chama a solução interna do problema
        xn, gs_lin = LP(c,A,b,xi,xs,n,m)

        # Teste de convergência do problema 
        if convergencia(x0,xn,dω,dω_xn,dV)
            break
        end

        # Roda e apita
        x2 .= x1
        x1 .= x0
        x0 .= xn

        # usando os valores das 3 ultimas iterações
        Δ2 = x1 .- x2
        Δ1 = x0 .- x1

        # Atualiza os limites 
        x_inf = xi 
        x_sup = xs

        dω_xn = dω

        # Escreve no arquivo de saída
        println("iter $iter")
        println("x   ",xn)
        println("δ   ",δ)
        # println(fd,"gs_lin   ",gs_lin)
        # println(fd,"dif_x   ",dif_x)
        # println(fd,"dif_fx   ",dif_fx)
        # println(fd,"g_real   ",g_real)
        # println(fd,"violação   ",violation)
        #println(fd,"------------------------------")

        
    end # loop externo
    
    # print convergenciu 
    print("Convergiu")
end
    





