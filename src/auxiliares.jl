#
# Retorna um vetor com o volume de cada elemento da malha para otimização
#
function Volumes(malha::LFrame.Malha)

    # Recupera os dados da malha
    L = malha.L 
    dados_ele = malha.dados_elementos
    dicionario_mat = malha.dicionario_materiais
    dicionario_geo = malha.dicionario_geometrias


    # Inicializa o vetor de volumes
    V = zeros(malha.ne)

    # Loop pelos elementos
    for ele = 1:malha.ne

        Iz, Iy, J0, A, α, E, G, geo = LFrame.Dados_fundamentais(ele, dados_ele, dicionario_mat, dicionario_geo)

        V[ele] = A * L[ele]
    end

    # Retorna o volume
    return V

end

#
# Critério de ótimo (loop local)
#
function OC(x0::Vector, dC::Vector, dV::Vector, V_sup::Float64, ne::Int64, μ1=0.0, μ2=1E5, δ=0.1, tol=1E-6,x_min = 1E-3)

    # Copia do x0
    x_estimado = copy(x0)

    # Aqui, V e dV são iguais
    V = dV

    # Loop da bisseção
    for k=1:1000  
       
       # Ponto médio do Intervalo 
       μ = (μ1 + μ2)/2
          
        # Loop pelos elementos
        for i=1:ne     
            
            # fator beta(aqui tinha um min)
            β = -dC[i]/(μ*dV[i])
        
            # Novo valor x estimado 
            x_e = x0[i] * (β^0.5)
            
            # Limites
            x_dir = min(x0[i] + δ, 1.0) 
            x_esq = max(x0[i] - δ, x_min) 
           
            # 
            #-------0[----|--x--|----]1---------  
            #           x-δ    x+δ 
            #
            
            # Verificando os Limites
            x_estimado[i] = max( min( x_e, x_dir), x_esq)

        end 

       # Volume novo
       V_atual = sum(x_estimado .* V)
 
       # Testando se atingiu a tolerancia
       if abs(μ2 - μ1)<=tol
            break
       end
        
       # Verificando se o volume atingiu o máximo
       if (V_atual > V_sup)
            μ1 =  μ
       else
            μ2 =  μ
       end
 
        
    end # k
 
    # Return o x_estimado
    return x_estimado


end

function precomputa_geometria(dados_secao)
    nn      = dados_secao.nn
    XY      = dados_secao.XY
    α       = dados_secao.α
    centroide = dados_secao.centroide
    ∇Φ      = dados_secao.∇Φ
    vizinhos = dados_secao.vizinhos
    area    = dados_secao.area
    Izl     = dados_secao.Izl
    Iyl     = dados_secao.Iyl
    Jeq     = dados_secao.Jeq

    zl_yl   = zeros(nn, 2)
    grad_xy = zeros(nn, 2)
    grad_n  = zeros(nn)
    # Pi_base[i] = Pi sem o αμ e sem o sign(T), separado em duas partes:
    #   Pi[i] = [1/area  0        -yl/Izl  zl/Iyl  ]
    #           [0       sign(T)*grad_n[i]/Jeq  0  0]
    # Guardamos a parte geométrica pura (sem sign(T) e sem αμ):
    Pi_geo = Vector{Matrix{Float64}}(undef, nn)

    for ino in 1:nn
        x, y = XY[ino, 1:2]
        zl, yl = Muda_coordenada(x, y, α, centroide[1], centroide[2])
        zl_yl[ino, :] = [zl, yl]

        vizi = vizinhos[ino]
        ∇xΦ = mean(∇Φ[vizi, 1])
        ∇yΦ = mean(∇Φ[vizi, 2])
        grad_xy[ino, :] = [∇xΦ, ∇yΦ]
        grad_n[ino]     = sqrt(∇xΦ^2 + ∇yΦ^2)

        Pi_geo[ino] = [1/area    0              -yl/Izl   zl/Iyl;
                       0         grad_n[ino]/Jeq  0        0     ]
    end

    return zl_yl, grad_xy, grad_n, Pi_geo
end