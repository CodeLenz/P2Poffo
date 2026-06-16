"""
    valida_dσ_FD(malha, x0, fkparam, fdkparam, P, iter;
                 h=1e-6, posfile=false)

Valida a derivada analítica `norma_dσ` usando diferenças finitas centrais.
"""

function valida_dσ_FD(arquivo::AbstractString, fkparam=P2Poffo.Kparam, fdkparam=P2Poffo.dKparam,
                      fσparam=P2Poffo.σparam, fdσparam=P2Poffo.dσparam,
                      P=8.0, iter=4; h=1e-3, s=2.0, σesc=150e6, posfile=true)


    ωn,U0,malha =Modal3D(arquivo)

    # densidade dos elementos
    x0 = 0.5 * ones(malha.ne) + 1E-2*randn(malha.ne)

    # Dimensão do vetor de variáveis de projeto
    ne = malha.ne
    
    # Solução no ponto base
    U,= Analise3D(malha, posfile, x0=x0,kparam=[fkparam], iter=iter)

    # Carrega arquivo de esforços
    nomeEsf = malha.nome_arquivo
    arquivoEsf = nomeEsf * "_iter$(iter).esf"
    cache_secoes = Dict()
    # Calcula as tensões de referência 
    Λ_tio, S, Pi, σe_tio = tensoes(arquivoEsf,malha,iter,posfile,cache_secoes)

    # Derivada analítica
    dσ_analitica = norma_dσ(Λ_tio, σe_tio, S, Pi, malha, U, x0, fkparam, fdkparam, fσparam, fdσparam, P, s, σesc)

    ncomp = 2 * malha.ne
    dσ_fd = zeros(Float64, ncomp, ne)

    for i in LinearIndices(x0)
        xb = x0[i]

        x0[i] = xb + h
        U_p, = Analise3D(malha, posfile, x0=x0, kparam=[fkparam], iter=5)
        Λ_p_tio, _, _, _ = tensoes(nomeEsf * "_iter5.esf", malha, 5, posfile,cache_secoes)

        x0[i] = xb - h
        U_m, = Analise3D(malha, posfile, x0=x0, kparam=[fkparam], iter=6)
        Λ_m_tio, _, _, _ = tensoes(nomeEsf * "_iter6.esf", malha, 6, posfile,cache_secoes)

        x0[i] = xb

        xtest_p = copy(x0); xtest_p[i] = xb + h
        xtest_m = copy(x0); xtest_m[i] = xb - h

        Λ_p = similar(Λ_p_tio)
        Λ_m = similar(Λ_m_tio)

        for idx in 1:ncomp
            ele = (idx - 1) ÷ 2 + 1
            fσ_p = fσparam(xtest_p[ele])
            fσ_m = fσparam(xtest_m[ele])
            Λ_p[idx] = fσ_p .* Λ_p_tio[idx]
            Λ_m[idx] = fσ_m .* Λ_m_tio[idx]
            σ_p_idx = (s/σesc) * norm(Λ_p[idx], P)
            σ_m_idx = (s/σesc) * norm(Λ_m[idx], P)
            dσ_fd[idx, i] = (σ_p_idx - σ_m_idx) / (2h)
        end

         
    end

    println("\n================== Validação de norma_dσ ==================")
    println("\n2 elementos × 2 nós = 4 restrições")
    println("\n4 restrições × 2 variáveis(elementos) = 8 derivadas\n") 
    @printf("%-14s %-18s %-18s %-14s\n","(rest,vari)", "Analítica", "DF Central", "Erro Rel.")
    println("---------------------------------------------------------------")

    for j in 1:ncomp
        for i in 1:ne

            a  = dσ_analitica[j,i]
            fd = dσ_fd[j,i]

            er = abs(a - fd) / (abs(fd) + 1e-12)

            # evita falso erro perto de zero
            er = (abs(a) < 1e-7 && abs(fd) < 1e-7) ? 0.0 : er

            @printf("(%2d,%2d)  % .6e   % .6e   % .3e     \n",
                    j, i, a, fd, er)
        end
        println("---------------------------------------------------------------")
    
    end

    return dσ_analitica, dσ_fd
end

function valida_dω_FD(malha, x0,fkparam=P2Poffo.Kparam, fdkparam=P2Poffo.dKparam, fmparam=P2Poffo.Mparam, fdmparam=P2Poffo.dMparam,n_modos=3, P=8.0;h=1e-4, posfile=false)

    # Número de variáveis de projeto 
    ne = length(x0)

    # Solução base
    ωn, U0, _ = Modal3D(malha, posfile,x0=x0,kparam=[fkparam],mparam=[fmparam])

    # recupera os modos
    ωnn = ωn[1:n_modos]
    U0n = U0[:, 1:n_modos]

    # Derivada utilizando a norma P
    dω_analitica = norma_dω(ωnn, U0n,malha, x0,fdkparam, fdmparam,fmparam, P)

    # Diferenças finitas centrais
    dω_fd = zeros(Float64, ne)

    # Loop pelas variáveis de projeto 
    for i in LinearIndices(x0)

        # Backup do valor atual 
        xb = x0[i]

        # Perturba para frente
        x0[i] += h

        # Calcula o problema para frente
        ωn_p, _, _ = Modal3D(malha, posfile,x0=x0,kparam=[fkparam],mparam=[fmparam])

        # frequência inferior dos modos no ponto a frente 
        f_plus = norm(ωn_p[1:n_modos], -P)

        # Perturba para trás 
        x0[i] = xb - h

        # frequências para trás
        ωn_m, _, _ = Modal3D(malha, posfile,x0=x0,kparam=[fkparam],mparam=[fmparam])

        # frequência inferior dos modos no ponto a trás
        f_minus = norm(ωn_m[1:n_modos], -P)

        # Diferença finita central
        dω_fd[i] = (f_plus - f_minus) / (2h)

        # Restaura a posição i 
        x0[i] = xb

    end

    println("\n=== Validação de norma_dω (Diferenças Finitas Centrais) ===")
    println(@sprintf("%-6s  %-14s  %-14s  %-12s",
                     "elem", "Norma P", "DF central", "erro rel."))
    println("-" ^ 58)

    for i in 1:ne
        
        a  = dω_analitica[i]
        fd = dω_fd[i]

        er = abs(a - fd) / (abs(fd) + 1e-30)

        flag = er > 1e-4 ? "  ← ERRO" : ""

        println(@sprintf(
            "%-6d  %-14.6e  %-14.6e  %-12.2e%s",
            i, a, fd, er, flag
        ))
    end

    return dω_analitica, dω_fd
end