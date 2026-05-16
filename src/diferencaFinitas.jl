"""
    valida_dσ_FD(malha, x0, fkparam, fdkparam, P, iter;
                 h=1e-6, posfile=false)

Valida a derivada analítica `norma_dσ` usando diferenças finitas centrais.
"""
function valida_dσ_FD(malha, x0, fkparam, fdkparam, P, iter;
                      h=1e-6, posfile=true)

    ne = length(x0)

    # Solução no ponto base
    U, = Analise3D(malha, posfile, x0=x0,kparam=[fkparam], iter=iter)

    nomeEsf = malha.nome_arquivo
    arquivoEsf = nomeEsf * "_iter$(iter).esf"

    σeqMaxima, SS, Pi, σe = tensoes(arquivoEsf, malha, P, iter, posfile)
    @show σeqMaxima, SS, Pi, σe
    # Derivada analítica
    dσ_analitica = norma_dσ(σeqMaxima, σe, SS, Pi,malha, U, x0,fkparam, fdkparam, P)

    # Diferenças finitas centrais
    dσ_fd = zeros(Float64, ne)

    for i in 1:ne

        x_plus  = copy(x0)
        x_minus = copy(x0)

        x_plus[i]  += h
        x_minus[i] -= h

        # f(x+h)
        U_p, = Analise3D(malha, posfile,x0=x_plus,kparam=[fkparam],iter=iter)

        σeqMaxima_p, _, _, _ = tensoes(arquivoEsf, malha, P, iter, posfile)

        f_plus = norm(σeqMaxima_p, P)

        # f(x-h)
        U_m, = Analise3D(malha, posfile,x0=x_minus,kparam=[fkparam],iter=iter)

        σeqMaxima_m, _, _, _ = tensoes(arquivoEsf, malha, P, iter, posfile)

        f_minus = norm(σeqMaxima_m, P)

        # Diferença finita central
        dσ_fd[i] = (f_plus - f_minus) / (2h)
    end

    println("\n=== Validação de norma_dσ (Diferenças Finitas Centrais) ===")
    println(@sprintf("%-6s  %-14s  %-14s  %-12s",
                     "elem", "analítica", "DF central", "erro rel."))
    println("-" ^ 58)

    for i in 1:ne
        a  = dσ_analitica[i]
        fd = dσ_fd[i]

        er = abs(a - fd) / (abs(fd) + 1e-30)

        flag = er > 1e-4 ? "  ← ERRO" : ""

        println(@sprintf("%-6d  %-14.6e  %-14.6e  %-12.2e%s",i, a, fd, er, flag))
    end

    return dσ_analitica, dσ_fd
end

function valida_dω_FD(malha, x0,fkparam, fdkparam,fmparam, fdmparam,n_modos, P;h=1e-4, posfile=false)

    ne = length(x0)

    # Solução base
    ωn, U0, _ = Modal3D(malha, posfile,x0=x0,kparam=[fkparam],mparam=[fmparam])

    # recupera os modos
    ωnn = ωn[1:n_modos]
    U0n = U0[:, 1:n_modos]

    # Derivada utilizando a norma P
    dω_analitica = norma_dω(ωnn, U0n,malha, x0,fdkparam, fdmparam,P)

    # Diferenças finitas centrais
    dω_fd = zeros(Float64, ne)

    for i in 1:ne

        x_plus  = copy(x0)
        x_minus = copy(x0)

        x_plus[i]  += h
        x_minus[i] -= h

        # f(x+h)
        ωn_p, _, _ = Modal3D(malha, posfile,x0=x_plus,kparam=[fkparam],mparam=[fmparam])

        # frequência inferior dos modos no ponto a frente 
        f_plus = norm(ωn_p[1:n_modos], -P)

        # f(x-h)
        ωn_m, _, _ = Modal3D(malha, posfile,x0=x_minus,kparam=[fkparam],mparam=[fmparam])

        # frequência inferior dos modos no ponto a trás
        f_minus = norm(ωn_m[1:n_modos], -P)

        # Diferença finita central
        dω_fd[i] = (f_plus - f_minus) / (2h)
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