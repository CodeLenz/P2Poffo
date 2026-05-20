"""
    valida_dσ_FD(malha, x0, fkparam, fdkparam, P, iter;
                 h=1e-6, posfile=false)

Valida a derivada analítica `norma_dσ` usando diferenças finitas centrais.
"""
function valida_dσ_FD(malha, x0, fkparam, fdkparam, P, iter;
                      h=1e-6, posfile=true)

    # Dimensão do vetor de variáveis de projeto
    ne = length(x0)

    # Solução no ponto base
    U, = Analise3D(malha, posfile, x0=x0,kparam=[fkparam], iter=iter)

    # Carrega arquivo de esforços
    nomeEsf = malha.nome_arquivo
    arquivoEsf = nomeEsf * "_iter$(iter).esf"

    # Calcula as tensões de referência 
    σeqMaxima, SS, Pi, σe = tensoes(arquivoEsf, malha, P, iter, posfile)

    #return typeof(σeqMaxima),typeof(σe),typeof(SS),typeof(Pi),typeof(malha),typeof(U),typeof(x0),typeof(fkparam),typeof(fdkparam),typeof(P)
    # Derivada analítica
    dσ_analitica = norma_dσ(σeqMaxima, σe, SS, Pi,malha, U, x0,fkparam, fdkparam, P)

    # Diferenças finitas centrais
    dσ_fd = zeros(Float64, ne)

    # Loop pelas variáveis de projeto
    for i in LinearIndices(x0)

        # Backup do valor atual 
        xb = x0[i]

        # Perturba para frente
        x0[i] += h
       
        # f(x+h) aq iteraçao 2 é um maneira de ele ler o arquivo de esforços gerado na frente
        U_p, = Analise3D(malha, posfile,x0=x0,kparam=[fkparam],iter=2)

        arquivoEsf_p = nomeEsf * "_iter2.esf"

        # Calcula as tensões com a perturbação para frente
        σeqMaxima_p, _, _, _ = tensoes(arquivoEsf_p, malha, P, iter, posfile)

        # Objetivo para frente
        f_plus = norm(σeqMaxima_p, P)

        # Perturba para trás 
        x0[i] = xb - h

        # Análise para trás
        U_m, = Analise3D(malha, posfile,x0=x0,kparam=[fkparam],iter=3)

        arquivoEsf_m = nomeEsf * "_iter3.esf"

        # Tensões para trás
        σeqMaxima_m, _, _, _ = tensoes(arquivoEsf_m, malha, P, iter, posfile)

        # Função objetivo para trás
        f_minus = norm(σeqMaxima_m, P)

        # Diferença finita central
        dσ_fd[i] = (f_plus - f_minus) / (2h)

        # Restaura a posição i 
        x0[i] = xb

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

    # Número de variáveis de projeto 
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