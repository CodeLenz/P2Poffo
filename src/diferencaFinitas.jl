"""
    valida_dσ_CS(arquivoEsf, malha, x0, U, fkparam, fdkparam, P, iter; h=1e-20)

Valida a derivada analítica `norma_dσ` usando complex step.
"""
function valida_dσ_CS(malha, x0, fkparam, fdkparam, P, iter;
                      h=1e-20, posfile=false)

    ne = length(x0)

    # Ponto base (real)
    U, = Analise3D(malha, posfile, x0=x0, kparam=[fkparam], iter=iter)
    nomeEsf = malha.nome_arquivo
    arquivoEsf = nomeEsf * "_iter$(iter).esf"
    σeqMaxima, SS, Pi, σe = tensoes(arquivoEsf, malha, P, iter, posfile)

    dσ_analitica = norma_dσ(σeqMaxima, σe, SS, Pi, malha, U, x0, fkparam, fdkparam, P)

    # Complex step
    dσ_cs = zeros(Float64, ne)

    for i in 1:ne
        xc = Complex{Float64}.(x0)
        xc[i] += im * h

        Uc, = Analise3D(malha, posfile, x0=xc, kparam=[fkparam], iter=iter)
        σeqMaxima_c, SS_c, Pi_c, σe_c = tensoes(arquivoEsf, malha, P, iter, posfile)

        # norma-P das tensões (função escalar)
        f_c = norm(σeqMaxima_c, P)   # ajuste para sua definição exata

        dσ_cs[i] = imag(f_c) / h
    end

    println("\n=== Validação de norma_dσ (Complex Step) ===")
    println(@sprintf("%-6s  %-14s  %-14s  %-12s", "elem", "analítica", "complex step", "erro rel."))
    println("-" ^ 52)

    for i in 1:ne
        a  = dσ_analitica[i]
        cs = dσ_cs[i]
        er = abs(a - cs) / (abs(cs) + 1e-30)
        flag = er > 1e-4 ? "  ← ERRO" : ""
        println(@sprintf("%-6d  %-14.6e  %-14.6e  %-12.2e%s", i, a, cs, er, flag))
    end

    return dσ_analitica, dσ_cs
end



"""
    valida_dω_CS(malha, x0, fkparam, fdkparam, fmparam, fdmparam, n_modos, P; h=1e-20)

Valida a derivada analítica `norma_dω` usando diferenças finitas complexas (complex step).
Compara elemento a elemento e imprime o erro relativo.
"""
function valida_dω_CS(malha, x0, fkparam, fdkparam, fmparam, fdmparam, n_modos, P;
                      h=1e-20, posfile=false)

    ne = length(x0)

    # Derivada analítica no ponto x0
    ωn, U0, _ = Modal3D(malha, posfile, x0=x0, kparam=[fkparam], mparam=[fmparam])
    ωnn = ωn[1:n_modos]
    U0n = U0[:, 1:n_modos]
    dω_analitica = norma_dω(ωnn, U0n, malha, x0, fdkparam, fdmparam, P)

    # Derivada por complex step
    dω_cs = zeros(Float64, ne)

    for i in 1:ne
        # Perturba somente o elemento i com parte imaginária
        xc = Complex{Float64}.(x0)
        xc[i] += im * h

        # Avalia a função objetivo com x complexo
        ωn_c, _, _ = Modal3D(malha, posfile, x0=xc, kparam=[fkparam], mparam=[fmparam])
        ωnn_c = ωn_c[1:n_modos]

        # Norma-P das frequências (valor complexo)
        f_c = norm(ωnn_c, -P)  # ou sua função norma_ω(ωnn_c, P)

        # Derivada = parte imaginária / h
        dω_cs[i] = imag(f_c) / h
    end

    # Comparação
    println("\n=== Validação de norma_dω (Complex Step) ===")
    println(@sprintf("%-6s  %-14s  %-14s  %-12s", "elem", "analítica", "complex step", "erro rel."))
    println("-" ^ 52)

    for i in 1:ne
        a  = dω_analitica[i]
        cs = dω_cs[i]
        er = abs(a - cs) / (abs(cs) + 1e-30)
        flag = er > 1e-4 ? "  ← ERRO" : ""
        println(@sprintf("%-6d  %-14.6e  %-14.6e  %-12.2e%s", i, a, cs, er, flag))
    end

    return dω_analitica, dω_cs
end