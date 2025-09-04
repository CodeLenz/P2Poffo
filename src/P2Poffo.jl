module P2Poffo

    using LinearAlgebra
    using SparseArrays
    using StaticArrays
    using Lgmsh

    # Carrega as rotinas
    include("conversor.jl")
    include("material.jl")
    include("elemento.jl")
    include("global.jl")
    include("apoios.jl")
    include("corpo.jl")
    include("main.jl")
    include("tensoes.jl")

    # Alguns valores para referência
    include("referencias.jl")

    # Exporta o arquivo .dat
    include("exportaarquivo.jl")
    
	# Exporta a rotina principal de análise
	export AnaliseTorcao

    # Exporta a rotina que cria o .dat
    export ExportaDat

	
end
