module Torcao

    using LinearAlgebra
    using SparseArrays
    using Lgmsh

    # Carrega as rotinas
    include("conversor.jl")
    include("material.jl")
    include("elemento.jl")
    include("global.jl")
    include("apoios.jl")
    include("corpo.jl")
    include("jeq.jl")
    include("momentosinercia.jl")
    include("main.jl")
    include("tensoes.jl")

    # Exporta o arquivo .dat
    include("exportaarquivo.jl")
    
	# Exporta a rotina principal de análise
	export AnaliseTorcao

    # Exporta a rotina que cria o .dat
    export ExportaDat

	
end
