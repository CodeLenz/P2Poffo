module P2Poffo

    using LinearAlgebra
    using SparseArrays
    using StaticArrays
    using Lgmsh
    using Gmsh

    # Carrega as rotinas
    include("conversor.jl")
    include("elemento.jl")
    include("global.jl")
    include("apoios.jl")
    include("corpo.jl")
    include("Pre.jl" )
    include("tensoes.jl")
    include("rotacao.jl")

    # Exporta o arquivo .dat
    include("exportaarquivo.jl")
    
	# Exporta a rotina principal de análise
	export Pre_processamento

    # Exporta a rotina que cria o .dat
    export ExportaDat

	export muda_coordenada
end
