module P2Poffo

    using LinearAlgebra
    using SparseArrays
    using StaticArrays
    using Statistics
    using Lgmsh
    using Gmsh

    # Carrega as rotinas
    include("conversor.jl")
    include("elemento.jl")
    include("global.jl")
    include("apoios.jl")
    include("corpo.jl")
    include("gradiente.jl")
    include("rotacao.jl")

    # Rotinas "principais"
    include("Pre.jl" )
    include("Pos.jl" )

    # Exporta o arquivo .dat
    include("exportaarquivo.jl")
    
	# Exporta a rotina principal de análise
	export Pre_processamento, Pos_processamento, ExportaDat 

end
