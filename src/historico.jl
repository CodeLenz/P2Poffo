mutable struct HistoricoOtim
    iters     ::Vector{Int}
    freq      ::Matrix{Float64}   # niter × n_modos
    volume    ::Vector{Float64}
    FS        ::Vector{Vector{Float64}}
    densidades::Matrix{Float64}
end

function HistoricoOtim(n_modos::Int, ne::Int, nFS::Int)
    HistoricoOtim(
        Int[],
        Matrix{Float64}(undef, 0, n_modos),
        Float64[],
        [Float64[] for _ in 1:nFS],
        Matrix{Float64}(undef, 0, ne),
    )
end