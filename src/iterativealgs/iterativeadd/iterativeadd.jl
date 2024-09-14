
include("ntt_mu.jl")



function _convert_single_data(data::String, ds::Vector{Int})
    v = [item-1 for item in parse.(Int, [a for a in data])]
    return prodmps(ds, v)
end

function frequency_counting(data::Vector{String})
    data_r = Dict{String, Int}()
    for item in data
        v = get(data_r, item, nothing)
        if v === nothing
            data_r[item] = 1
        else
            data_r[item] += 1
        end
    end
    return data_r
end

function cross_approximation(data::Dict{String, Int}, ds::Vector{Int}; kwargs...)
    L = length(ds)
    N = 0
    for (item, v) in data
        (length(item)==L) || error("data size mismatch with number of sites.")
        N += v
    end
    mpsxs = [_convert_single_data(k, ds) * (v/N) for (k, v) in data]
    y = ntt_mu_add(mpsxs; kwargs...)
    return renormalize_tt!(y)
end

cross_approximation(data::Vector{String}, ds::Vector{Int}; kwargs...) = cross_approximation(
frequency_counting(data), ds; kwargs...)


cross_approximation(data::Vector{Vector{Int}}, ds::Vector{Int}; kwargs...) = cross_approximation(
[join(string.(item)) for item in data], ds; kwargs...)
