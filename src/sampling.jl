


function tt_sampling(mps::AbstractMPS)
    L = length(mps)
    hstorage = Vector{Any}(undef, L+1)
    hstorage[L+1] = [1.]
    for i in L:-1:1
        hstorage[i] = update_tt_right(hstorage[i+1], mps[i])
    end
    isapprox(sum(hstorage[1]), 1.) || error("total probability is not 1.")
    r = Int[]
    for i in 1:L
        phi_x = contract(contract(hstorage[i], mps[i], ((1,), (1,))), hstorage[i+1], ((2,), (1,)))
        (maximum(abs.(imag.(phi_x))) <= 1.0e-8) || error("imaginary part of probability distribution is too large.")
        phi_x = real.(phi_x)
        isapprox(sum(phi_x), 1.) || error("total probability is not 1.")
        rj = discrete_sample(phi_x)
        push!(r, rj)
        hstorage[i+1] = contract(hstorage[i], mps[i][:, rj, :] / phi_x[rj] , ((1,), (1,)))
    end
    return r
end

function tt_probability(mps::AbstractMPS, v::Vector{Int})
    (length(mps) == length(v)) || error("tt mismatch number of basis.")
    return dot(prodmps(physical_dimensions(mps), v .- 1), mps)
end

tt_probability(mps::AbstractMPS, v::String) = tt_probability(mps, [parse(Int, item) for item in v])


function discrete_sample(l::Vector{Float64})
    isempty(l) && error("no results.")
    s = sum(l)
    L = length(l)
    l1 = Vector{Float64}(undef, L+1)
    l1[1] = 0
    for i=1:L
        l1[i+1] = l1[i] + l[i]/s
    end
    s = rand(Float64)
    for i = 1:L
        if (s >= l1[i] && s < l1[i+1])
            return i
        end
    end
    return L
end