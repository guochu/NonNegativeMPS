module NonNegativeMPS


export renormalize_tt!, createrandomtt, tt_sum, tt_sampling, cross_approximation
export tt_probability, dm_to_prop, prop_to_dm, marginal_distribution, ntt_mu_add
export dosweep!, dosweeps!, compute!


using LinearAlgebra, TensorOperations

include("auxiliary/truncation.jl")
include("auxiliary/tensorops.jl")
include("auxiliary/distance.jl")

include("states/states.jl")

include("constants.jl")
include("transfer.jl")
include("constructor.jl")
include("sampling.jl")
include("conversion.jl")
include("marginal.jl")

include("iterativealgs/iterativealgs.jl")



end