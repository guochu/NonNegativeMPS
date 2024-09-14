

function _single_proj(v)
    m = kron(v, v')/2
    Iv = one(zeros(size(m)))
    return reshape(Iv, length(Iv))' * kron(Iv, m)
end

"""
    The second dimension is the physical dimension.
"""
function create_povm()
    v1 = [1., 0]
    v2 = [1/sqrt(3), sqrt(2/3)]
    v3 = [1/sqrt(3), sqrt(2/3)*exp( (2*pi*im)/3 )]
    v4 = [1/sqrt(3), sqrt(2/3)*exp( (4*pi*im)/3 )]
    vs = _single_proj.([v1, v2, v3, v4])
    return vcat(vs...)
end


function dm_to_prop(mps::AbstractMPS)
    all(physical_dimensions(mps) .== 4) || error("only spin half dm are supported")
    isapprox(norm_open(mps), 1.0, atol=1.0e-6) || error("trace norm of dm is not 1")
    povm = create_povm()
    return MPS([permute(contract(povm, item, ((2,), (2,))), (2,1,3)) for item in mps])
end

function prop_to_dm(mps::AbstractMPS)
    all(physical_dimensions(mps) .== 4) || error("only spin half prop are supported")
    isapprox(tt_sum(mps), 1.0, atol=1.0e-6) || error("tt sum is not 1")
    povm = inv(create_povm())
    return MPS([permute(contract(povm, item, ((2,), (2,))), (2,1,3)) for item in mps])
end
