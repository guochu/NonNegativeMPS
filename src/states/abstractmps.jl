const ValidIndices = Union{Integer,AbstractRange{Int64}, Colon}

const MPSTensor{T} = AbstractArray{T, 3} where {T<:Number}
abstract type AbstractMPS end


"""
	The raw mps data as a list of 3-dimension tensors
	This is not supposed to be directly used by users
"""
raw_data(psi::AbstractMPS) = psi.data

scalartype(psi::AbstractMPS) = scalartype(typeof(psi))
Base.length(psi::AbstractMPS) = length(raw_data(psi))
Base.isempty(psi::AbstractMPS) = isempty(raw_data(psi))
Base.size(psi::AbstractMPS, i...) = size(raw_data(psi), i...)
Base.getindex(psi::AbstractMPS, i::ValidIndices) = getindex(raw_data(psi), i)
Base.lastindex(psi::AbstractMPS) = lastindex(raw_data(psi))
Base.firstindex(psi::AbstractMPS) = firstindex(raw_data(psi))

function Base.setindex!(psi::AbstractMPS, v, i::ValidIndices)
	return setindex!(raw_data(psi), v, i)
end 

space_l(m::MPSTensor) = size(m, 1)
space_r(m::MPSTensor) = size(m, 3)
space_l(psi::AbstractMPS) = size(psi[1], 1)
space_r(psi::AbstractMPS) = size(psi[end], 3)

r_RR(psiA::AbstractMPS, psiB::AbstractMPS) = _eye(promote_type(eltype(psiA), eltype(psiB)), space_r(psiA), space_r(psiB))
r_RR(psi::AbstractMPS) = r_RR(psi, psi)
l_LL(psiA::AbstractMPS, psiB::AbstractMPS) = _eye(promote_type(eltype(psiA), eltype(psiB)), space_l(psiA), space_l(psiB))
l_LL(psi::AbstractMPS) = l_LL(psi, psi)

"""
	bond_dimension(psi::MPS[, bond::Int])
	bond_dimension(h::MPO[, bond::Int])
return bond dimension at the given bond, or return the largest bond dimension of all bonds.
"""
bond_dimension(psi::AbstractMPS, bond::Int) = begin
	((bond >= 1) && (bond <= length(psi))) || throw(BoundsError())
	space_r(psi[bond])
end 
bond_dimensions(psi::AbstractMPS) = [bond_dimension(psi, i) for i in 1:length(psi)]
bond_dimension(psi::AbstractMPS) = maximum(bond_dimensions(psi))

"""
	physical_dimensions(psi::MPS)
	physical_dimensions(psi::MPO) 
Return all the physical spaces of MPS or MPO
"""
physical_dimensions(psi::AbstractMPS) = [size(item, 2) for item in raw_data(psi)]



function isleftcanonical(psij::MPSTensor; kwargs...)
	@tensor r[-1, -2] := conj(psij[1,2,-1]) * psij[1,2,-2]
	return isapprox(r, one(r); kwargs...) 
end
function isrightcanonical(psij::MPSTensor; kwargs...)
	@tensor r[-1, -2] := conj(psij[-1,1,2]) * psij[-2,1,2]
	return isapprox(r, one(r); kwargs...) 
end