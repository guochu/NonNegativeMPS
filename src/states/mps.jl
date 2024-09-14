
struct MPS{T<:Number} <: AbstractMPS
	data::Vector{Array{T, 3}}

function MPS{T}(data::Vector) where {T<:Number}
	_check_mps_space(data)
	new{T}(convert(Vector{Array{T, 3}}, data))
end
end
MPS(data::Vector{<:MPSTensor{T}}) where {T <: Number} = MPS{T}(data)
scalartype(::Type{MPS{T}}) where {T} = T

Base.copy(psi::MPS) = MPS(copy(raw_data(psi)))
Base.conj(psi::MPS) = MPS(conj.(raw_data(psi)))

Base.vcat(psiA::MPS, psiB::MPS) = MPS(vcat(raw_data(psiA), raw_data(psiB)))
# Base.conj(psi::MPS) = MPS(conj.(raw_data(psi)), raw_singular_matrices(psi))

Base.complex(psi::MPS) = MPS(complex.(raw_data(psi)))
isleftcanonical(a::MPS; kwargs...) = all(x->isleftcanonical(x; kwargs...), a.data)
isrightcanonical(a::MPS; kwargs...) = all(x->isrightcanonical(x; kwargs...), a.data)

function max_bond_dimensions(physpaces::Vector{Int}, D::Int) 
	L = length(physpaces)
	left = 1
	right = 1
	virtualpaces = Vector{Int}(undef, L+1)
	virtualpaces[1] = left
	for i in 2:L
		virtualpaces[i] = min(virtualpaces[i-1] * physpaces[i-1], D)
	end
	virtualpaces[L+1] = right
	for i in L:-1:2
		virtualpaces[i] = min(virtualpaces[i], physpaces[i] * virtualpaces[i+1])
	end
	return virtualpaces
end
max_bond_dimensions(psi::MPS, D::Int) = max_bond_dimensions(physical_dimensions(psi), D)


function increase_bond!(psi::MPS; D::Int)
	if bond_dimension(psi) < D
		virtualpaces = max_bond_dimensions(physical_dimensions(psi), D)
		for i in 1:length(psi)
			sl = max(min(virtualpaces[i], D), size(psi[i], 1))
			sr = max(min(virtualpaces[i+1], D), size(psi[i], 3))
			m = zeros(eltype(psi), sl, size(psi[i], 2), sr)
			m[1:size(psi[i], 1), :, 1:size(psi[i], 3)] .= psi[i]
			psi[i] = m
		end
	end
	return psi
end


function _check_mps_space(mpstensors::Vector)
	L = length(mpstensors)
	for i in 1:L-1
		(space_r(mpstensors[i]) == space_l(mpstensors[i+1])) || throw(DimensionMismatch())
	end

	# just require the left boundary to be a single sector
	(space_l(mpstensors[1]) == 1) || throw(DimensionMismatch("left boundary should be size 1."))
	# (size(mpstensors[L], 3) == 1) || throw(DimensionMismatch("right boundary should be size 1."))
	return true
end
