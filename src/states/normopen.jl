function norm_open(mps::AbstractMPS, idens::Vector)
	isempty(mps) && error("mps must not be empty.")
	L = length(mps)
	hold = ones(scalartype(mps), 1)
	for i = 1:L
		hold = updateCleftOpenh1h2(hold, nothing, mps[i], idens[i])
	end
	return only(hold)
end

_iden_state(d::Int) = reshape(one(zeros(d, d)), d*d)
_get_iden_state(d::Int) = _iden_state(round(Int, sqrt(d)))
norm_open(mps::AbstractMPS) = norm_open(mps, _get_iden_state.([shape(item, 2) for item in mps]))

function updateCleftOpenh1h2(hold::AbstractVector, obj, mpsj::AbstractArray{<:Number, 3}, iden::AbstractVector)
	Hnew = contract(hold, mpsj, ((1,), (1,)))
	if obj !== nothing
		Hnew = contract(obj, Hnew, ((2,), (1,)))
	end
	Hnew = contract(iden, Hnew, ((1,), (1,)))
	return Hnew
end