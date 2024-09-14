

function _marginal_distribution_impl(mps::AbstractMPS, site::Int)
	isempty(mps) && error("mps must not be empty.")
	L = length(mps)
	(site <= L) || error("site out of bound.")
	# m = _tracej(mps[site], 2, diagonals)
	m = dropdims(sum(mps[site], dims=2), dims=2)
	(L==1) && return only(m)
	mpsout = typeof(mps)(L-1)
	if site==1
		mpsj = contract(m, mps[2], ((2,), (1,)))
		mpsout[site] = mpsj
		for i = (site+1):(L-1)
			mpsout[i] = mps[i+1]
		end
	else
		mpsj = contract(mps[site-1], m, ((3,), (1,)))
		for i = 1:(site-2)
			mpsout[i] = mps[i]
		end
		mpsout[site-1] = mpsj
		for i = site:(L-1)
			mpsout[i] = mps[i+1]
		end
	end
	return mpsout
end

"""
    marginal_distribution(mps::AbstractMPS, sites::Vector{Int})
compute the marginal distribution on sites
"""
function marginal_distribution(mps::AbstractMPS; trace_out::Vector{Int})
	sites = sort(trace_out, rev=true)
	for i = 1:length(sites)
		mps = _marginal_distribution_impl(mps, sites[i])
	end
	return mps
end
