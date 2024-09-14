
abstract type AbstractIterativeTT end

get_D(s::AbstractIterativeTT) = s.paras[:maxbonddimension]
get_bonddimension(s::AbstractIterativeTT) = get_D(s)
set_bonddimension!(s::AbstractIterativeTT, D::Int) = s.paras[:maxbonddimension] = D

get_kvals(x::AbstractIterativeTT) = x.kvals
get_bonds(x::AbstractIterativeTT) = x.bonds
get_errors(x::AbstractIterativeTT) = x.errors
get_mps(x::AbstractIterativeTT) = x.mps

function onesitesweep!(m::AbstractIterativeTT; verbose::Int=1, kwargs...)
	Evals1 = onesitesweepleft!(m; verbose=verbose, kwargs...)
	Evals2 = onesitesweepright!(m; verbose=verbose, kwargs...)
	r = vcat(Evals1, Evals2)
	append!(get_kvals(m), r)
	(verbose >= 2) && println("energy change from $(r[1]) to $(r[end]) after onesite sweep.")
	return r
end

function dosweep!(m::AbstractIterativeTT; single_site::Bool=true, kwargs...)
    single_site || error("two site algorithm not implemented.")
    return onesitesweep!(m; kwargs...)
end

function dosweeps!(s::AbstractIterativeTT, n::Int=10; kwargs...)
	(n >= 1) || error("number of sweep must be larger than 1.")
	r = dosweep!(s; kwargs...)
	for i in 2:n
	    tmp = dosweep!(s; kwargs...)
	    append!(r, tmp)
	end
	return r
end

compute_tolerance(s::AbstractIterativeTT, v::Vector) = iterative_error_2(v)
compute!(dmrg::AbstractIterativeTT; kwargs...) = dmrg_compute!(dmrg; kwargs...)


iterative_error_2(m::AbstractVector) = std(m) / abs(mean(m))

function dmrg_compute!(dmrg; maxitr::Int=DEFAULT_DMRG_MAXITER,
	tol::Real=DEFAULT_DMRG_CONVERGE_TOL, verbose::Int=1, kwargs...)
	err = 1.
	i = 0
	converged = false
	while (i < maxitr)
		(verbose >= 2) && println("we are at the $i-th sweep...")
		Evals = dosweep!(dmrg; verbose=verbose, kwargs...)
		i += 1
		err = compute_tolerance(dmrg, Evals)
		if (err < tol)
			converged = true
			(verbose > 2) && println("DMRG converges in $i sweeps")
			break
		end
	end
	if !converged
		i = 0
		(verbose >= 1) && @warn "DMRG fails to converge in $maxitr sweeps"
	end
	return (i, err)
end