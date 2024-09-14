
abstract type AbstractTTIterativeAdd <: AbstractIterativeTT end

get_mpsxs(x::AbstractTTIterativeAdd) = x.mpsxs
get_hstorage(x::AbstractTTIterativeAdd) = x.hstorage
get_qstorage(x::AbstractTTIterativeAdd) = x.qstorage

struct TTIterativeAdd_MU{M} <: AbstractTTIterativeAdd
    mps::MPS{Float64}
    mpsxs::Vector{M}
    hstorage::Vector{Vector{Array{Float64, 2}}}
    qstorage::Vector{Array{Float64, 2}}
    paras::Dict{Symbol, Any}

	# results
	kvals::Vector{Float64}
	bonds::Vector{Int}
	errors::Vector{Float64}
end

function TTIterativeAdd_MU(mpsxs::Vector{<:AbstractMPS}; maxbonddimension::Int=DEFAULT_DMRG_MAX_BOND_DIMENSION,
	svdcutoff::Float64=DEFAULT_DMRG_SVD_TRUNCATION)
    isempty(mpsxs) && error("no data.")
    mps = createrandomtt(physical_dimensions(mpsxs[1]), D=maxbonddimension)
    canonicalize!(mps, alg=Orthogonalize(normalize=true))
    for i in 1:length(mps)
        mps[i] = abs.(mps[i])
    end
    qstorage = init_cstorage_right(mps, mps)
    hstorage = [init_cstorage_right(item, mps) for item in mpsxs]
    kvals = Float64[]
    bonds = Int[]
    errors = Float64[]
    paras = Dict{Symbol, Any}(:maxbonddimension=>maxbonddimension, :svdcutoff=>svdcutoff)
    return TTIterativeAdd_MU(mps, mpsxs, hstorage, qstorage, paras, kvals, bonds, errors)
end


function onesitesweepleft!(m::TTIterativeAdd_MU; verbose::Int=1)
    Evals = Float64[]
    mps = get_mps(m)
    mpsxs = get_mpsxs(m)
    hstorage = get_hstorage(m)
    qstorage = get_qstorage(m)
    # N = length(mpsxs)
    for site in 1:length(mps)-1
        (verbose > 2) && println("we are sweeping from left to right at site: $site")
        heff = reduceD_single_site(mps[site], qstorage[site], qstorage[site+1])
        xeff = zero(heff)
        for (xj, hj) in zip(mpsxs, hstorage)
            xeff .+= reduceD_single_site(xj[site], hj[site], hj[site+1])
        end
        mps[site] .*= (xeff ./ heff)
        # println("two norm is $(two_norm(mps[site])).")
        err = dot(reduceD_single_site(mps[site], qstorage[site], qstorage[site+1]), mps[site]) - 2*dot(xeff, mps[site])
        push!(Evals, err)
        (verbose > 2) && println("residual is $(Evals[end])...")
        qstorage[site+1] = updateleft(qstorage[site], mps[site], mps[site])
        for n in 1:length(hstorage)
            hstorage[n][site+1] = updateleft(hstorage[n][site], mpsxs[n][site], mps[site])
        end
    end

    return Evals
end

function onesitesweepright!(m::TTIterativeAdd_MU; verbose::Int=1)
    Evals = Float64[]
    mps = get_mps(m)
    mpsxs = get_mpsxs(m)
    hstorage = get_hstorage(m)
    qstorage = get_qstorage(m)
    # N = length(mpsxs)
    for site in length(mps):-1:2
        (verbose > 2) && println("we are sweeping from right to left at site: $site")
        heff = reduceD_single_site(mps[site], qstorage[site], qstorage[site+1])
        xeff = zero(heff)
        for (xj, hj) in zip(mpsxs, hstorage)
            xeff .+= reduceD_single_site(xj[site], hj[site], hj[site+1])
        end
        mps[site] .*= (xeff ./ heff)
        # println("two norm is $(two_norm(mps[site])).")
        err = dot(reduceD_single_site(mps[site], qstorage[site], qstorage[site+1]), mps[site]) - 2*dot(xeff, mps[site])
        push!(Evals, err)
        (verbose > 2) && println("residual is $(Evals[end])...")
        qstorage[site] = updateright(qstorage[site+1], mps[site], mps[site])
        for n in 1:length(hstorage)
            hstorage[n][site] = updateright(hstorage[n][site+1], mpsxs[n][site], mps[site])
        end
    end

    # qstorage .= init_cstorage_right(mps, mps)
    # hstorage .= [init_cstorage_right(item, mps) for item in mpsxs]
    return Evals
end

function ntt_mu_add(mpsxs::Vector{<:AbstractMPS}; maxbonddimension::Int=DEFAULT_DMRG_MAX_BOND_DIMENSION,
    svdcutoff::Float64=DEFAULT_DMRG_SVD_TRUNCATION,
    maxitr::Int=DEFAULT_DMRG_MAXITER, tol::Float64=DEFAULT_DMRG_CONVERGE_TOL,
    single_site::Bool=true, verbose::Int=1)
    dmrg = TTIterativeAdd_MU(mpsxs, maxbonddimension=maxbonddimension, svdcutoff=svdcutoff)
    i, err = compute!(dmrg, maxitr=maxitr, tol=tol, single_site=single_site, verbose=verbose)
    return get_mps(dmrg)
end

function reduceD_single_site(A::MPSTensor, Cleft::AbstractMatrix, Cright::AbstractMatrix)
    @tensor tmp[1 3; 5] := Cleft[1,2] * A[2,3,4] * Cright[5,4]
    return tmp
end


function init_cstorage_right(a::AbstractMPS, b::AbstractMPS) 
    L = length(a)
    T = promote_type(scalartype(a), scalartype(b))
    Cstorage = Vector{Matrix{T}}(undef, L+1)
    Cstorage[L+1] = ones(T, 1, 1)
    for i = L:-1:1
        Cstorage[i] = updateright(Cstorage[i+1], a[i], b[i])
    end     
    return Cstorage
end