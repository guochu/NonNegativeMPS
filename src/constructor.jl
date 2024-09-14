

function _compute_tt_D(d::Vector{Int}, Ds::Vector{Int})
	(length(d)==1) && return [1, 1]
	L = length(d)
	Ds_n = [1 for i =1:(L+1)]
	Lhalf = div(L, 2)
	for i in 1:Lhalf
		Ds_n[i+1] = min(Ds[i], Ds_n[i]*d[i])
	end
	for i in L:-1:(Lhalf+1)
		Ds_n[i] = min(Ds[i-1], Ds_n[i+1]*d[i])
	end
	return Ds_n
end

function renormalize_tt!(x::AbstractMPS)
    hold = [1.]
    for i = length(x):-1:1
        hold = update_tt_right(hold, x[i])
        if (i%3 == 0)
            s = sum(hold)
            x[i] /= s
            hold /= s
        end
    end
    sc = sum(hold)
    x[1] /= sc
    return x
end

function createrandomtt(ds::Vector{Int}, bonds::Vector{Int})
    (length(ds) == length(bonds)+1) || error("number of physical dimensions mismatch with bond dimensions.")
    L = length(ds)
    Ds = _compute_tt_D(ds, bonds)
    mps = Vector{Array{Float64, 3}}(undef, L)
    for i in 1:L
        mps[i] = rand(Float64, Ds[i], ds[i], Ds[i+1])
        mps[i] ./= sqrt(length(mps[i]))
    end
    return renormalize_tt!(MPS(mps))
end

createrandomtt(ds::Vector{Int}; D::Int) = createrandomtt(ds, [D for i in 1:length(ds)-1])
createrandomtt(L::Int; d::Int, D::Int) = createrandomtt([d for i in 1:L], D=D)

function tt_sum(x::AbstractMPS)
    hold = [1.]
    for i in length(x):-1:1
        hold = update_tt_right(hold, x[i])
    end
    return sum(hold)
end
