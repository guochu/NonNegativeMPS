

function update_tt_right(hold::AbstractArray{<:Number, 1}, mpsj::AbstractArray{<:Number, 3})
    tmp = contract(mpsj, hold, ((3,), (1,)))
    return dropdims(sum(tmp, dims=2), dims=2)
end

function update_tt_left(hold::AbstractArray{<:Number, 1}, mpsj::AbstractArray{<:Number, 3})
    tmp = contract(hold, mpsj, ((1,), (1,)))
    return dropdims(sum(tmp, dims=1), dims=1)
end
