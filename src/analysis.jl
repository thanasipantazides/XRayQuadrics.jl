"""
    bin(data<:Vector{<AbstractFloat}, ticks:Vector{AbstractFloat})

Makes low-to-high histogram of `data` split into bins demarked by `ticks` spanning the range of `data`. Returns the weight of each bin (counts of data per bin), a list of data per bin, and the indices into input `data` set for each bin.
"""
function bin(data::Vector{<:AbstractFloat}, ticks::Vector{<:AbstractFloat})
    sortI = sortperm(data)
    sort!(ticks)
    n = length(ticks)

    weights = zeros(n-1)
    binvals = fill([], n-1)
    bindices = fill([], n-1)
    for i = 1:n-1
        inbin = findall(x->x>ticks[i] && x<ticks[i+1], data)
        weights[i] = length(inbin)
        bindices[i] = inbin
        binvals[i] = data[inbin]
    end

    return (weights, binvals, bindices)
end

bin(data::Vector{<:AbstractFloat}, ticks::LinRange) = bin(data, Vector(ticks))
bin(data::Vector{<:AbstractFloat}, n::Integer) = bin(data, LinRange(min(data...), max(data...), n))

