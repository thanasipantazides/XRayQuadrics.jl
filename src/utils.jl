import HDF5

function get_reflection_data(file::String)
    fid = HDF5.h5open(file, "r")
    angle = HDF5.read(fid["angle"])
    energy = HDF5.read(fid["energy"])
    reflectivity = HDF5.read(fid["reflectivity/ir"])
    close(fid)

    keydata = cat(angle, reflectivity, dims=2)
    keydata = cat([0.0 energy'], keydata, dims=1)

    return keydata
end

struct Parabola{T}
    a::T
    c::T
end

struct Hyperbola{T}
    a::T
    b::T
end

@doc raw"""
    function fitp(x::Vector{<:Real}, y::Vector{<:Real})

Fit a parabola to data x and y data. Must provide 2 or more points in each argument. Will result in the form

``y = ax^2 + c``

Coefficients ``a``, ``c`` will be returned in a `Parabola` `struct`.
"""
function fitp(x::Vector{<:Real}, y::Vector{<:Real})
    n = min(length(x), length(y))
    if n < 2
        error("Need at least 2 points in each argument!")
    end

    A = zeros(n, 2)
    b = zeros(n)
    for i in 1:n
        A[i, :] = [x[i]^2 1]
        b[i] = y[i]
    end
    res = A\b
    out = Parabola(res[1], res[2])
    return out
end

@doc raw"""
    function fith(x::Vector{<:Real}, y::Vector{<:Real})

Fit a hyperbola to data x and y data. Must provide 2 or more points in each argument. Will result in the form

``1 = \frac{(x - \sqrt{a^2 + b^2})^2}{a^2} - \frac{y^2}{b^2}``

Coefficients ``a``, ``b`` will be returned in a `Hyperbola` `struct`.
"""
function fith(x::Vector{<:Real}, y::Vector{<:Real})
    n = min(length(x), length(y))
    if n < 2
        error("Need at least 2 points in each argument!")
    end

    A = zeros(n, 2)
    b = zeros(n)
    for i in 1:n
        A[i, :] = [x[i]^2 + 1 x]
        b[i] = y[i]^2
    end
    res = A\b

    out = Hyperbola(res[1], res[2], res[3])
    return out
end

function focus(object::Union{Parabola, Hyperbola})
    if typeof(object) === Parabola
        return object.c + 1/4/object.a
    end
    if typeof(object) === Hyperbola

    end
end