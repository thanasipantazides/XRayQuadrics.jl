import HDF5, CSV

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

@doc raw"""
    function get_photon_data(file::String)
        
Return a Matrix 
"""
function get_photon_data(file::String)
    data = CSV.File(file, header=true)
    
    energies = []
    angles = []
    for row in data
        push!(energies, row[3])
        push!(angles, row[2])
    end
    
    return Dict("energies" => energies, "angles" => angles)
end

function get_mass_attenuation_data(file::String)
    data = CSV.File(file, header=true)
    
    energies = []
    attenuations = []
    for row in data
        push!(energies, row[1])
        push!(attenuations, row[2])
    end
    
    return Dict("energies" => energies, "attenuations" => attenuations)
end

function get_empirical_attenuation_data(file::String)
    data = CSV.File(file, header=true)
    energies = []
    modeled_attenuations = []
    measured_attenuations = []
    for row in data
        push!(energies, row[1])
        push!(modeled_attenuations, row[3])
        push!(measured_attenuations, row[2])
    end
    return Dict("energies" => energies, "measured_attenuations" => measured_attenuations, "modeled_attenuations" => modeled_attenuations)
end

function parse_mca(file::String)
    
    content = read(file, String)
    
    reg_cal = r"<<CALIBRATION>>(?s)(.*)<<ROI>>"
    reg_data = r"<<DATA>>(?s)(.*)<<END>>"
    reg_time = r"(?<=REAL_TIME - )(.*)"
    
    match_cal = match(reg_cal, content)
    match_data = match(reg_data, content)
    match_time = match(reg_time, content)
    
    livetime = parse(Float64, match_time.match)
    
    cal = split(match_cal.match, '\n')
    data = split(match_data.match, '\n')
    caladc = Float64[]
    calenergy = Float64[]
    for line in cal
        if occursin("<<", line) || '-' in line
            continue
        end
        val = split(line, ' ')
        push!(caladc, parse(Float64, val[1]))
        push!(calenergy, parse(Float64, val[2]))
    end
    
    adc = Float64[]
    for line in data
        if occursin("<<", line)
            continue
        end
        push!(adc, parse(Float64, line))
    end
    
    n = length(caladc)
    A = [ones(n) caladc]
    b = calenergy
    reg = inv(A'*A)*A'*b
    println(reg)
    
    ret = Dict("fit_const"=>reg[1], "fit_slope"=>reg[2], "adc"=>adc, "livetime"=>livetime)
    return ret
end

function interpolateattenuation(E::T, table::Dict{String, Vector{Any}}; level=0) where T<:Real
    if level == 0 && (E >= table["energies"][end] || E < table["energies"][1])
        @error "looking up energy outside bounds!"
        println(E)
        return
    end
    
    lower_k = findlast(x -> x <= E, table["energies"])
    upper_k = lower_k + 1
    
    cE = table["attenuations"][lower_k] + (table["attenuations"][upper_k] - table["attenuations"][lower_k])*(E - table["energies"][lower_k]) / (table["energies"][upper_k] - table["energies"][lower_k])
    
    return cE
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

function LinearAlgebra.cross(x::AbstractVector)
    if length(x) > 3
        throw(Exception("too long!"))
    end

    return [0 -x[3] x[2];
        x[3] 0 -x[1];
        -x[2] x[1] 0]
end

function randr()::Matrix{<:Real}
    ax = rand(3) .- 0.5
    ax = ax / norm(ax)
    cang = 2*rand() - 1

    return I * cang + (1 - cang) * ax * ax' + cross(ax) * sqrt(1 - cang^2)
end

function r_euler3(ang::Real)
    return [cos(ang) -sin(ang) 0;
        sin(ang) cos(ang) 0;
        0 0 1]
end
function r_euler2(ang::Real)
    return [cos(ang) 0 sin(ang);
        -sin(ang) 0 cos(ang);
        0 1 0]
end
function r_euler1(ang::Real)
    return [1 0 0;
        0 cos(ang) -sin(ang);
        0 sin(ang) cos(ang)]
end
function r_min_arc(x_A::AbstractVector{<:Real}, x_B::AbstractVector{<:Real})
    ax = LinearAlgebra.cross(x_A, x_B)
    ax = ax / norm(ax)
    cang = x_A' * x_B / norm(x_A) / norm(x_B)

    return I * cang + (1 - cang) * ax * ax' + cross(ax) * sqrt(1 - cang^2)
end