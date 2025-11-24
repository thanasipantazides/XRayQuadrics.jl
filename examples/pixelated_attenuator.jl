using XRayQuadrics
using LinearAlgebra
using GLMakie
import CairoMakie
using Random
using ProgressMeter
import Printf, CSV
import Profile, PProf

function merge_csv(loadfiles::Vector{String})
    n = 0
    offset = 0
    savepathprefix = "" 
    for (k,f) in enumerate(loadfiles)
        n_str_in= split(f, '.')[1]
        n_str_reg = r"(?<=_n)(.*)"
        match_o = match(n_str_reg, n_str_in)
        n_str = split(match_o.match, '_')[1]
        #mantissa, exponent = split(n_str, 'e')
        n += parse(Float64, n_str)
        if k == 1
            offset = match_o.offset - 2
            savepathprefix = f[1:offset]
        end
    end
    println("merged data for " * string(n) * " photons")
    savepath = savepathprefix*"merged_n"*Printf.@sprintf("%.0e", n)*".csv"
    merge_csv(loadfiles, savepath)
    return savepath
end

function merge_csv(loadfiles::Vector{String}, savepath::String)
    
    M = zeros(2, 6, length(loadfiles))
    names = ""
    for (f,file) in enumerate(loadfiles)
        data = CSV.File(file, header=true)
        if f == 1
            names = [string(name) for name in data.names]
        end
        n = length(data)
        
        # resize M to fit the new data
        if n > length(M[:,1,f])
            P = M
            M = zeros(n, length(M[1,:,f]), length(M[1,1,:]))
            M[1:length(P[:,1,f]), :, :] = P
        end
        
        # populate the table M with each source CSV's data. Just stacking the CSVs.
        for (r,row) in enumerate(data)
            for (c,col) in enumerate(row)
                M[r,c,f] = col
            end
        end
    end
    
    N = zeros(length(M[:,1,1]), length(M[1,:,1]))
    # same energy bins as source data
    N[:,1] = M[:,1,1]
    N[:,2] = M[:,2,1]
    for k in eachindex(M[:,1,1])
        # total count in each energy bin
        N[k,3] = sum(M[k,3,:]) 
        N[k,4] = sum(M[k,4,:])
        
        # average transmission in each bin
        N[k,5] = sum(dot(M[k,5,:], M[k,3,:])) / N[k,3]
        N[k,6] = sum(dot(M[k,6,:], M[k,4,:])) / N[k,4]
    end
    
    println("\nprephotons:\t", sum(N[:,3]))
    println("postphotons:\t", sum(N[:,4]))
    
    csvstr = join(names, ",")
    # csvstr *= join(N', ",")
    for k in eachindex(N[:,1])
        csvstr *= "\n" * join(N[k,:], ",")
    end
    write(savepath, csvstr)
end

function do_merge()
    merge_csv(
        [
            "/home/thanasi/projects/ray/results/2025/nov12/f32/transmit_data_n1e+07_r1.csv",
            "/home/thanasi/projects/ray/results/2025/nov12/f32/transmit_data_n1e+07_r2.csv",
            "/home/thanasi/projects/ray/results/2025/nov12/f32/transmit_data_n1e+07_r3.csv",
            "/home/thanasi/projects/ray/results/2025/nov12/f32/transmit_data_n1e+07_r4.csv",
        ],
        "/home/thanasi/projects/ray/results/2025/nov12/f32/foxsi4_pix_atten_merge_sim_transmission.csv"
    )
end

function compare_matten(savepath::String)
    matten_nist_f = "src/data/ehsteve_matten.csv"
    matten_lbnl_f = "src/data/LBNL_attenlength_Si.csv"
    matten_nist_phot_f = "src/data/NIST_XCOM_attenlength_si.csv"
    matten_nist_scat_f = "src/data/NIST_XCOM_attenlength_scatter_si.csv"
    matten_nist_t = get_mass_attenuation_data(matten_nist_f)
    matten_lbnl_t = get_mass_attenuation_data(matten_lbnl_f)
    matten_nist_phot_t = get_mass_attenuation_data(matten_nist_phot_f)
    matten_nist_scat_t = get_mass_attenuation_data(matten_nist_scat_f)
    
    energy = (1:0.1:28) .* 1e3
    energy_range = min(energy...)..max(energy...)
    matten_nist = [interpolateattenuation(E, matten_nist_t) for E in energy]
    matten_lbnl = [interpolateattenuation(E, matten_lbnl_t) for E in energy]
    matten_nist_phot = [interpolateattenuation(E, matten_nist_phot_t) for E in energy]
    matten_nist_scat = [interpolateattenuation(E, matten_nist_scat_t) for E in energy]
    
    CairoMakie.activate!()
    fig = CairoMakie.Figure(size=(600,400))
    # ax_raw = Axis(fig[1,1], xlabel=L"\text{Energy [eV]}", ylabel=L"\text{Attenuation length [m}^-1]")
    ax_compare = Axis(fig[2:3,1], xlabel=L"\text{Energy [eV]}", ylabel=L"\frac{\mu_1(E)}{\mu_2(E)} - 1", yminorticks=IntervalsBetween(5), yminorticksvisible=true)
    
    raw_nist_overlap_range = findall(in(energy_range),matten_nist_t["energies"])
    raw_lbnl_overlap_range = findall(in(energy_range),matten_lbnl_t["energies"])
    raw_nist_phot_overlap_range = findall(in(energy_range),matten_nist_phot_t["energies"])
    raw_nist_scat_overlap_range = findall(in(energy_range),matten_nist_scat_t["energies"])
    
    # scatter!(ax_raw, matten_nist_t["energies"][raw_nist_overlap_range], matten_nist_t["attenuations"][raw_nist_overlap_range], marker='o', markersize=10, color=:black)
    # scatter!(ax_raw, matten_lbnl_t["energies"][raw_lbnl_overlap_range], matten_lbnl_t["attenuations"][raw_lbnl_overlap_range], marker='+', markersize=10, color=:red)
    
    # lines!(ax_raw, energy, matten_nist, color=:black, linewidth=1)
    # lines!(ax_raw, energy, matten_lbnl, color=:red, linewidth=1)
    # lines!(ax_raw, energy, matten_nist_phot, color=:blue, linewidth=1)
    # lines!(ax_raw, energy, matten_nist_scat, color=:green, linewidth=1)
    
    # lines!(ax_compare, energy, matten_nist ./ matten_lbnl .- 1, linewidth=1, color=:black, label="roentgen vs. LBNL")
    
    lines!(ax_compare, energy, matten_nist ./ matten_nist_scat .- 1, linewidth=1, label="Roentgen vs. NIST total")
    lines!(ax_compare, energy, matten_nist ./ matten_lbnl .- 1, linewidth=1, label="Roentgen vs. LBNL")
    lines!(ax_compare, energy, matten_nist_scat ./ matten_lbnl .- 1, linewidth=1, label="NIST total vs. LBNL")
    leg = Legend(fig[4,1], ax_compare)
    ylims!(ax_compare, [-0.25,0.25])
    
    # scatter!(ax_raw, matten_lbnl_t["energies"][raw_lbnl_overlap_range], matten_nist2_t["attenuations"][raw_lbnl_overlap_range], marker='+', markersize=5, color=:red)
    
    save(savepath, fig)
end

function make_monolithic()
    top_r = 10.5e-6
    bottom_r = 5e-6
    # top_h = 120e-6
    top_h = 104.63673925931994e-6
    bottom_h = 175e-6
    gap = 105e-6
    pitch = 60e-6
    nx = 0
    ny = 0
    
    corner1 = [0 - pitch/2, 0 - pitch/2, 0]
    corner2 = [1*pitch + pitch/2, 1*pitch + pitch/2, bottom_h + gap + top_h]
    bbox = (corner1, corner2)
    
    z = [0;0;1]
    tang = 0*pi/100000000
    transform = [1 0 0; 0 cos(tang) -sin(tang); 0 sin(tang) cos(tang)]
    bottom_plane = Plane(zeros(3), z)
    bottom_cyl_plane = Plane((bottom_h)*z, z)
    top_cyl_plane = Plane((bottom_h+gap)*z, z)
    top_plane = Plane((bottom_h+gap+top_h)*z, z)
    truncated_cylinders = TruncatedQuadric[]
    for i in 1:nx
        for j in 1:ny
            this_bottom_c = Cylinder(
                bottom_r,
                [(i-1)*pitch, (j-1)*pitch, 0],
                transform*z
            )
            this_top_c = Cylinder(
                top_r,
                [(i-1)*pitch, (j-1)*pitch, bottom_h + gap],
                transform*z
            )
            push!(truncated_cylinders, TruncatedQuadric(this_bottom_c,  [bottom_plane, bottom_cyl_plane], true))
            push!(truncated_cylinders, TruncatedQuadric(this_top_c, [top_cyl_plane, top_plane], true))
        end
    end
    
    return PixelatedAttenuator(
        truncated_cylinders,    # list of features
        max(top_r, bottom_r),   # biggest hole size
        2390.0,                 # density
        top_plane,   # top point
        bottom_plane,                      # bottom point
        z,                      # normal
        # get_mass_attenuation_data("src/data/LBNL_attenlength_Si.csv"),
        get_mass_attenuation_data("src/data/NIST_XCOM_attenlength_scatter_si.csv"),
        bbox,
        pitch
    ) 
end

function make_monolithic_al()
    top_r = 10.5e-6
    bottom_r = 5e-6
    top_h =     0.127e-3/3
    bottom_h =  0.127e-3/3
    gap =       0.127e-3/3
    pitch = 60e-6
    nx = 0
    ny = 0
    
    corner1 = [0 - pitch/2, 0 - pitch/2, 0]
    corner2 = [1*pitch + pitch/2, 1*pitch + pitch/2, bottom_h + gap + top_h]
    bbox = (corner1, corner2)
    
    z = [0;0;1]
    tang = 0*pi/100000000
    transform = [1 0 0; 0 cos(tang) -sin(tang); 0 sin(tang) cos(tang)]
    bottom_plane = Plane(zeros(3), z)
    bottom_cyl_plane = Plane((bottom_h)*z, z)
    top_cyl_plane = Plane((bottom_h+gap)*z, z)
    top_plane = Plane((bottom_h+gap+top_h)*z, z)
    truncated_cylinders = TruncatedQuadric[]
    for i in 1:nx
        for j in 1:ny
            this_bottom_c = Cylinder(
                bottom_r,
                [(i-1)*pitch, (j-1)*pitch, 0],
                transform*z
            )
            this_top_c = Cylinder(
                top_r,
                [(i-1)*pitch, (j-1)*pitch, bottom_h + gap],
                transform*z
            )
            push!(truncated_cylinders, TruncatedQuadric(this_bottom_c,  [bottom_plane, bottom_cyl_plane], true))
            push!(truncated_cylinders, TruncatedQuadric(this_top_c, [top_cyl_plane, top_plane], true))
        end
    end
    
    return PixelatedAttenuator(
        truncated_cylinders,    # list of features
        max(top_r, bottom_r),   # biggest hole size
        2699.0,                 # density
        top_plane,
        bottom_plane,
        z,                      # normal
        # get_mass_attenuation_data("src/data/LBNL_attenlength_Al.csv"),
        get_mass_attenuation_data("src/data/NIST_XCOM_attenlength_scatter_al.csv"),
        bbox,
        pitch
    ) 
end

function make_attenuator()
    # define these to instantiate PixelatedAttenuator:
    #   holes (a list of cylinders)
    #   largeradius (biggest hole radius, as a length scale)
    #   density
    #   toppoint
    #   bottompoint
    #   normal
    #   massattenuation

    top_r = 10.5e-6
    bottom_r = 5e-6
    top_h = 120e-6
    bottom_h = 175e-6
    gap = 105e-6
    pitch = 60e-6
    nx = 3
    ny = 3
    
    corner1 = [0 - pitch/2, 0 - pitch/2, 0]
    corner2 = [(nx - 1)*pitch + pitch/2, (ny - 1)*pitch + pitch/2, bottom_h + gap + top_h]
    bbox = (corner1, corner2)
    
    # The inner tracing.jl functions need to be adapted to handle TruncatedQuadrics.
    # But first just construct the attenuator and plot it for sanity check.
    
    z = [0;0;1]
    tang = 0*pi/100000000
    transform = [1 0 0; 0 cos(tang) -sin(tang); 0 sin(tang) cos(tang)]
    bottom_plane = Plane(zeros(3), z)
    bottom_cyl_plane = Plane((bottom_h)*z, z)
    top_cyl_plane = Plane((bottom_h+gap)*z, z)
    top_plane = Plane((bottom_h+gap+top_h)*z, z)
    truncated_cylinders = []
    for i in 1:nx
        for j in 1:ny
            this_bottom_c = Cylinder(
                bottom_r,
                [(i-1)*pitch, (j-1)*pitch, 0],
                transform*z
            )
            this_top_c = Cylinder(
                top_r,
                [(i-1)*pitch, (j-1)*pitch, bottom_h + gap],
                transform*z
            )
            push!(truncated_cylinders, TruncatedQuadric(this_bottom_c,  [bottom_plane, bottom_cyl_plane], true))
            push!(truncated_cylinders, TruncatedQuadric(this_top_c, [top_cyl_plane, top_plane], true))
        end
    end
    
    return PixelatedAttenuator(
        truncated_cylinders,    # list of features
        max(top_r, bottom_r),   # biggest hole size
        2390.0,                 # density
        top_plane,
        bottom_plane,
        z,                      # normal
        # get_mass_attenuation_data("src/data/LBNL_attenlength_Si.csv"),
        get_mass_attenuation_data("src/data/NIST_XCOM_attenlength_scatter_si.csv"),
        bbox,
        pitch
    ) 
end

function make_attenuator_with_hole()
    # define these to instantiate PixelatedAttenuator:
    #   holes (a list of cylinders)
    #   largeradius (biggest hole radius, as a length scale)
    #   density
    #   toppoint
    #   bottompoint
    #   normal
    #   massattenuation

    top_r = 10.5e-6
    bottom_r = 5e-6
    top_h = 120e-6
    bottom_h = 175e-6
    gap = 105e-6
    pitch = 60e-6
    
    midpitch = 2
    nx = 3
    ny = 3
    
    corner1 = [0 - pitch/2, 0 - pitch/2, 0]
    corner2 = [(nx - 1)*pitch + pitch/2, (ny - 1)*pitch + pitch/2, bottom_h + gap + top_h]
    bbox = (corner1, corner2)
    
    # The inner tracing.jl functions need to be adapted to handle TruncatedQuadrics.
    # But first just construct the attenuator and plot it for sanity check.
    
    z = [0;0;1]
    tang = 0*pi/100000000
    transform = [1 0 0; 0 cos(tang) -sin(tang); 0 sin(tang) cos(tang)]
    bottom_plane = Plane(zeros(3), z)
    bottom_cyl_plane = Plane((bottom_h)*z, z)
    top_cyl_plane = Plane((bottom_h+gap)*z, z)
    top_plane = Plane((bottom_h+gap+top_h)*z, z)
    truncated_cylinders = []
    for i in 1:nx
        for j in 1:ny
            if i == midpitch && j == midpitch # if we are the middle(-ish) of the attenuator,
                # add a cylinder straight through the whole thing.
                this_c = Cylinder(
                    max(bottom_r, top_r),
                    [(i-1)*pitch, (j-1)*pitch, 0],
                    transform*z
                )
                push!(truncated_cylinders, TruncatedQuadric(this_c,  [bottom_plane, top_plane], true))
                continue
            end
            this_bottom_c = Cylinder(
                bottom_r,
                [(i-1)*pitch, (j-1)*pitch, 0],
                transform*z
            )
            this_top_c = Cylinder(
                top_r,
                [(i-1)*pitch, (j-1)*pitch, bottom_h + gap],
                transform*z
            )
            push!(truncated_cylinders, TruncatedQuadric(this_bottom_c,  [bottom_plane, bottom_cyl_plane], true))
            push!(truncated_cylinders, TruncatedQuadric(this_top_c, [top_cyl_plane, top_plane], true))
        end
    end
    
    return PixelatedAttenuator(
        truncated_cylinders,    # list of features
        max(top_r, bottom_r),   # biggest hole size
        2390.0,                 # density
        top_plane,
        bottom_plane,
        z,                      # normal
        # get_mass_attenuation_data("src/data/LBNL_attenlength_Si.csv"),
        get_mass_attenuation_data("src/data/NIST_XCOM_attenlength_scatter_si.csv"),
        bbox,
        pitch
    ) 
end

function setup_n_photons(att::PixelatedAttenuator, path::Union{String, Nothing}, n::Int; transform::AbstractMatrix=I(3), hires=false)
    # the photons should be populated over a unit cell of the pixelated attenuator---any 60 µm x 60 µm region. 
    # Then they should be traced back along their negative velocity direction to starting points (so they don't start inside the volume).
    # This way we have full angular coverage and a uniform sampling for the attenuator.
    
    angles = zeros(n)
    if !isnothing(path)
        d = get_photon_data(path)
        angles = d["angles"]
    end
    
    if hires
        angles = angles[((angles .> 0.9) .& (angles .< 0.95)) .| ((angles .> 1.02) .& (angles .< 1.07))]
    end
    println("Source photon angles [deg]: ", sum(angles)/length(angles))
    
    # Shuffle the incident photons, so their is no correlation between their spatial arrangement and the original list order.
    # randi = randperm(length(angles))
    # angles = angles[randi]
     
    npitch = 1
    midatten = (att.bbox[1] .+ att.bbox[2]) ./ 2
    midatten[3] = att.topplane.c[3]
    
    out_energies = zeros(n)
    out_angles = zeros(n)
    photons = Vector{Particle}(undef, n)
    # out_energies = []
    # out_angles = []
    # photons = Particle[]
    @showprogress desc = "setting up photons...\t" for i in 1:n
        # energy = rand()*28e3 + 1e3 # random energy from 1 keV to 29 keV.
        energy = rand()*28e3 + 1e3 # random energy from 1 keV to 11 keV.
        z = [0,0,1]
        anga = rand(angles)*pi/180
        angc = rand()*2*pi
        
        euler3 = [cos(angc) -sin(angc) 0;
                  sin(angc) cos(angc) 0;
                  0  0         1]
        euler1 = [1  0         0;
                  0  cos(anga) -sin(anga);
                  0  sin(anga) cos(anga)]
        # v = euler3*[0, sin(angles[i]*pi/180), -cos(angles[i]*pi/180)]
        v = transform*euler3*euler1*(-z)
        v = v/norm(v)
        r = Point3f([rand(), rand(), 0].*att.pitch*npitch .+ midatten .- [npitch*att.pitch/2, npitch*att.pitch/2, 0])
        photon = Particle(
            r - 600e-6*v,
            v,
            energy,
            i
        )
        photons[i] = photon
        out_energies[i] = energy
        out_angles[i] = anga*180/pi
    end
    return photons, out_energies, out_angles
end

function setup_n_photons!(photons::Vector{Particle}, out_energies::Vector{<:Real}, out_angles::Vector{<:Real}, att::PixelatedAttenuator, path::Union{String, Nothing}, n::Int; transform::AbstractMatrix=I(3), hires=false)
    # the photons should be populated over a unit cell of the pixelated attenuator---any 60 µm x 60 µm region. 
    # Then they should be traced back along their negative velocity direction to starting points (so they don't start inside the volume).
    # This way we have full angular coverage and a uniform sampling for the attenuator.
    
    angles = zeros(n)
    if !isnothing(path)
        d = get_photon_data(path)
        angles = d["angles"]
    end
    
    if hires
        angles = angles[((angles .> 0.9) .& (angles .< 0.95)) .| ((angles .> 1.02) .& (angles .< 1.07))]
    end
    println("Source photon angles [deg]: ", sum(angles)/length(angles))
    
    # Shuffle the incident photons, so their is no correlation between their spatial arrangement and the original list order.
    # randi = randperm(length(angles))
    # angles = angles[randi]
     
    npitch = 1
    midatten = (att.bbox[1] .+ att.bbox[2]) ./ 2
    midatten[3] = att.topplane.c[3]
    
    # out_energies = []
    # out_angles = []
    # photons = Particle[]
    @showprogress desc = "setting up photons...\t" for i in 1:n
        # energy = rand()*28e3 + 1e3 # random energy from 1 keV to 29 keV.
        energy = rand()*28e3 + 1e3 # random energy from 1 keV to 11 keV.
        z = [0,0,1]
        anga = rand(angles)*pi/180
        angc = rand()*2*pi
        
        euler3 = [cos(angc) -sin(angc) 0;
                  sin(angc) cos(angc) 0;
                  0  0         1]
        euler1 = [1  0         0;
                  0  cos(anga) -sin(anga);
                  0  sin(anga) cos(anga)]
        # v = euler3*[0, sin(angles[i]*pi/180), -cos(angles[i]*pi/180)]
        v = transform*euler3*euler1*(-z)
        v = v/norm(v)
        r = Point3f([rand(), rand(), 0].*att.pitch*npitch .+ midatten .- [npitch*att.pitch/2, npitch*att.pitch/2, 0])
        photon = Particle(
            r - 600e-6*v,
            v,
            energy,
            i
        )
        photons[i] = photon
        out_energies[i] = energy
        out_angles[i] = anga*180/pi
    end
end

function plot_attenuator!(ax::GLMakie.LScene, attenuator::PixelatedAttenuator)
    for cyl in attenuator.holes
        XRayQuadrics.plot!(ax, cyl)
    end
    p1 = attenuator.bbox[1]
    p2 = attenuator.bbox[2]
    perm_border = [
        Point3f(p1),
        Point3f(p1[1], p2[2], p1[3]),
        Point3f(p2[1], p2[2], p1[3]),
        Point3f(p2[1], p1[2], p1[3]),
        Point3f(p1),
        Point3f(p1[1], p1[2], p2[3]),
        Point3f(p1[1], p2[2], p2[3]),
        Point3f(p2[1], p2[2], p2[3]),
        Point3f(p2[1], p1[2], p2[3]),
        Point3f(p1[1], p1[2], p2[3]),
        Point3f(p2[1], p1[2], p2[3]),
        Point3f(p2[1], p1[2], p1[3]),
        Point3f(p2[1], p2[2], p1[3]),
        Point3f(p2[1], p2[2], p2[3]),
        Point3f(p1[1], p2[2], p2[3]),
        Point3f(p1[1], p2[2], p1[3]),
    ]
    lines!(ax, perm_border, color=:black)
end

function save_products(path::String, data::Dict{String, Vector}; do_gsfc_meas=false, index=nothing)
    CairoMakie.activate!()
    fig = CairoMakie.Figure(size=(600, 600), backgroundcolor=:transparent)
    
    gsfc_data = get_empirical_attenuation_data("src/data/20240607_foxsi4_transmission.csv")
    gsfc_energy = gsfc_data["energies"] .* 1e3
    gsfc_mod = gsfc_data["modeled_attenuations"]
    gsfc_meas = gsfc_data["measured_attenuations"]
    
    bins = data["bins"]   # each energy bin, to put photon data into
    pretransmits = data["pretransmits"] # by each photon
    posttransmits = data["posttransmits"]
    prephotons = data["prephotons"]
    postphotons = data["postphotons"]
    
    # count photons for reference
    nphotons = length(prephotons)
    n_str = Printf.@sprintf("%.0e", nphotons)
    mantissa, exponent = split(n_str, 'e')
    n_str_latex = L"%$mantissa \times 10^{%$exponent}"
    
    preenergies = [p.E for p in prephotons]
    postenergies = [p.E for p in postphotons]
    
    prew, prev, preI = bin(preenergies, bins)
    postw, postv, postI = bin(postenergies, bins)
    gsfcw, gsfcv, gsfcmodI = bin(gsfc_energy, bins)
    gsfcw, gsfcv, gsfcmeasI = bin(gsfc_energy, bins)
    
    prebintransmits = zeros(length(bins) - 1)
    postbintransmits = zeros(length(bins) - 1)
    gsfcmodbintransmits = zeros(length(bins) - 1)
    gsfcmeasbintransmits = zeros(length(bins) - 1)
    for i in eachindex(preI)
        rolling = 0
        for j in 1:length(preI[i])
            rolling += pretransmits[preI[i][j]]
        end
        prebintransmits[i] = rolling / length(preI[i])
    end
    for i in eachindex(postI)
        rolling = 0
        for j in 1:length(postI[i])
            rolling += posttransmits[postI[i][j]]
        end
        postbintransmits[i] = rolling / length(postI[i])
    end
    for i in eachindex(gsfcmodI)
        rolling = 0
        for j in 1:length(gsfcmodI[i])
            rolling += gsfc_mod[gsfcmodI[i][j]]
        end
        gsfcmodbintransmits[i] = rolling / length(gsfcmodI[i])
    end
    for i in eachindex(gsfcmeasI)
        rolling = 0
        for j in 1:length(gsfcmeasI[i])
            rolling += gsfc_mod[gsfcmeasI[i][j]]
        end
        gsfcmeasbintransmits[i] = rolling / length(gsfcmeasI[i])
    end
    
    csvdata = "binstart(eV),binstop(eV),precount,postcount,pretransmission,posttransmission"
    for k in 1:length(bins) - 1
        csvdata *= "\n" * 
            string(bins[k]) * "," * 
            string(bins[k + 1]) * "," *
            string(prew[k]) * "," *
            string(postw[k]) * "," *
            string(prebintransmits[k]) * "," *
            string(postbintransmits[k])
    end
    
    ax_each = CairoMakie.Axis(
        fig[1,1], 
        xlabel=L"\text{Energy [keV]}", 
        ylabel=L"\text{Transmission}", 
        title=L"\text{Raw transmission}",
        # xtickformat = values -> [L"%$(value)" for value in values],
        xticklabelsize = 12.f0,
        yticklabelsize = 12.f0
    )
    ax_compare = CairoMakie.Axis(
        fig[2,1], 
        xlabel=L"\text{Energy [keV]}", 
        ylabel=L"\frac{T_1}{T_2} - 1", 
        title=L"\text{Transmission comparison, } n=%$n_str_latex",
        xticklabelsize = 12.f0,
        yticklabelsize = 12.f0
    )
    
    stairs!(ax_each, bins[1:end-1] ./ 1e3, gsfcmeasbintransmits, linewidth=1, color=:purple, label=L"\text{GSFC measurement}")
    stairs!(ax_each, gsfc_energy ./ 1e3, gsfc_mod, linewidth=1, color=:green, label=L"\text{GSFC model}")
    stairs!(ax_each, bins[1:end-1] ./ 1e3, prebintransmits, linewidth=1, color=:blue, label=L"\text{model: parallel rays}")
    stairs!(ax_each, bins[1:end-1] ./ 1e3, postbintransmits, linewidth=1, color=:red, label=L"\text{model: post-optic rays}")
    leg_each = CairoMakie.Legend(fig[1,2], ax_each, frame_visible=false, labelsize=10.f0)
    
    stairs!(ax_compare, bins[1:end-1] ./ 1e3, prebintransmits ./ gsfcmodbintransmits .- 1, linewidth=1, label=L"\text{Thanasi parallel rays vs. GSFC model}")
    stairs!(ax_compare, bins[1:end-1] ./ 1e3, postbintransmits ./ prebintransmits .- 1, linewidth=1, label=L"\text{Thanasi parallel rays vs. post-optic rays}")
    stairs!(ax_compare, bins[1:end-1] ./ 1e3, postbintransmits ./ gsfcmodbintransmits .- 1, linewidth=1, label=L"\text{Thanasi post-optic rays vs. GSFC model}")
    leg_compare = CairoMakie.Legend(fig[2,2], ax_compare, frame_visible=false, labelsize=10.f0)
    
    
    
    ylims!(ax_compare, [-1, 1])
    xlims!(ax_compare, [1, 29] .* 1e0)
    xlims!(ax_each, [1, 29] .* 1e0)
    ylims!(ax_each, [-0.05, 1.1])
    
    fname_index = ""
    if !isnothing(index)
        fname_index = "_r"*string(index)
    end
    
    csvfname = joinpath(path, "transmit_data_n" * n_str * fname_index * ".csv")
    write(csvfname, csvdata)
    save(joinpath(path, "diff_transmit_n" * n_str * fname_index * ".pdf"), fig)
    
    return csvfname
end

function save_products(path::String, source_data::String; do_gsfc_meas=false)
    data = CSV.File(source_data, header=true)
    binstop = []
    binstart = []
    precount = []
    postcount = []
    pretransmits = []
    posttransmits = []
    for row in data
        push!(binstop, row[1])
        push!(binstart, row[2])
        push!(precount, row[3])
        push!(postcount, row[4])
        push!(pretransmits, row[5])
        push!(posttransmits, row[6])
    end
    
    CairoMakie.activate!()
    fig = CairoMakie.Figure(size=(600, 600), backgroundcolor=:transparent)
    
    gsfc_data = get_empirical_attenuation_data("src/data/20240607_foxsi4_transmission.csv")
    gsfc_energy = gsfc_data["energies"] .* 1e3
    gsfc_mod = gsfc_data["modeled_attenuations"]
    gsfc_meas = gsfc_data["measured_attenuations"]
    
    bins = sort!(unique!([binstart; binstop]))
    bins = convert(Vector{Float64}, bins)
    
    gsfcmodbintransmits = zeros(length(bins) - 1)
    gsfcmeasbintransmits = zeros(length(bins) - 1)
    gsfcw, gsfcv, gsfcmodI = bin(gsfc_energy, bins)
    gsfcw, gsfcv, gsfcmeasI = bin(gsfc_energy, bins)
    for i in eachindex(gsfcmodI)
        rolling = 0
        for j in 1:length(gsfcmodI[i])
            rolling += gsfc_mod[gsfcmodI[i][j]]
        end
        gsfcmodbintransmits[i] = rolling / length(gsfcmodI[i])
    end
    for i in eachindex(gsfcmeasI)
        rolling = 0
        for j in 1:length(gsfcmeasI[i])
            rolling += gsfc_meas[gsfcmeasI[i][j]]
        end
        gsfcmeasbintransmits[i] = rolling / length(gsfcmeasI[i])
    end
    
    # count photons for reference
    n_str_in= split(source_data, '.')[1]
    n_str_reg = r"(?<=_n)(.*)"
    n_str = match(n_str_reg, n_str_in).match
    mantissa, exponent = split(n_str, 'e')
    n_str_latex = L"%$mantissa \times 10^{%$exponent}"
    
    prew = precount
    postw = postcount
    prebintransmits = pretransmits
    postbintransmits = posttransmits
    
    ax_each = CairoMakie.Axis(
        fig[1,1], 
        xlabel=L"\text{Energy [keV]}", 
        ylabel=L"\text{Transmission}", 
        title=L"\text{Raw transmission}",
        xticklabelsize = 12.f0,
        yticklabelsize = 12.f0
    )
    ax_compare = CairoMakie.Axis(
        fig[2,1], 
        xlabel=L"\text{Energy [keV]}", 
        ylabel=L"\frac{T_1}{T_2} - 1", 
        title=L"\text{Transmission comparison, } n=%$n_str_latex", 
        yminorticks=IntervalsBetween(5), 
        yminorticksvisible=true,
        xticklabelsize = 12.f0,
        yticklabelsize = 12.f0
    )
    
    # stairs!(ax_each, bins[1:end-1], gsfcmeasbintransmits, linewidth=1, color=:purple, label=L"\text{GSFC measurement}")
    stairs!(ax_each, gsfc_energy ./ 1e3, gsfc_mod, linewidth=1, color=:green, label=L"\text{GSFC model}")
    stairs!(ax_each, bins[1:end-1] ./ 1e3, prebintransmits, linewidth=1, color=:blue, label=L"\text{model: parallel rays}")
    stairs!(ax_each, bins[1:end-1] ./ 1e3, postbintransmits, linewidth=1, color=:red, label=L"\text{model: post-optic rays}")
    
    
    stairs!(ax_compare, bins[1:end-1] ./ 1e3, prebintransmits ./ gsfcmodbintransmits .- 1, linewidth=1, label=L"\text{Thanasi parallel rays vs. GSFC model}")
    stairs!(ax_compare, bins[1:end-1] ./ 1e3, postbintransmits ./ prebintransmits .- 1, linewidth=1, label=L"\text{Thanasi parallel rays vs. post-optic rays}")
    stairs!(ax_compare, bins[1:end-1] ./ 1e3, postbintransmits ./ gsfcmodbintransmits .- 1, linewidth=1, label=L"\text{Thanasi post-optic rays vs. GSFC model}")
    # stairs!(ax_compare, bins[1:end-1], gsfcmeasbintransmits ./ gsfcmodbintransmits .- 1, linewidth=1, label=L"\text{GSFC measurement vs. GSFC model}")
    
    if do_gsfc_meas
        bg = parse_mca("src/data/20250515_163020_A1026d2_30kV_132uA_unfilt2.mca")
        fg = parse_mca("src/data/20250515_160110_A1026d2_30kV_132uA_center.mca")
        
        bg_en = bg["fit_const"] .+ bg["fit_slope"].*Vector(1:length(bg["adc"]))        
        fg_en = fg["fit_const"] .+ fg["fit_slope"].*Vector(1:length(fg["adc"]))
        bg_en *= 1e3
        fg_en *= 1e3
        println("range: ", min(bg_en...), ", ", max(bg_en...))
        
        # then bin bg_en, fg_en, then do binwise transmission.
        gsfcbgw, gsfcbgv, gsfcmeasbgI = bin(bg_en, bins)
        gsfcfgw, gsfcfgv, gsfcmeasfgI = bin(fg_en, bins)
        
        gsfcmeasbinbgtransmits = zeros(length(bins) - 1)
        gsfcmeasbinfgtransmits = zeros(length(bins) - 1)
        for i in eachindex(gsfcmeasbgI)
            rolling = 0
            for j in 1:length(gsfcmeasbgI[i])
                rolling += bg["adc"][gsfcmeasbgI[i][j]]
            end
            gsfcmeasbinbgtransmits[i] = rolling / length(gsfcmeasbgI[i])
        end
        for i in eachindex(gsfcmeasfgI)
            rolling = 0
            for j in 1:length(gsfcmeasfgI[i])
                rolling += fg["adc"][gsfcmeasfgI[i][j]]
            end
            gsfcmeasbinfgtransmits[i] = rolling / length(gsfcmeasfgI[i])
        end
        
        gsfclongmeastransmits = (gsfcmeasbinfgtransmits/fg["livetime"]) ./ (gsfcmeasbinbgtransmits/bg["livetime"])
        
        stairs!(ax_each, bins[1:end-1] ./ 1e3, gsfclongmeastransmits, linewidth=1, label=L"\text{GSFC long measurement}", color=:purple)
        
        stairs!(ax_compare, bins[1:end-1] ./ 1e3, gsfclongmeastransmits ./ gsfcmodbintransmits .- 1, linewidth=1, label=L"\text{GSFC long measurement vs. GSFC model}")
        
        # stairs!(ax_compare, bins[1:end-1], postbintransmits ./ gsfclongmeastransmits .- 1, linewidth=1, label=L"\text{Thanasi post-optic rays vs. GSFC long measurement}")
    
        # gsfcmeastransmits2 = (fg_en/fg["livetime"])./(bg_en/bg["livetime"])
    end
    leg_each = CairoMakie.Legend(fig[1,2], ax_each, frame_visible=false, labelsize=10.f0)
    leg_compare = CairoMakie.Legend(fig[2,2], ax_compare, frame_visible=false, labelsize=10.f0)
    
    ylims!(ax_compare, [-0.5, 0.5])
    xlims!(ax_compare, [1, 29] .* 1e0)
    xlims!(ax_each, [1, 29] .* 1e0)
    ylims!(ax_each, [-0.05, 1.1])
    
    save(joinpath(path, "diff_transmit_n" * n_str * ".pdf"), fig)
end

#######################################################################
### ----------------- for many, self-generated photons ------------ ###
#######################################################################

function main_nphoton()
    
    mon = make_monolithic()
    att = make_attenuator()
    println("thickness: ", att.topplane.a'*(att.topplane.c - att.bottomplane.c))
    n = Int(1e5)

    # if we want to tip all photons off-axis:
    tang = 10*pi/180
    offax_transform = [1 0 0; 0 cos(tang) -sin(tang); 0 sin(tang) cos(tang)]
    
    angphotons, angenergies, angangles = setup_n_photons(att, "src/data/milo_input.csv", n; transform=offax_transform)
    prephotons, preenergies, preangles = setup_n_photons(att, nothing, n)
    println("Mean angle [deg]: ", sum(angangles)/length(angangles))
    
    GLMakie.activate!()
    f = GLMakie.Figure(size=(900, 600))
    layout3d = GLMakie.GridLayout(f[1,1])
    layout2d = GLMakie.GridLayout(f[1,2])
    scene = GLMakie.LScene(
        layout3d[1,1],
        show_axis=false
    )
    
    plot_attenuator!(scene, att)
    XRayQuadrics.plot!(scene, angphotons, length_scale=1200e-6)
    # XRayQuadrics.plot!(scene, badphotons, length_scale=1200e-6)
    # XRayQuadrics.plot!(scene, parphotons, length_scale=1200e-6)
    
    display(f)
    
    preenergies = [p.E for p in prephotons]
    angenergies = [p.E for p in angphotons]
    
    # ax_energyangle = Axis(layout2d[1,1], ylabel="Angle [deg]", xlabel="Energy [eV]", title="Incident photon distribution")
    # ax_transmit = Axis(layout2d[2,1], ylabel="Transmission", xlabel="Energy [eV]", title="Pixelated attenuator transmission")
    # ax_compare = Axis(layout2d[3,1], ylabel=L"\frac{T(\theta) - T(0)}{T(0)}", xlabel="Energy [eV]", title="Comparative transmission")
    # ax_count = Axis(layout2d[4,1], ylabel="Total photons\nper bin", xlabel="Energy [eV]", title="Bin size")
    # ax_bin = Axis(layout2d[3,1], ylabel="Transmission", xlabel="Energy [eV]", title="Pixelated attenuator transmission")
    
    # scatter!(ax_energyangle, angenergies, angangles, markersize=1, color=:black)
    # scatter!(ax_energyangle, preenergies, preangles, markersize=1, color=:green)
    
    # println("Tracing post-optic photons through monolithic attenuator...")
    # basetransmit = batch_photons_through_attenuator(angphotons, mon)
    println("Tracing parallel photons through pixelated attenuator...")
    pretransmits = batch_photons_through_attenuator(prephotons, att)
    println("Tracing post-optic photons through pixelated attenuator...")
    angtransmits = batch_photons_through_attenuator(angphotons, att)
    
    # scatter!(ax_transmit, preenergies, pretransmits, markersize=1, color=:green, label="pre-optic\ntransmission")
    # scatter!(ax_transmit, angenergies, angtransmits, markersize=1, color=:black, label="post-optic\ntransmission")
    
    savedata = Dict(
        "bins" => Vector(1:0.25:29).*1e3,
        "pretransmits" => pretransmits,
        "posttransmits" => angtransmits,
        "prephotons" => prephotons,
        "postphotons" => angphotons
    )
    
    save_products("/home/thanasi/projects/ray/results/2025/nov12/f32", savedata, do_gsfc_meas=true)
end


function batch_random_photons(savepath::String, nphoton::Int, nrep::Int)
    # mon = make_monolithic()
    # att = make_monolithic()
    # att = make_monolithic_al()
    att = make_attenuator()
    # att = make_attenuator_with_hole()
    println("thickness: ", att.topplane.a'*(att.topplane.c - att.bottomplane.c))
    n = nphoton
    
    # if we want to tip all photons off-axis:
    tang = 0*pi/180
    offax_transform = [1 0 0; 0 cos(tang) -sin(tang); 0 sin(tang) cos(tang)]
    
    csvnames = String[]
    angphotons = Vector{Particle}(undef, n)
    prephotons = Vector{Particle}(undef, n)
    angenergies = zeros(n)
    preenergies = zeros(n)
    angangles = zeros(n)
    preangles = zeros(n)
    angtransmits = zeros(n)
    pretransmits = zeros(n)
    savedata = Dict(
        "bins" => Vector(1:0.25:29).*1e3,
        "pretransmits" => pretransmits,
        "posttransmits" => angtransmits,
        "prephotons" => prephotons,
        "postphotons" => angphotons
    )
    for k in 1:nrep
        setup_n_photons!(angphotons, angenergies, angangles, att, "src/data/milo_input.csv", n; transform=offax_transform, hires=true)
        setup_n_photons!(prephotons, preenergies, preangles, att, nothing, n)
        
        # preenergies = [p.E for p in prephotons]
        # angenergies = [p.E for p in angphotons]
        
        println("Tracing parallel photons through pixelated attenuator...")
        pretransmits = batch_photons_through_attenuator(prephotons, att)
        println("Tracing post-optic photons through pixelated attenuator...")
        angtransmits = batch_photons_through_attenuator(angphotons, att)
        
        println("Saving...")
        savedata["pretransmits"] = pretransmits
        savedata["posttransmits"] = angtransmits
        savedata["prephotons"] = prephotons
        savedata["postphotons"] = angphotons
        
        fname = save_products(savepath, savedata; do_gsfc_meas=true, index=k)
        push!(csvnames, fname)
    end
    # println(csvnames)
    mergefile = merge_csv(csvnames)
    
    save_products(savepath, mergefile; do_gsfc_meas=true)
    
    GC.gc()
end

function get_previous_photon_data(file::String)
    data = CSV.file(file, header=true)
    binstop = []
    binstart = []
    precount = []
    postcount = []
    pretransmission = []
    posttransmission = []
    for row in data
        push!(binstop, row[1])
        push!(binstart, row[2])
        push!(precount, row[3])
        push!(postcount, row[4])
        push!(pretransmission, row[5])
        push!(posttransmission, row[6])
    end
    return Dict(
        "binstop"=>binstop,
        "binstart"=>binstart,
        "precount"=>precount,
        "postcount"=>postcount,
        "pretransmits"=>pretransmission,
        "posttransmits"=>posttransmission
    )
end

function just_plot()
   save_products("/home/thanasi/projects/ray/results/2025/nov14/pix/lbnl/","/home/thanasi/projects/ray/results/2025/nov14/pix/lbnl/transmit_data_merged_n4e+07.csv"; do_gsfc_meas=true)
    
end

#batch_random_photons("/home/thanasi/projects/ray/results/2025/nov14/pix/nist/nagoya/trial", Int(1e4), 3)

