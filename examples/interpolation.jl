using XRayQuadrics
using GLMakie
using LinearAlgebra

function interpolation()
    data = get_mass_attenuation_data("src/data/LBNL_attenlength_Si.csv")
    rawenergies = data["energies"]
    rawattenuations = data["attenuations"]
    
    interpenergies = Vector(1:0.01:29)*1e3
    interpattenuation = zeros(length(interpenergies))
    for (k,energy) in enumerate(interpenergies)
        interpattenuation[k] = interpolateattenuation(energy, data)
    end
    
    GLMakie.activate!()
    fig = GLMakie.Figure(size=(500, 800))
    ax = GLMakie.Axis(fig[1,1])
    scatter!(ax, rawenergies, rawattenuations, marker='+', markersize=10, color=:black)
    lines!(ax, interpenergies, interpattenuation, linewidth=1, color=:blue)
    
    display(fig)
end