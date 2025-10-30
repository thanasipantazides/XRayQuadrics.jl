using XRayQuadrics
using YAML

# optic table:
#   lookup shell
#       lookup station (in, out, inter)

function main()
    file = "src/data/prescription_sim.yml"
    optics_data = YAML.load_file(file)

    # shells = [8, 10]
    shells = 1:10
    y = [optics_data["axial"]["outer"]; optics_data["axial"]["inter"]; optics_data["axial"]["inner"]]
    for shell in shells
        this_lookup = optics_data["radial"][shell]
        x_p = [this_lookup["outer"]; this_lookup["inter"]]
        x_h = [this_lookup["inter"]; this_lookup["inner"]]
        resultp = fitp(x_p, y[1:2])
        resulth = fith(x_h, y[2:3])
        # println("p: ", resultp)
        # println("h: ", resulth)

        println("f_p: ", resultp.c + 1/4/resultp.a)
        println("f_h: ", resulth.c * (1 + 1/resultp.a))
    end
end