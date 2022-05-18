using XRayQuadrics
using GLMakie
using Random
using Distributions
import CSV,DataFrames

print(@__DIR__)

df = DataFrames.DataFrame(CSV.File(joinpath(@__DIR__,"../src/data/milo_input.csv")))

savdata = readsav(joinpath(@__DIR__,"../src/data/peaks_det101.sav"))

angles = df.angle*π/180     # in rad
energies = df.energy*1000   # in eV
N = length(angles)

α = 2*tan(mean(angles))
println(2*tan(mean(angles)))

sc = Scene()

fig,ax,scat = scatter(
    angles*180/π, 
    energies/1000, 
    markersize=0.5,
    axis=(
        xlabel="angle (deg)",
        ylabel="energy (keV)",
        title="Post-optic energy/angle distribution for flat input spectrum"
    )
)

Isort = sortperm(angles)
angles = angles[Isort]
energies = energies[Isort]

# fig,ax,line = lines(
#     angles*180/π,
#     energies/1000,
#     axis=(
#         xlabel="angle (deg)",
#         ylabel="energy (keV)",
#     )
# )

cam2d!(sc)

optax = [0;0;1]
focus = [0;0;0]

randorder = randperm(N)
angles = angles[randorder]
energies = energies[randorder]

# MAKE SURE angles IS RANDOMLY PERMUTED SO THIS IS NOT BIASING
negatives = [-1 for i in 1:floor(N/2)]
positives = [1 for i in ceil(N/2):N]
signs = cat(negatives, positives, dims=1)
angles = angles.*signs

v = zeros(3,N)
for i in 1:N
    v[:,i] = [sin(angles[i]); 0; cos(angles[i])]
end

# which positions along optic axis to sample rays
# zplanes = [0.1, 0.2, 0.3, 0.4]
exprange = -(1:1)
zplanes = 2.0.^exprange
zplanes = [0.1]

# store distribution of ray intersections (across X), for zplanes positions
spreads = zeros(length(zplanes), N)
for plane in 1:length(zplanes)
    spreads[plane, :] = v[1,:]*(zplanes[plane] - focus[3])
end

# scpfig,scpax,scp = scatter(ones(N), spreads[1,:])
hfig = Figure(resolution=(2000,1000))
hax = Axis(hfig[1,1], xlabel="position", ylabel="counts per bin")

# slider = Slider(hfig[2,1], range=1:1000, startvalue=1)

slidergrid = SliderGrid(
    hfig[2,1], (
        label="number of bins",   
        range=10.0.^(0:0.05:3) .+ 1, 
        startvalue=5,
        format=x -> string(round(x), " bins")
    )
)

bincount = lift(slidergrid.sliders[1].value) do x
    Int(round(x))
end

zv = 0.1
spread = v[1,:]*(zv - focus[3])

ticks = lift(bincount) do bincount
    Array(LinRange(min(spread...), max(spread...), bincount+1))
end

# centers = lift(ticks) do ticks
#     println(length(ticks))
#     ticks[1:end-1] .+ (ticks[2] - ticks[1])/2
# end

# weights = lift(ticks) do ticks
#     println(length(ticks))
#     bin(spread, ticks)[1]
# end

# weights = lift(bincount) do bincount
#     println(bincount)
#     bin(spread, bincount+1)[1]
# end

container = lift(ticks) do ticks
    [ticks[1:end-1] .+ (ticks[2] - ticks[1])/2, bin(spread, ticks)[1]]
end

# ticks[]
# centers[]
# weights[]


# @show weights
# @show ticks
# @show centers

xn = spread[spread .< 0]
xp = spread[spread .> 0]

σn = var(xn)
σp = var(xp)

μn = mean(xn)
μp = mean(xp)

histogram = hist!(
    spread,
    normalization=:pdf,
    bins=bincount
)

np = 1/σp/sqrt(2*π).*exp.(-((sort(xp) .- µp).^2.0)/2/σp^2)
nn = 1/σn/sqrt(2*π).*exp.(-((sort(xn) .- µn).^2.0)/2/σn^2)

findmax(np)
findmax(nn)

# PEAK-PEAK height: 0.0042

# lines!(sort(xp), np)

# bar = barplot!(
#     centers,
#     weights,
# )

# bar = barplot!(
#     x=container,
#     y=container,
# )

limits!(hax, min(spread...), max(spread...), 0, 1000)

# hfig,hax,hplt = hist(spreads[1,:],bins=100)
# save("rays_projected_along_optic_axis.png", hfig, pt_per_unit = 10)

# ISSUES:
#   Working with photons from ideal optics. If I throw them all at one focal point, 
#   as soon as I move away from that point there will be a hole in the distribution.
#   So the distribution is either degenerate (at the focus) or unrealistic (anywhere else).
#   If there is a model for photon exit position relative to the optic axis, I could
#   reconstruct a distribution.
hfig