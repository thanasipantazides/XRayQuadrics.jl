# XRayQuadrics
Physical simulation of x-rays interacting with quadric surfaces.

Quadric surfaces look like this:
![Image of five quadric surfaces, cylinder, cone, paraboloid, hyperboloid, ellipsoid.](docs/src/assets/figures/shapes.png)

## Development
### Install Julia
Install the Julia language via `juliaup` (the installer), available [here](https://julialang.org/install/).

### Clone this repository
```bash
    git clone https://github.com/thanasipantazides/XRayQuadrics.jl.git
```

### Add required packages
Run
```bash
% julia --project=.
```
which launches the Julia REPL. From here, activate the project environment:
```bash
julia> ]
pkg> activate .
```
### Run unit tests
```bash
pkg> test
```

### Make a sample plot
```bash
... hit backspace to get out of the Pkg prompt...
julia> using Revise
julia> include("examples/plotshapes.jl")
julia> main()
```