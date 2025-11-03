# XRayQuadrics
Physical simulation of x-rays interacting with quadric surfaces.

Quadric surfaces look like this:
![Image of five quadric surfaces, cylinder, cone, paraboloid, hyperboloid, ellipsoid.](docs/src/assets/figures/shapes.png)

## Development
### Install Julia

### Clone this repository

### Add packages
```bash
% julia --project=.
julia> ]
pkg> activate .
```
### Run unit tests
```bash
pkg> test
```

### Make sample plot
```bash
... hit backspace to get out of Pkg ...
julia> using Revise
julia> include("examples/plotshapes.jl")
```