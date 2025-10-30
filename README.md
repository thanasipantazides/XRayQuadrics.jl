# XRayQuadrics
Physical simulation of x-rays interacting with 

## Development
### Install Julia

### Clone this repository

### Add packages
```bash
% julia
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
julia> include("test/mplot.jl")
```
because I have written the plotting functionality at the file-level, not wrapped it in a module.