# Introduction

## Notation
`XRayQuadrics` is developed to analyze interactions between particles (such as photons) and physical surfaces, and to model the outcome of those interactions. In this documentation, we follow these conventions:
* Scalars are denoted with italic letters: ``a``, ``Z``, ``\gamma``.
* Vectors are denoted with bold lowercase: ``\boldsymbol{q}``, ``\boldsymbol{r}``.
* Matrices are denoted with bold uppercase: ``\boldsymbol{Q}``, ``\boldsymbol{A}``.
* The identity matrix is denoted: ``\boldsymbol{1}``.
* The zero matrix and zero vector are both denoted: ``\boldsymbol{0}``.

## Particles
For our purposes, a particle ``p`` is defined by an initial position ``\boldsymbol{r}_0``, a direction of travel (unit vector) ``\boldsymbol{v}``, and an energy ``E``:
```math
p \equiv \{\boldsymbol{r}_0, \boldsymbol{v}, E\}.
```
This set is stored as a `struct` in `Particle`. Then, any point along a particle's path is defined by some time-offset from the initial position, ``\Delta t`` (reflections will be dealt with later):
```math
\boldsymbol{r}(\Delta t) = \boldsymbol{r}_0 + \boldsymbold{v}\Delta t.
```
We will typically assume the particle is a photon.

## Surfaces
There are a few elementary shapes typically used in optical applications: lines, parabolas, hyperbolas, ellipses. These are all examples of *conic sections*, or slices of cones.
<!-- put figure -->
If any of these surfaces are rotated about their symmetry axis, they form their three-dimensional analogs: planes, paraboloids, hyperboloids, ellipsoids. All these shapes are classified as *quadric surfaces*. These elementary shapes are used in the construction of many practical optical systems. Conveniently, all quadric surfaces can be uniquely defined by a single ``4 \times 4`` matrix ``\boldsymbol{Q}``. The page [Quadric Surfaces](@ref) documents mathematical details of these quadric surfaces.