```@meta
CurrentModule = SawyerEliassenSolver
```

# SawyerEliassenSolver.jl

A fast pseudospectral solver for the Sawyer-Eliassen equation written in Julia. 

## About 

- Solves the [forced Sawyer-Eliassen equation](@ref "The Sawyer-Eliassen equation").
- Utilises [pseudospectral discretisation](@ref "Pseudospectral discretisation") and [implicit time-stepping](@ref "Diagonally Implicit Runge-Kutta Nyström").
- Run relatively large simulations (e.g. 1024x1024 grid) in a matter of minutes on a laptop.

## Installation

SawyerEliassenSolver.jl can be installed using the Julia package manager.

```julia
julia> using Pkg

julia> Pkg.add("SawyerEliassenSolver")
```

!!! compat "Required Julia version"
    Julia 1.8 or later is required

!!! info "Tested Julia versions"
    SawyerEliassenSolver.jl is tested on Julia 1.8--1.11.

## Examples

Check out the examples to see how it all works:
```@contents
Pages = [page.second for page in Main.example_pages]
Depth = 1
```



