# SawyerEliassenSolver.jl

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Jamie-Hilditch.github.io/SawyerEliassenSolver.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Jamie-Hilditch.github.io/SawyerEliassenSolver.jl/dev/)
[![Build Status](https://github.com/Jamie-Hilditch/SawyerEliassenSolver/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Jamie-Hilditch/SawyerEliassenSolver/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)


A fast pseudospectral solver for the Sawyer-Eliassen equation. 

## About 

- Solves the forced Sawyer-Eliassen equation.
- Utilises pseudospectral discretisation and implicit time-stepping.
- Run relatively large simulations (e.g. 1024x1024 grid) in a matter of minutes on a laptop.

## Installation

SawyerEliassenSolver.jl can be installed using the Julia package manager.

```julia
julia> using Pkg

julia> Pkg.add("SawyerEliassenSolver")
```

> [!IMPORTANT]
> Julia 1.8 or later is required
>
> SawyerEliassenSolver.jl is tested on Julia 1.8--1.11.

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
