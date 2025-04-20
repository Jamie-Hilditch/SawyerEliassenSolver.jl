const TEST_FORCING_GRID_NX::Int = 4
const TEST_FORCING_GRID_NZ::Int = 4

function create_test_forcing_domain(FT)
    grid = Grid(FT, TEST_FORCING_GRID_NX, TEST_FORCING_GRID_NZ, 1, (zero(FT), π))
    return Domain(grid)
end

# forcings functions
pointwise_physical_forcing_with_parameters(x, z, t, p) = p.amp * sin(z)
pointwise_physical_forcing_without_parameters(x, z, t) = sin(z)
function global_physical_forcing_with_parameters(F, t, p)
    x, z = gridpoints(F)
    @inbounds @. F = p.amp * sin(z)
    return nothing
end
function global_physical_forcing_without_parameters(F, t)
    x, z = gridpoints(F)
    @inbounds @. F = sin(z)
    return nothing
end
function pointwise_spectral_forcing_with_parameters(kx, kz, t, p)
    return TEST_FORCING_GRID_NX *
           TEST_FORCING_GRID_NZ *
           p.amp *
           ifelse(kz == one(kz), one(kz), zero(kz))
end
function pointwise_spectral_forcing_without_parameters(kx, kz, t)
    return TEST_FORCING_GRID_NX *
           TEST_FORCING_GRID_NZ *
           ifelse(kz == one(kz), one(kz), zero(kz))
end
function global_spectral_forcing_with_parameters(F, t, p)
    F .= zero(eltype(F))
    @inbounds F[1, 1] = p.amp * prod(size(domain))
    return nothing
end
function global_spectral_forcing_without_parameters(F, t)
    F .= zero(eltype(F))
    @inbounds F[1, 1] = one(eltype(F)) * prod(size(domain))
    return nothing
end

# expected outputs
function physical_forcing_expected_output(domain::Domain{T}) where {T}
    x, z = gridpoints(domain)
    out = XZVariable(domain)
    @inbounds @. out = sin(z)
    return out
end

function spectral_forcing_expected_output(domain::Domain{T}) where {T}
    out = FSVariable(domain)
    @inbounds out[1, 1] = one(T) * prod(size(domain))
    return out
end

# test functions
function test_forcing_constructors(FT)
    domain = create_test_forcing_domain(FT)
    @test PointwisePhysicalForcing(
        domain, pointwise_physical_forcing_with_parameters, (; amp=one(FT))
    ) isa PointwisePhysicalForcing{FT}
    @test PointwisePhysicalForcing(
        domain, pointwise_physical_forcing_without_parameters
    ) isa PointwisePhysicalForcing{FT,Nothing}
    @test GlobalPhysicalForcing(
        domain, global_physical_forcing_with_parameters, (; amp=one(FT))
    ) isa GlobalPhysicalForcing{FT}
    @test GlobalPhysicalForcing(domain, global_physical_forcing_without_parameters) isa
        GlobalPhysicalForcing{FT,Nothing}
    @test PointwiseSpectralForcing(
        domain, pointwise_spectral_forcing_with_parameters, (; amp=one(FT))
    ) isa PointwiseSpectralForcing{FT}
    @test PointwiseSpectralForcing(
        domain, pointwise_spectral_forcing_without_parameters
    ) isa PointwiseSpectralForcing{FT,Nothing}
    @test GlobalSpectralForcing(
        domain, global_spectral_forcing_with_parameters, (; amp=one(FT))
    ) isa GlobalSpectralForcing{FT}
    @test GlobalSpectralForcing(domain, global_spectral_forcing_without_parameters) isa
        GlobalSpectralForcing{FT,Nothing}
end

function test_evaluate_physical_forcing_pointwise_physical_forcing(
    domain::Domain{T}
) where {T}
    out = XZVariable(domain)
    F = PointwisePhysicalForcing(
        domain, pointwise_physical_forcing_with_parameters, (; amp=one(T))
    )
    @test evaluate_physical_forcing!(F, out, one(T)) === nothing
    @test out ≈ physical_forcing_expected_output(domain)
end

function test_evaluate_physical_forcing_global_physical_forcing(domain::Domain{T}) where {T}
    out = XZVariable(domain)
    F = GlobalPhysicalForcing(
        domain, global_physical_forcing_with_parameters, (; amp=one(T))
    )
    @test evaluate_physical_forcing!(F, out, one(T)) === nothing
    @test out ≈ physical_forcing_expected_output(domain)
end

function test_evaluate_physical_forcing_pointwise_spectral_forcing(
    domain::Domain{T}
) where {T}
    out = XZVariable(domain)
    F = PointwiseSpectralForcing(
        domain, pointwise_spectral_forcing_with_parameters, (; amp=one(T))
    )
    @test_throws MethodError evaluate_physical_forcing!(F, out, one(T))
end

function test_evaluate_physical_forcing_global_spectral_forcing(domain::Domain{T}) where {T}
    out = XZVariable(domain)
    F = GlobalSpectralForcing(
        domain, global_spectral_forcing_with_parameters, (; amp=one(T))
    )
    @test_throws MethodError evaluate_physical_forcing!(F, out, one(T))
end

function test_evaluate_physical_forcing(FT)
    domain = create_test_forcing_domain(FT)
    test_evaluate_physical_forcing_pointwise_physical_forcing(domain)
    test_evaluate_physical_forcing_global_physical_forcing(domain)
    test_evaluate_physical_forcing_pointwise_spectral_forcing(domain)
    return test_evaluate_physical_forcing_global_spectral_forcing(domain)
end

function test_evaluate_ζ_forcing_pointwise_physical_forcing(domain::Domain{T}) where {T}
    out = FSVariable(domain)
    F = PointwisePhysicalForcing(
        domain, pointwise_physical_forcing_with_parameters, (; amp=one(T))
    )
    @test evaluate_ζ_forcing!(F, out, one(T), XSVariable(domain), XZVariable(domain)) ===
        nothing
    @test out ≈ spectral_forcing_expected_output(domain)
end

function test_evaluate_ζ_forcing(FT)
    domain = create_test_forcing_domain(FT)
    return test_evaluate_ζ_forcing_pointwise_physical_forcing(domain)
end

@testset "Forcing constructors" for FT in FLOAT_TYPES
    @info "\tTesting forcing constructors for $FT"
    test_forcing_constructors(FT)
end

@testset "Evaluate physical forcing" for FT in FLOAT_TYPES
    @info "\tTesting evaluate physical forcing for $FT"
    test_evaluate_physical_forcing(FT)
end

@testset "Evaluate ζ forcing" for FT in FLOAT_TYPES
    @info "\tTesting evaluate ζ forcing for $FT"
    test_evaluate_ζ_forcing(FT)
end
