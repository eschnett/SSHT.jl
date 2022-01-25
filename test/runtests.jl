using Random
using Test
using ssht

@testset "Grid point counts" begin
    for L in 0:10
        n = ssht.sampling_dh_n(L)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)
        @test nphi == 2 * L - 1
        @test ntheta == 2 * L
        @test n == nphi * ntheta
    end
end

@testset "Coordinates" begin
    for L in 1:10
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)
        for p in 1:nphi
            phi = ssht.sampling_dh_p2phi(p, L)
            @test phi ≈ 2π * (p - 1) / nphi
        end
        for t in 1:ntheta
            theta = ssht.sampling_dh_t2theta(t, L)
            @test theta ≈ π * (t - 1 / 2) / ntheta
        end
    end
end

Random.seed!(100)
@testset "Simple real transforms: $mode" for mode in [:one, :x, :y, :z]
    verbosity = 0
    for iter in 1:10
        L = rand(2:20)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        f = Array{Float64}(undef, nphi, ntheta)
        for p in 1:nphi, t in 1:ntheta
            phi = ssht.sampling_dh_p2phi(p, L)
            theta = ssht.sampling_dh_t2theta(t, L)
            mode ≡ :one && (f[p, t] = 1)
            mode ≡ :x && (f[p, t] = sin(theta) * cos(phi))
            mode ≡ :y && (f[p, t] = sin(theta) * sin(phi))
            mode ≡ :z && (f[p, t] = cos(theta))
        end

        flm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov_real!(flm, f, L, verbosity)
        @assert L ≥ 2
        if mode ≡ :one
            @test flm[1:4] ≈ [sqrt(4π), 0, 0, 0]
        elseif mode ≡ :x
            @test flm[1:4] ≈ [0, sqrt(2π / 3), 0, -sqrt(2π / 3)]
        elseif mode ≡ :y
            @test flm[1:4] ≈ [0, sqrt(2π / 3) * im, 0, sqrt(2π / 3) * im]
        elseif mode ≡ :z
            @test flm[1:4] ≈ [0, 0, sqrt(4π / 3), 0]
        end
        if L > 2
            @test all(isapprox(0; atol=100eps()), flm[5:end])
        end

        f′ = similar(f)
        ssht.core_dh_inverse_sov_real!(f′, flm, L, verbosity)
        @test isapprox(f′, f; atol=1000eps())
    end
end

@testset "Simple complex transforms: $mode" for mode in [:one, :x, :y, :z]
    verbosity = 0
    spin = 0
    for iter in 1:10
        L = rand(2:20)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        f = Array{Complex{Float64}}(undef, nphi, ntheta)
        for p in 1:nphi, t in 1:ntheta
            phi = ssht.sampling_dh_p2phi(p, L)
            theta = ssht.sampling_dh_t2theta(t, L)
            mode ≡ :one && (f[p, t] = 1)
            mode ≡ :x && (f[p, t] = sin(theta) * cis(phi))
            mode ≡ :y && (f[p, t] = sin(theta) * cis(-phi))
            mode ≡ :z && (f[p, t] = cos(theta))
        end

        flm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov!(flm, f, L, spin, verbosity)
        @assert L ≥ 2
        if mode ≡ :one
            @test flm[1:4] ≈ [sqrt(4π), 0, 0, 0]
        elseif mode ≡ :x
            @test flm[1:4] ≈ [0, sqrt(2π / 3), 0, -sqrt(2π / 3)]
        elseif mode ≡ :y
            @test flm[1:4] ≈ [0, sqrt(2π / 3) * im, 0, sqrt(2π / 3) * im]
        elseif mode ≡ :z
            @test flm[1:4] ≈ [0, 0, sqrt(4π / 3), 0]
        end
        if L > 2
            @test all(isapprox(0; atol=100eps()), flm[5:end])
        end

        f′ = similar(f)
        ssht.core_dh_inverse_sov!(f′, flm, L, spin, verbosity)
        @test isapprox(f′, f; atol=1000eps())
    end
end

@testset "Linearity of transforms (real/complex)" begin end

@testset "nonzero spin" begin end

@testset "derivatives (eth, eth-bar)" begin end
