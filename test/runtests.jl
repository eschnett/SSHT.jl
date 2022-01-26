using Random
using Test
using ssht

chop(x) = abs2(x) < eps(x) ? zero(x) : x
chop(x::Complex) = Complex(chop(real(x)), chop(imag(x)))

function integrate(f::AbstractArray{T,2}, g::AbstractArray{U,2}, L::Integer) where {T,U}
    nphi = ssht.sampling_dh_nphi(L)
    ntheta = ssht.sampling_dh_ntheta(L)
    @assert size(f) == size(g) == (nphi, ntheta)

    s = zero(T) * zero(U) * zero(Float64)
    for p in 1:nphi, t in 1:ntheta
        phi = ssht.sampling_dh_p2phi(p, L)
        theta = ssht.sampling_dh_t2theta(t, L)
        dtheta = ssht.sampling_weight_dh(theta, L)
        # dtheta = sin(theta) * π / ntheta
        dphi = 2π / nphi
        s += conj(f[p, t]) * g[p, t] * dtheta * dphi
    end

    return s
end

################################################################################

@testset "Grid point counts" begin
    for L in 1:10
        n = ssht.sampling_dh_n(L)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)
        @test nphi == 2 * L - 1
        @test ntheta == 2 * L
        @test n == nphi * ntheta
    end
end

@testset "Mode layout and counts" begin
    for L in 1:10
        mode_seen = falses(L^2)
        for l in 0:(L - 1), m in (-l):l
            ind = ssht.sampling_elm2ind(l, m)
            @test 1 ≤ ind ≤ L^2
            @test ssht.sampling_ind2elm(ind) == (l, m)
            @test !mode_seen[ind]
            mode_seen[ind] = true
        end
        @test all(mode_seen)
    end
end

@testset "Coordinates" begin
    for L in 1:100
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
modes = [(name="(0,0)", fun=(p, t) -> sqrt(1 / 4π), modes=[1, 0, 0, 0]),
         (name="(0,1)", fun=(theta, phi) -> -sqrt(3 / 2π) * sin(theta) * cos(phi), modes=[0, -1, 0, 1]),
         (name="(0,-1)", fun=(theta, phi) -> sqrt(3 / 2π) * sin(theta) * sin(phi), modes=[0, im, 0, im]),
         (name="(1,0)", fun=(theta, phi) -> sqrt(3 / 4π) * cos(theta), modes=[0, 0, 1, 0])]
@testset "Simple real transforms: $(mode.name)" for mode in modes
    verbosity = 0
    for iter in 1:100
        L = rand(2:20)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        f = Array{Float64}(undef, nphi, ntheta)
        for p in 1:nphi, t in 1:ntheta
            phi = ssht.sampling_dh_p2phi(p, L)
            theta = ssht.sampling_dh_t2theta(t, L)
            f[p, t] = mode.fun(theta, phi)
        end

        flm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov_real!(flm, f, L, verbosity)
        @assert L ≥ 2
        @test flm[1:4] ≈ mode.modes
        if L > 2
            @test all(isapprox(0; atol=100eps()), flm[5:end])
        end

        f′ = similar(f)
        ssht.core_dh_inverse_sov_real!(f′, flm, L, verbosity)
        @test isapprox(f′, f; atol=1000eps())
    end
end

Random.seed!(100)
modes = [(name="(0,0,0)", spin=0, fun=(p, t) -> sqrt(1 / 4π), modes=[1, 0, 0, 0]),
         (name="(0,0,1)", spin=0, fun=(theta, phi) -> -sqrt(3 / 8π) * sin(theta) * cis(phi), modes=[0, 0, 0, 1]),
         (name="(0,0,-1)", spin=0, fun=(theta, phi) -> sqrt(3 / 8π) * sin(theta) * cis(-phi), modes=[0, 1, 0, 0]),
         (name="(0,1,0)", spin=0, fun=(theta, phi) -> sqrt(3 / 4π) * cos(theta), modes=[0, 0, 1, 0]),
         (name="(1,0,1)", spin=1, fun=(theta, phi) -> -sqrt(3 / 16π) * (1 - cos(theta)) * cis(phi), modes=[0, 0, 0, 1]),
         (name="(1,0,-1)", spin=1, fun=(theta, phi) -> -sqrt(3 / 16π) * (1 + cos(theta)) * cis(-phi), modes=[0, 1, 0, 0]),
         (name="(1,1,0)", spin=1, fun=(theta, phi) -> sqrt(3 / 8π) * sin(theta), modes=[0, 0, 1, 0])]
@testset "Simple complex transforms: $(mode.name) spin=$(mode.spin)" for mode in modes
    verbosity = 0
    for iter in 1:100
        L = rand(2:20)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        f = Array{Complex{Float64}}(undef, nphi, ntheta)
        for p in 1:nphi, t in 1:ntheta
            phi = ssht.sampling_dh_p2phi(p, L)
            theta = ssht.sampling_dh_t2theta(t, L)
            f[p, t] = mode.fun(theta, phi)
        end

        flm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov!(flm, f, L, mode.spin, verbosity)
        @assert L ≥ 2
        @test flm[1:4] ≈ mode.modes
        if L > 2
            @test all(isapprox(0; atol=100eps()), flm[5:end])
        end

        f′ = similar(f)
        ssht.core_dh_inverse_sov!(f′, flm, L, mode.spin, verbosity)
        @test isapprox(f′, f; atol=1000eps())
    end
end

Random.seed!(100)
@testset "Linearity of real transforms" begin
    verbosity = 0
    for iter in 1:100
        L = rand(1:100)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        f = randn(Float64, nphi, ntheta)
        g = randn(Float64, nphi, ntheta)
        α = randn(Float64)
        h = f + α * g

        flm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov_real!(flm, f, L, verbosity)
        glm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov_real!(glm, g, L, verbosity)
        hlm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov_real!(hlm, h, L, verbosity)

        @test flm + α * glm ≈ hlm
    end
end

Random.seed!(100)
@testset "Linearity of complex transforms" begin
    verbosity = 0
    for iter in 1:100
        L = rand(1:100)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        spin = rand(-4:4)

        f = randn(Complex{Float64}, nphi, ntheta)
        g = randn(Complex{Float64}, nphi, ntheta)
        α = randn(Complex{Float64})
        h = f + α * g

        flm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov!(flm, f, L, spin, verbosity)
        glm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov!(glm, g, L, spin, verbosity)
        hlm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov!(hlm, h, L, spin, verbosity)

        @test flm + α * glm ≈ hlm

        flm = randn(Complex{Float64}, L^2)
        glm = randn(Complex{Float64}, L^2)
        hlm = flm + α * glm

        ssht.core_dh_inverse_sov!(f, flm, L, spin, verbosity)
        ssht.core_dh_inverse_sov!(g, glm, L, spin, verbosity)
        ssht.core_dh_inverse_sov!(h, hlm, L, spin, verbosity)

        @test f + α * g ≈ h
    end
end

Random.seed!(100)
@testset "Orthonormality of complex transforms" begin
    verbosity = 0
    for iter in 1:100
        L = rand(1:100)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        spin = 0rand(-4:4)

        flm = zeros(Complex{Float64}, L^2)
        glm = zeros(Complex{Float64}, L^2)

        lf = rand(0:(L - 1))
        mf = rand((-lf):lf)
        lg = rand(0:(L - 1))
        mg = rand((-lg):lg)

        flm[ssht.sampling_elm2ind(lf, mf)] = 1
        glm[ssht.sampling_elm2ind(lg, mg)] = 1

        f = Array{Complex{Float64}}(undef, nphi, ntheta)
        ssht.core_dh_inverse_sov!(f, flm, L, spin, verbosity)
        g = Array{Complex{Float64}}(undef, nphi, ntheta)
        ssht.core_dh_inverse_sov!(g, glm, L, spin, verbosity)

        @test isapprox(integrate(f, f, L), 1; atol=1 / L^2)
        @test isapprox(integrate(f, g, L), (lf == lg) * (mf == mg); atol=1 / L^2)

        h = conj(f) .* f
        hlm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov!(hlm, h, L, spin, verbosity)
        @test isapprox(hlm[ssht.sampling_elm2ind(0, 0)], sqrt(1 / 4π); atol=sqrt(eps()))

        h = conj(f) .* g
        ssht.core_dh_forward_sov!(hlm, h, L, spin, verbosity)
        @test isapprox(hlm[ssht.sampling_elm2ind(0, 0)], (lf == lg) * (mf == mg) * sqrt(1 / 4π); atol=sqrt(eps()))
    end
end

@testset "derivatives (eth, eth-bar)" begin end
