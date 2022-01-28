using Random
using Test
using ssht

bitsign(b::Bool) = b ? -1 : 1
bitsign(i::Integer) = bitsign(isodd(i))

chop(x) = abs2(x) < 100eps(x) ? zero(x) : x
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
        dphi = 2π / nphi
        s += conj(f[p, t]) * g[p, t] * dtheta * dphi
    end

    return s
end

################################################################################

# Half angle formulae:
#     sin(θ/2)^2 = (1-cos(θ))/2
#     cos(θ/2)^2 = (1+cos(θ))/2

# These sYlm are taken from Wikipedia and black-holes.org:
# <https://en.wikipedia.org/wiki/Spherical_harmonics> and
# <https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics>
# <https://www.black-holes.org/SpinWeightedSphericalHarmonics.nb>
function sYlm(s::Integer, l::Integer, m::Integer, θ::Real, ϕ::Real)
    @assert abs(s) ≤ l
    @assert -l ≤ m ≤ l

    # Parity:
    s < 0 && return bitsign(s + m) * conj(sYlm(-s, l, -m, θ, ϕ))

    (s, l, m) == (0, 0, 0) && return sqrt(1 / 4π)
    (s, l, m) == (0, 1, -1) && return sqrt(3 / 8π) * sin(θ) * cis(-ϕ)
    (s, l, m) == (0, 1, 0) && return sqrt(3 / 4π) * cos(θ)
    (s, l, m) == (0, 1, +1) && return -sqrt(3 / 8π) * sin(θ) * cis(ϕ)
    (s, l, m) == (0, 2, -2) && return sqrt(15 / 2π) * cos(θ / 2)^2 * sin(θ / 2)^2 * cis(-2ϕ)
    (s, l, m) == (0, 2, -1) && return -sqrt(15 / 2π) * cos(θ / 2) * sin(θ / 2) * (-cos(θ / 2)^2 + sin(θ / 2)^2) * cis(-ϕ)
    (s, l, m) == (0, 2, 0) && return sqrt(5 / 4π) * (cos(θ / 2)^4 - 4 * cos(θ / 2)^2 * sin(θ / 2)^2 + sin(θ / 2)^4)
    (s, l, m) == (0, 2, +1) && return -sqrt(15 / 2π) * cos(θ / 2) * sin(θ / 2) * (cos(θ / 2)^2 - sin(θ / 2)^2) * cis(ϕ)
    (s, l, m) == (0, 2, +2) && return sqrt(15 / 2π) * cos(θ / 2)^2 * sin(θ / 2)^2 * cis(2ϕ)

    (s, l, m) == (+1, 1, -1) && return -sqrt(3 / 16π) * (1 + cos(θ)) * cis(-ϕ)
    (s, l, m) == (+1, 1, 0) && return sqrt(3 / 8π) * sin(θ)
    (s, l, m) == (+1, 1, +1) && return -sqrt(3 / 16π) * (1 - cos(θ)) * cis(ϕ)
    (s, l, m) == (+1, 2, -2) && return -sqrt(5 / π) * cos(θ / 2)^3 * sin(θ / 2) * cis(-2ϕ)
    (s, l, m) == (+1, 2, -1) && return -sqrt(5 / 4π) * cos(θ / 2)^2 * (cos(θ / 2)^2 - 3 * sin(θ / 2)^2) * cis(-ϕ)
    (s, l, m) == (+1, 2, 0) && return sqrt(15 / 2π) * cos(θ / 2) * sin(θ / 2) * (cos(θ / 2)^2 - sin(θ / 2)^2)
    (s, l, m) == (+1, 2, +1) && return sqrt(5 / 4π) * sin(θ / 2)^2 * (-3 * cos(θ / 2)^2 + sin(θ / 2)^2) * cis(ϕ)
    (s, l, m) == (+1, 2, +2) && return sqrt(5 / π) * cos(θ / 2) * sin(θ / 2)^3 * cis(2ϕ)

    (s, l, m) == (+2, 2, -2) && return sqrt(5 / 4π) * cos(θ / 2)^4 * cis(-2ϕ)
    (s, l, m) == (+2, 2, -1) && return -sqrt(5 / π) * cos(θ / 2)^3 * sin(θ / 2) * cis(-ϕ)
    (s, l, m) == (+2, 2, 0) && return sqrt(15 / 2π) * cos(θ / 2)^2 * sin(θ / 2)^2
    (s, l, m) == (+2, 2, +1) && return -sqrt(5 / π) * cos(θ / 2) * sin(θ / 2)^3 * cis(ϕ)
    (s, l, m) == (+2, 2, +2) && return sqrt(5 / 4π) * sin(θ / 2)^4 * cis(2ϕ)

    @assert false
end

# ð F = - (sin θ)^s (∂_θ + i / sin(θ) ∂_ϕ) (sin θ)^-s F
# (Calculated manually from sYlm above)
function ðsYlm(s::Integer, l::Integer, m::Integer, θ::Real, ϕ::Real)
    @assert abs(s) ≤ l
    @assert -l ≤ m ≤ l

    (s, l, m) == (0, 0, 0) && return 0
    (s, l, m) == (0, 1, -1) && return -sqrt(3 / 8π) * (1 + cos(θ)) * cis(-ϕ)
    (s, l, m) == (0, 1, 0) && return sqrt(3 / 4π) * sin(θ)
    (s, l, m) == (0, 1, +1) && return -sqrt(3 / 8π) * (1 - cos(θ)) * cis(ϕ)
    (s, l, m) == (0, 2, -2) && return -sqrt(15 / 8π) * (1 + cos(θ)) * sin(θ) * cis(-2ϕ)
    (s, l, m) == (0, 2, -1) && return -sqrt(15 / 2π) * cos(θ / 2)^2 * (-1 + 2 * cos(θ)) * cis(-ϕ)
    (s, l, m) == (0, 2, 0) && return sqrt(45 / 16π) * sin(2θ)
    (s, l, m) == (0, 2, +1) && return -sqrt(15 / 2π) * (1 + 2 * cos(θ)) * sin(θ / 2)^2 * cis(ϕ)
    (s, l, m) == (0, 2, +2) && return -sqrt(15 / 8π) * (-1 + cos(θ)) * sin(θ) * cis(2ϕ)

    (s, l, m) == (+1, 1, -1) && return 0
    (s, l, m) == (+1, 1, 0) && return 0
    (s, l, m) == (+1, 1, +1) && return 0
    (s, l, m) == (+1, 2, -2) && return sqrt(5 / π) * cos(θ / 2)^4 * cis(-2ϕ)
    (s, l, m) == (+1, 2, -1) && return -sqrt(5 / π) * cos(θ / 2)^2 * sin(θ) * cis(-ϕ)
    (s, l, m) == (+1, 2, 0) && return sqrt(15 / 8π) * sin(θ)^2
    (s, l, m) == (+1, 2, +1) && return -sqrt(5 / π) * sin(θ / 2)^2 * sin(θ) * cis(ϕ)
    (s, l, m) == (+1, 2, +2) && return sqrt(5 / π) * sin(θ / 2)^4 * cis(2ϕ)

    (s, l, m) == (+2, 2, -2) && return 0
    (s, l, m) == (+2, 2, -1) && return 0
    (s, l, m) == (+2, 2, 0) && return 0
    (s, l, m) == (+2, 2, +1) && return 0
    (s, l, m) == (+2, 2, +2) && return 0

    @assert false
end

# ð̄ F = - (sin θ)^-s (∂_θ - i / sin(θ) ∂_ϕ) (sin θ)^s F
# (Calculated manually from sYlm above)
function ð̄sYlm(s::Integer, l::Integer, m::Integer, θ::Real, ϕ::Real)
    @assert abs(s) ≤ l
    @assert -l ≤ m ≤ l

    (s, l, m) == (0, 0, 0) && return 0
    (s, l, m) == (0, 1, -1) && return sqrt(3 / 8π) * (1 - cos(θ)) * cis(-ϕ)
    (s, l, m) == (0, 1, 0) && return sqrt(3 / 4π) * sin(θ)
    (s, l, m) == (0, 1, +1) && return sqrt(3 / 8π) * (1 + cos(θ)) * cis(ϕ)
    (s, l, m) == (0, 1, +1) && return -sqrt(3 / 8π) * (1 - cos(θ)) * cis(ϕ)
    (s, l, m) == (0, 2, -2) && return -sqrt(15 / 8π) * (-1 + cos(θ)) * sin(θ) * cis(-2ϕ)
    (s, l, m) == (0, 2, -1) && return sqrt(15 / 2π) * (1 + 2 * cos(θ)) * sin(θ / 2)^2 * cis(-ϕ)
    (s, l, m) == (0, 2, 0) && return sqrt(45 / 16π) * sin(2θ)
    (s, l, m) == (0, 2, +1) && return sqrt(15 / 2π) * cos(θ / 2)^2 * (-1 + 2 * cos(θ)) * cis(ϕ)
    (s, l, m) == (0, 2, +2) && return -sqrt(15 / 8π) * (1 + cos(θ)) * sin(θ) * cis(2ϕ)

    (s, l, m) == (+1, 1, -1) && return -sqrt(3 / 4π) * sin(θ) * cis(-ϕ)
    (s, l, m) == (+1, 1, 0) && return -sqrt(3 / 2π) * cos(θ)
    (s, l, m) == (+1, 1, +1) && return sqrt(3 / 4π) * sin(θ) * cis(ϕ)
    (s, l, m) == (+1, 2, -2) && return -sqrt(45 / 16π) * sin(θ)^2 * cis(-2ϕ)
    (s, l, m) == (+1, 2, -1) && return -sqrt(45 / 16π) * sin(2θ) * cis(-ϕ)
    (s, l, m) == (+1, 2, 0) && return -sqrt(15 / 32π) * (1 + 3 * cos(2θ))
    (s, l, m) == (+1, 2, +1) && return sqrt(45 / 16π) * sin(2θ) * cis(ϕ)
    (s, l, m) == (+1, 2, +2) && return -sqrt(45 / 16π) * sin(θ)^2 * cis(2ϕ)

    (s, l, m) == (+2, 2, -2) && return sqrt(5 / 4π) * (1 + cos(θ)) * sin(θ) * cis(-2ϕ)
    (s, l, m) == (+2, 2, -1) && return sqrt(5 / 4π) * (cos(θ) + cos(2θ)) * cis(-ϕ)
    (s, l, m) == (+2, 2, 0) && return -sqrt(15 / 2π) * cos(θ) * sin(θ)
    (s, l, m) == (+2, 2, +1) && return sqrt(5 / 4π) * (cos(θ) - cos(2θ)) * cis(ϕ)
    (s, l, m) == (+2, 2, +2) && return -sqrt(5 / π) * sin(θ / 2)^2 * sin(θ) * cis(2ϕ)

    @assert false
end

################################################################################

@testset "Grid point counts (Driscoll & Healy)" begin
    for L in 1:10
        n = ssht.sampling_dh_n(L)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)
        @test nphi == 2 * L - 1
        @test ntheta == 2 * L
        @test n == nphi * ntheta
    end
end

@testset "Coordinates (Driscoll & Healy)" begin
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

@testset "Grid point counts (McEwen & Wiaux)" begin
    for L in 1:2:10
        n = ssht.sampling_mw_n(L)
        nphi = ssht.sampling_mw_nphi(L)
        ntheta = ssht.sampling_mw_ntheta(L)
        @test nphi == 2 * L - 1
        @test ntheta == L
        # The south pole is covered only by a single point
        @test n == nphi * (ntheta - 1) + 1
    end
end

@testset "Coordinates (McEwen & Wiaux)" begin
    for L in 1:100
        nphi = ssht.sampling_mw_nphi(L)
        ntheta = ssht.sampling_mw_ntheta(L)
        for p in 1:nphi
            phi = ssht.sampling_mw_p2phi(p, L)
            @test phi ≈ 2π * (p - 1) / nphi
        end
        for t in 1:ntheta
            theta = ssht.sampling_mw_t2theta(t, L)
            @test theta ≈ π * (2 * t - 1) / (2 * ntheta - 1)
        end
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

################################################################################

# TODO: Test phase and parity of sYlm
#
# {\displaystyle Y_{\ell }^{m}(\theta ,\phi )\to Y_{\ell }^{m}(\pi -\theta ,\pi +\phi )=(-1)^{\ell }Y_{\ell }^{m}(\theta ,\phi )}
#
# <math display="block">Y_\ell^m(\theta,\phi) \to Y_\ell^m(\pi-\theta,\pi+\phi) = (-1)^\ell Y_\ell^m(\theta,\phi)</math>
#
# {}_s\bar Y_{l m} &= \left(-1\right)^{s+m}{}_{-s}Y_{l(-m)}\\
#
# {}_sY_{l m}(\pi-\theta,\phi+\pi) &= \left(-1\right)^l {}_{-s}Y_{l m}(\theta,\phi).

Random.seed!(100)
modes = [(name="(l=0,m=0)", fun=(θ, ϕ) -> sYlm(0, 0, 0, θ, ϕ), modes=[1, 0, 0, 0]),
         (name="(l=1,m=-1)", fun=(θ, ϕ) -> im * sYlm(0, 1, +1, θ, ϕ) + im * sYlm(0, 1, -1, θ, ϕ), modes=[0, im, 0, im]),
         (name="(l=1,m=0)", fun=(θ, ϕ) -> sYlm(0, 1, 0, θ, ϕ), modes=[0, 0, 1, 0]),
         (name="(l=1,m=+1)", fun=(θ, ϕ) -> sYlm(0, 1, +1, θ, ϕ) - sYlm(0, 1, -1, θ, ϕ), modes=[0, -1, 0, 1])]
@testset "Simple real transforms: $(mode.name)" for mode in modes
    for L in 2:20
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        f = Array{Float64}(undef, nphi, ntheta)
        for p in 1:nphi, t in 1:ntheta
            ϕ = ssht.sampling_dh_p2phi(p, L)
            θ = ssht.sampling_dh_t2theta(t, L)
            f[p, t] = mode.fun(θ, ϕ)
        end

        flm = ssht.core_dh_forward_sov_real(f, L)
        @assert L ≥ 2
        @test flm[1:4] ≈ mode.modes
        if L > 2
            @test all(isapprox(0; atol=100eps()), flm[5:end])
        end

        f′ = ssht.core_dh_inverse_sov_real(flm, L)
        @test isapprox(f′, f; atol=1000eps())
    end
end

Random.seed!(100)
modes = [(name="(s=$s,l=$l,m=$m)", spin=s, el=l, fun=(θ, ϕ) -> sYlm(s, l, m, θ, ϕ),
          modes=L -> [l == l′ && m == m′ for l′ in 0:(L - 1) for m′ in (-l′):l′]) for s in -2:+2 for l in abs(s):2 for m in (-l):l]
@testset "Simple complex transforms: $(mode.name)" for mode in modes
    for L in (mode.el + 1):20
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        f = Array{Complex{Float64}}(undef, nphi, ntheta)
        for p in 1:nphi, t in 1:ntheta
            ϕ = ssht.sampling_dh_p2phi(p, L)
            θ = ssht.sampling_dh_t2theta(t, L)
            f[p, t] = mode.fun(θ, ϕ)
        end

        flm = ssht.core_dh_forward_sov(f, L, mode.spin)
        @test flm ≈ mode.modes(L)

        f′ = ssht.core_dh_inverse_sov(flm, L, mode.spin)
        @test isapprox(f′, f; atol=1000eps())
    end
end

Random.seed!(100)
@testset "Linearity of real transforms" begin
    for iter in 1:100
        L = rand(1:100)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        f = randn(Float64, nphi, ntheta)
        g = randn(Float64, nphi, ntheta)
        α = randn(Float64)
        h = f + α * g

        flm = ssht.core_dh_forward_sov_real(f, L)
        glm = ssht.core_dh_forward_sov_real(g, L)
        hlm = ssht.core_dh_forward_sov_real(h, L)

        @test flm + α * glm ≈ hlm
    end
end

Random.seed!(100)
@testset "Linearity of complex transforms" begin
    for iter in 1:100
        L = rand(1:100)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        spin = rand(-4:4)

        f = randn(Complex{Float64}, nphi, ntheta)
        g = randn(Complex{Float64}, nphi, ntheta)
        α = randn(Complex{Float64})
        h = f + α * g

        flm = ssht.core_dh_forward_sov(f, L, spin)
        glm = ssht.core_dh_forward_sov(g, L, spin)
        hlm = ssht.core_dh_forward_sov(h, L, spin)

        @test flm + α * glm ≈ hlm

        flm = randn(Complex{Float64}, L^2)
        glm = randn(Complex{Float64}, L^2)
        hlm = flm + α * glm

        f = ssht.core_dh_inverse_sov(flm, L, spin)
        g = ssht.core_dh_inverse_sov(glm, L, spin)
        h = ssht.core_dh_inverse_sov(hlm, L, spin)

        @test f + α * g ≈ h
    end
end

Random.seed!(100)
@testset "Orthonormality of complex transforms" begin
    for iter in 1:100
        L = rand(1:100)
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        smax = min(4, L - 1)
        spin = rand((-smax):smax)

        flm = zeros(Complex{Float64}, L^2)
        glm = zeros(Complex{Float64}, L^2)

        lf = rand(abs(spin):(L - 1))
        mf = rand((-lf):lf)
        lg = rand(abs(spin):(L - 1))
        mg = rand((-lg):lg)

        flm[ssht.sampling_elm2ind(lf, mf)] = 1
        glm[ssht.sampling_elm2ind(lg, mg)] = 1

        f = ssht.core_dh_inverse_sov(flm, L, spin)
        g = ssht.core_dh_inverse_sov(glm, L, spin)

        @test isapprox(integrate(f, f, L), 1; atol=1 / L^2)
        @test isapprox(integrate(f, g, L), (lf == lg) * (mf == mg); atol=1 / L^2)

        h = conj(f) .* f
        hlm = ssht.core_dh_forward_sov(h, L, 0)
        @test isapprox(hlm[ssht.sampling_elm2ind(0, 0)], sqrt(1 / 4π); atol=sqrt(eps()))

        h = conj(f) .* g
        hlm = ssht.core_dh_forward_sov(h, L, 0)
        @test isapprox(hlm[ssht.sampling_elm2ind(0, 0)], (lf == lg) * (mf == mg) * sqrt(1 / 4π); atol=sqrt(eps()))
    end
end

Random.seed!(100)
modes = [(name="(s=$s,l=$l,m=$m)", spin=s, el=l, fun=(θ, ϕ) -> sYlm(s, l, m, θ, ϕ), ðfun=(θ, ϕ) -> ðsYlm(s, l, m, θ, ϕ),
          ð̄fun=(θ, ϕ) -> ð̄sYlm(s, l, m, θ, ϕ)) for s in 0:+2 for l in abs(s):2 for m in (-l):l]
@testset "Simple derivatives (eth, eth-bar): $(mode.name)" for mode in modes
    for L in (mode.el + 1):20
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        f = Array{Complex{Float64}}(undef, nphi, ntheta)
        ðf₀ = Array{Complex{Float64}}(undef, nphi, ntheta)
        ð̄f₀ = Array{Complex{Float64}}(undef, nphi, ntheta)
        for p in 1:nphi, t in 1:ntheta
            ϕ = ssht.sampling_dh_p2phi(p, L)
            θ = ssht.sampling_dh_t2theta(t, L)
            f[p, t] = mode.fun(θ, ϕ)
            ðf₀[p, t] = mode.ðfun(θ, ϕ)
            ð̄f₀[p, t] = mode.ð̄fun(θ, ϕ)
        end

        flm = ssht.core_dh_forward_sov(f, L, mode.spin)

        ðflm = ssht.eth(flm, L, mode.spin)
        ðf = ssht.core_dh_inverse_sov(ðflm, L, mode.spin + 1)
        @test isapprox(ðf, ðf₀; atol=10000eps())

        ð̄flm = ssht.ethbar(flm, L, mode.spin)
        ð̄f = ssht.core_dh_inverse_sov(ð̄flm, L, mode.spin - 1)
        @test isapprox(ð̄f, ð̄f₀; atol=10000eps())
    end
end
