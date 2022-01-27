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

# These Ylm are taken from Wikipedia:
# <https://en.wikipedia.org/wiki/Spherical_harmonics> and
# <https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics>
function sYlm(s::Integer, l::Integer, m::Integer, θ::Real, ϕ::Real)
    @assert abs(s) ≤ l
    @assert -l ≤ m ≤ l

    # Parity:
    s < 0 && return bitsign(s + m) * conj(sYlm(-s, l, -m, θ, ϕ))

    (s, l, m) == (0, 0, 0) && return sqrt(1 / 4π)
    (s, l, m) == (0, 1, -1) && return sqrt(3 / 8π) * sin(θ) * cis(-ϕ)
    (s, l, m) == (0, 1, 0) && return sqrt(3 / 4π) * cos(θ)
    (s, l, m) == (0, 1, +1) && return -sqrt(3 / 8π) * sin(θ) * cis(ϕ)

    (s, l, m) == (+1, 1, -1) && return -sqrt(3 / 16π) * (1 + cos(θ)) * cis(-ϕ)
    (s, l, m) == (+1, 1, 0) && return sqrt(3 / 8π) * sin(θ)
    (s, l, m) == (+1, 1, +1) && return -sqrt(3 / 16π) * (1 - cos(θ)) * cis(ϕ)
    (s, l, m) == (+1, 2, -2) && return -sqrt(5 / 16π) * (1 + cos(θ)) * sin(θ) * cis(-2ϕ)
    (s, l, m) == (+1, 2, -1) && return -sqrt(5 / 16π) * (cos(θ) + cos(2θ)) * cis(-ϕ)
    (s, l, m) == (+1, 2, 0) && return sqrt(15 / 8π) * cos(θ) * sin(θ)
    (s, l, m) == (+1, 2, +1) && return -sqrt(5 / 16π) * (-cos(θ) + cos(2θ)) * cis(ϕ)
    (s, l, m) == (+1, 2, +2) && return -sqrt(5 / 16π) * (-1 + cos(θ)) * sin(θ) * cis(2ϕ)

    # (s, l, m) == (+2, 2, -2)

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
    (s, l, m) == (+1, 1, -1) && return 0
    (s, l, m) == (+1, 1, 0) && return 0
    (s, l, m) == (+1, 1, +1) && return 0

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
    (s, l, m) == (+1, 1, -1) && return 0
    (s, l, m) == (+1, 1, 0) && return 0
    (s, l, m) == (+1, 1, +1) && return 0

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
modes = [(name="(0,0)", fun=(θ, ϕ) -> sYlm(0, 0, 0, θ, ϕ), modes=[1, 0, 0, 0]),
         (name="(1,-1)", fun=(θ, ϕ) -> im * sYlm(0, 1, +1, θ, ϕ) + im * sYlm(0, 1, -1, θ, ϕ), modes=[0, im, 0, im]),
         (name="(1,0)", fun=(θ, ϕ) -> sYlm(0, 1, 0, θ, ϕ), modes=[0, 0, 1, 0]),
         (name="(1,+1)", fun=(θ, ϕ) -> sYlm(0, 1, +1, θ, ϕ) - sYlm(0, 1, -1, θ, ϕ), modes=[0, -1, 0, 1])]

@testset "Simple real transforms: $(mode.name)" for mode in modes
    verbosity = 0
    for L in 2:20
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        f = Array{Float64}(undef, nphi, ntheta)
        for p in 1:nphi, t in 1:ntheta
            ϕ = ssht.sampling_dh_p2phi(p, L)
            θ = ssht.sampling_dh_t2theta(t, L)
            f[p, t] = mode.fun(θ, ϕ)
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
modes = [(name="(0,0,0)", spin=0, fun=(θ, ϕ) -> sYlm(0, 0, 0, θ, ϕ), modes=[1, 0, 0, 0]),
         (name="(0,0,-1)", spin=0, fun=(θ, ϕ) -> sYlm(0, 1, -1, θ, ϕ), modes=[0, 1, 0, 0]),
         (name="(0,1,0)", spin=0, fun=(θ, ϕ) -> sYlm(0, 1, 0, θ, ϕ), modes=[0, 0, 1, 0]),
         (name="(0,0,+1)", spin=0, fun=(θ, ϕ) -> sYlm(0, 1, +1, θ, ϕ), modes=[0, 0, 0, 1]),
         (name="(+1,1,-1)", spin=+1, fun=(θ, ϕ) -> sYlm(+1, 1, -1, θ, ϕ), modes=[0, 1, 0, 0]),
         (name="(+1,1,0)", spin=+1, fun=(θ, ϕ) -> sYlm(+1, 1, 0, θ, ϕ), modes=[0, 0, 1, 0]),
         (name="(+1,1,+1)", spin=+1, fun=(θ, ϕ) -> sYlm(+1, 1, +1, θ, ϕ), modes=[0, 0, 0, 1]),
         (name="(-1,1,-1)", spin=-1, fun=(θ, ϕ) -> sYlm(-1, 1, -1, θ, ϕ), modes=[0, 1, 0, 0]),
         (name="(-1,1,0)", spin=-1, fun=(θ, ϕ) -> sYlm(-1, 1, 0, θ, ϕ), modes=[0, 0, 1, 0]),
         (name="(-1,1,+1)", spin=-1, fun=(θ, ϕ) -> sYlm(-1, 1, +1, θ, ϕ), modes=[0, 0, 0, 1])]
@testset "Simple complex transforms: $(mode.name) spin=$(mode.spin)" for mode in modes
    verbosity = 0
    for L in 2:20
        nphi = ssht.sampling_dh_nphi(L)
        ntheta = ssht.sampling_dh_ntheta(L)

        f = Array{Complex{Float64}}(undef, nphi, ntheta)
        for p in 1:nphi, t in 1:ntheta
            ϕ = ssht.sampling_dh_p2phi(p, L)
            θ = ssht.sampling_dh_t2theta(t, L)
            f[p, t] = mode.fun(θ, ϕ)
        end

        flm = Array{Complex{Float64}}(undef, L^2)
        ssht.core_dh_forward_sov!(flm, f, L, mode.spin, verbosity)
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

Random.seed!(100)
modes = [(name="(0,0,0)", spin=0, fun=(θ, ϕ) -> sYlm(0, 0, 0, θ, ϕ), ðfun=(θ, ϕ) -> ðsYlm(0, 0, 0, θ, ϕ),
          ð̄fun=(θ, ϕ) -> ð̄sYlm(0, 0, 0, θ, ϕ)),
         (name="(0,1,-1)", spin=0, fun=(θ, ϕ) -> sYlm(0, 1, -1, θ, ϕ), ðfun=(θ, ϕ) -> ðsYlm(0, 1, -1, θ, ϕ),
          ð̄fun=(θ, ϕ) -> ð̄sYlm(0, 1, -1, θ, ϕ)),
         (name="(0,1,0)", spin=0, fun=(θ, ϕ) -> sYlm(0, 1, 0, θ, ϕ), ðfun=(θ, ϕ) -> ðsYlm(0, 1, 0, θ, ϕ),
          ð̄fun=(θ, ϕ) -> ð̄sYlm(0, 1, 0, θ, ϕ)),
         (name="(0,1,+1)", spin=0, fun=(θ, ϕ) -> sYlm(0, 1, +1, θ, ϕ), ðfun=(θ, ϕ) -> ðsYlm(0, 1, +1, θ, ϕ),
          ð̄fun=(θ, ϕ) -> ð̄sYlm(0, 1, +1, θ, ϕ)),
         (name="(+1,1,-1)", spin=1, fun=(θ, ϕ) -> sYlm(+1, 1, -1, θ, ϕ), ðfun=(θ, ϕ) -> ðsYlm(+1, 1, -1, θ, ϕ),
          ð̄fun=(θ, ϕ) -> ð̄sYlm(+1, 1, -1, θ, ϕ)),
         (name="(+1,1,0)", spin=1, fun=(θ, ϕ) -> sYlm(+1, 1, 0, θ, ϕ), ðfun=(θ, ϕ) -> ðsYlm(+1, 1, 0, θ, ϕ),
          ð̄fun=(θ, ϕ) -> ð̄sYlm(+1, 1, 0, θ, ϕ)),
         (name="(+1,1,+1)", spin=1, fun=(θ, ϕ) -> sYlm(+1, 1, +1, θ, ϕ), ðfun=(θ, ϕ) -> ðsYlm(+1, 1, +1, θ, ϕ),
          ð̄fun=(θ, ϕ) -> ð̄sYlm(+1, 1, +1, θ, ϕ))]
@testset "Simple derivatives (eth, eth-bar): $(mode.name)" for mode in modes
    verbosity = 0
    for L in 2:20
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

        flm = ssht.core_dh_forward_sov(f, L, mode.spin, verbosity)

        ðflm = ssht.eth(flm, L, mode.spin)
        ðf = ssht.core_dh_inverse_sov(ðflm, L, mode.spin + 1, verbosity)
        @test isapprox(ðf, ðf₀; atol=10000eps())

        ð̄flm = ssht.ethbar(flm, L, mode.spin)
        ð̄f = ssht.core_dh_inverse_sov(ð̄flm, L, mode.spin - 1, verbosity)
        @test isapprox(ð̄f, ð̄f₀; atol=10000eps())
    end
end
