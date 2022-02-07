module SSHT

using ssht_jll

# File ssht_adjoint.h

# File ssht_core.h

# L is number of modes; 0 ≤ l < L
# verbosity: 0 ≤ v ≤ 5

"""
    SSHT.core_dh_inverse_sov!(f::AbstractArray{Complex{Float64},2},
                              flm::AbstractVector{Complex{Float64}},
                              L::Integer,
                              spin::Integer,
                              verbosity::Integer=0)

Evaluate spin-weighted spherical harmonic coefficients `flm` with
spin-weight `spin`. `L = lmax+1` is the number of modes in `flm`.

The arrays `flm` and `f` must have the following sizes:

    nphi = SSHT.sampling_dh_nphi(L)
    ntheta = SSHT.sampling_dh_ntheta(L)
    length(flm) == L^2
    size(f) == (nphi, ntheta)

See the [SSHT
reference](https://astro-informatics.github.io/ssht/c/html/ssht__core_8h.html)
for details.

See also: [`SSHT.core_dh_inverse_sov`](@ref),
[`SSHT.core_dh_inverse_sov_real!`](@ref),
[`SSHT.core_dh_forward_sov!`](@ref).
"""
function core_dh_inverse_sov!(f::AbstractArray{Complex{Float64},2}, flm::AbstractVector{Complex{Float64}}, L::Integer,
                              spin::Integer, verbosity::Integer=0)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    any(size(f) .< (nphi, ntheta)) && throw(DimensionMismatch())
    length(flm) < L^2 && throw(DimensionMismatch())
    ccall((:ssht_core_dh_inverse_sov, libssht), Cvoid, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Cint, Cint, Cint), f, flm, L,
          spin, verbosity)
    return f
end

"""
    f = SSHT.core_dh_inverse_sov(flm::AbstractVector{Complex{Float64}},
                                 L::Integer,
                                 spin::Integer,
                                 verbosity::Integer=0)
    f::Array{Complex{Float64},2}

Evaluate spin-weighted spherical harmonic coefficients `flm` with
spin-weight `spin`. `L = lmax+1` is the number of modes in `flm`.

The array `flm` must have the length `L^2`. The result `f` will have
the following size:

    nphi = SSHT.sampling_dh_nphi(L)
    ntheta = SSHT.sampling_dh_ntheta(L)
    size(f) == (nphi, ntheta)

See the [SSHT
reference](https://astro-informatics.github.io/ssht/c/html/ssht__core_8h.html)
for details.

See also: [`SSHT.core_dh_inverse_sov!`](@ref),
[`SSHT.core_dh_inverse_sov_real`](@ref),
[`SSHT.core_dh_forward_sov`](@ref).
"""
function core_dh_inverse_sov(flm::AbstractVector{Complex{Float64}}, L::Integer, spin::Integer, verbosity::Integer=0)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    f = Array{Complex{Float64}}(undef, nphi, ntheta)
    core_dh_inverse_sov!(f, flm, L, spin, verbosity)
    return f
end

"""
    SSHT.core_dh_inverse_sov_real!(f::AbstractArray{Float64,2},
                                   flm::AbstractVector{Complex{Float64}},
                                   L::Integer,
                                   spin::Integer,
                                   verbosity::Integer=0)

Evaluate spin-weighted spherical harmonic coefficients `flm` with
spin-weight `spin`. `L = lmax+1` is the number of modes in `flm`. The
result `f` is real-valued.

The arrays `flm` and `f` must have the following sizes:

    nphi = SSHT.sampling_dh_nphi(L)
    ntheta = SSHT.sampling_dh_ntheta(L)
    length(flm) == L^2
    size(f) == (nphi, ntheta)

See the [SSHT
reference](https://astro-informatics.github.io/ssht/c/html/ssht__core_8h.html)
for details.

See also: [`SSHT.core_dh_inverse_sov_real`](@ref),
[`SSHT.core_dh_inverse_sov!`](@ref),
[`SSHT.core_dh_forward_sov_real!`](@ref).
"""
function core_dh_inverse_sov_real!(f::AbstractArray{Float64,2}, flm::AbstractVector{Complex{Float64}}, L::Integer,
                                   verbosity::Integer=0)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    any(size(f) .< (nphi, ntheta)) && throw(DimensionMismatch())
    length(flm) < L^2 && throw(DimensionMismatch())
    ccall((:ssht_core_dh_inverse_sov_real, libssht), Cvoid, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Cint, Cint), f, flm, L,
          verbosity)
    return f
end

"""
    f = SSHT.core_dh_inverse_sov_real(flm::AbstractVector{Complex{Float64}},
                                      L::Integer,
                                      spin::Integer,
                                      verbosity::Integer=0)
    f::Array{Float64,2}

Evaluate spin-weighted spherical harmonic coefficients `flm` with
spin-weight `spin`. `L = lmax+1` is the number of modes in `flm`. The
result `f` is real-valued.

The array `flm` must have the length `L^2`. The result `f` will have
the following size:

    nphi = SSHT.sampling_dh_nphi(L)
    ntheta = SSHT.sampling_dh_ntheta(L)
    size(f) == (nphi, ntheta)

See the [SSHT
reference](https://astro-informatics.github.io/ssht/c/html/ssht__core_8h.html)
for details.

See also: [`SSHT.core_dh_inverse_sov_real!`](@ref),
[`SSHT.core_dh_inverse_sov`](@ref),
[`SSHT.core_dh_forward_sov_real`](@ref).
"""
function core_dh_inverse_sov_real(flm::AbstractVector{Complex{Float64}}, L::Integer, verbosity::Integer=0)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    f = Array{Float64}(undef, nphi, ntheta)
    core_dh_inverse_sov_real!(f, flm, L, verbosity)
    return f
end

"""
    SSHT.core_dh_forward_sov!(flm::AbstractVector{Complex{Float64}},
                              f::AbstractArray{Complex{Float64},2},
                              L::Integer,
                              spin::Integer,
                              verbosity::Integer=0)

Calculate spin-weighted spherical harmonic coefficients with
spin-weight `spin` `flm` from grid point values `f`. `L = lmax+1` is
the number of modes in `flm`.

The arrays `flm` and `f` must have the following sizes:

    nphi = SSHT.sampling_dh_nphi(L)
    ntheta = SSHT.sampling_dh_ntheta(L)
    length(flm) == L^2
    size(f) == (nphi, ntheta)

Use [`SSTH.sampling_elm2ind`](@ref) to access individual modes in the
`flm` array.

See the [SSHT
reference](https://astro-informatics.github.io/ssht/c/html/ssht__core_8h.html)
for details.

See also: [`SSHT.core_dh_forward_sov`](@ref),
[`SSHT.core_dh_forward_sov_real!`](@ref),
[`SSHT.core_dh_inverse_sov!`](@ref).
"""
function core_dh_forward_sov!(flm::AbstractVector{Complex{Float64}}, f::AbstractArray{Complex{Float64},2}, L::Integer,
                              spin::Integer, verbosity::Integer=0)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    length(flm) < L^2 && throw(DimensionMismatch())
    any(size(f) .< (nphi, ntheta)) && throw(DimensionMismatch())
    ccall((:ssht_core_dh_forward_sov, libssht), Cvoid, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Cint, Cint, Cint), flm, f, L,
          spin, verbosity)
    return flm
end

"""
    flm = SSHT.core_dh_forward_sov(fl:AbstractVector{Complex{Float64}},
                                   L::Integer,
                                   spin::Integer,
                                   verbosity::Integer=0)
    flm::Vector{Complex{Float64}}

Calculate spin-weighted spherical harmonic coefficients `flm` with
spin-weight `spin` from grid point values `f`. `L = lmax+1` is the
number of modes in `flm`.

The array `flm` will have length `L^2`. The input `f` must have the
following size:

    nphi = SSHT.sampling_dh_nphi(L)
    ntheta = SSHT.sampling_dh_ntheta(L)
    size(f) == (nphi, ntheta)

Use [`SSTH.sampling_elm2ind`](@ref) to access individual modes in the
`flm` array.

See the [SSHT
reference](https://astro-informatics.github.io/ssht/c/html/ssht__core_8h.html)
for details.

See also: [`SSHT.core_dh_forward_sov!`](@ref),
[`SSHT.core_dh_forward_sov_real`](@ref),
[`SSHT.core_dh_inverse_sov`](@ref).
"""
function core_dh_forward_sov(f::AbstractArray{Complex{Float64},2}, L::Integer, spin::Integer, verbosity::Integer=0)
    flm = Array{Complex{Float64}}(undef, L^2)
    core_dh_forward_sov!(flm, f, L, spin, verbosity)
    return flm
end

"""
    SSHT.core_dh_forward_sov_real!(flm::AbstractVector{Complex{Float64}},
                                   f::AbstractArray{Float64,2},
                                   L::Integer,
                                   spin::Integer,
                                   verbosity::Integer=0)

Calculate spin-weighted spherical harmonic coefficients `flm` with
spin-weight `spin` from real-valued grid point values `f`. `L =
lmax+1` is the number of modes in `flm`.

The arrays `flm` and `f` must have the following sizes:

    nphi = SSHT.sampling_dh_nphi(L)
    ntheta = SSHT.sampling_dh_ntheta(L)
    length(flm) == L^2
    size(f) == (nphi, ntheta)

Use [`SSTH.sampling_elm2ind`](@ref) to access individual modes in the
`flm` array.

See the [SSHT
reference](https://astro-informatics.github.io/ssht/c/html/ssht__core_8h.html)
for details.

See also: [`SSHT.core_dh_forward_sov_real`](@ref),
[`SSHT.core_dh_forward_sov!`](@ref),
[`SSHT.core_dh_inverse_sov_real!`](@ref).
"""
function core_dh_forward_sov_real!(flm::AbstractVector{Complex{Float64}}, f::AbstractArray{Float64,2}, L::Integer,
                                   verbosity::Integer=0)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    length(flm) < L^2 && throw(DimensionMismatch())
    any(size(f) .< (nphi, ntheta)) && throw(DimensionMismatch())
    ccall((:ssht_core_dh_forward_sov_real, libssht), Cvoid, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Cint, Cint), flm, f, L,
          verbosity)
    return flm
end

"""
    flm = SSHT.core_dh_inverse_sov_real(f::AbstractVector{Float64},
                                        L::Integer,
                                        spin::Integer,
                                        verbosity::Integer=0)
    f::Array{Complex{Float64},2}

Calculate spin-weighted spherical harmonic coefficients `flm` with
spin-weight `spin` from real-valued grid point values `f`. `L =
lmax+1` is the number of modes in `flm`.

The array `flm` must have the length `L^2`. The result `f` will have
the following size:

    nphi = SSHT.sampling_dh_nphi(L)
    ntheta = SSHT.sampling_dh_ntheta(L)
    size(f) == (nphi, ntheta)

Use [`SSTH.sampling_elm2ind`](@ref) to access individual modes in the
`flm` array.

See the [SSHT
reference](https://astro-informatics.github.io/ssht/c/html/ssht__core_8h.html)
for details.

See also: [`SSHT.core_dh_forward_sov_real!`](@ref),
[`SSHT.core_dh_forward_sov`](@ref),
[`SSHT.core_dh_inverse_sov_real`](@ref).
"""
function core_dh_forward_sov_real(f::AbstractArray{Float64,2}, L::Integer, verbosity::Integer=0)
    flm = Array{Complex{Float64}}(undef, L^2)
    core_dh_forward_sov_real!(flm, f, L, verbosity)
    return flm
end

# File ssht_dl.h

@enum dl_size_t begin
    DL_QUARTER = 0
    SSHT_DL_QUARTER_EXTENDED
    SSHT_DL_HALF
    SSHT_DL_FULL
end

@enum dl_method_t begin
    SSHT_DL_RISBO = 0
    SSHT_DL_TRAPANI
end

# File ssht_sampling.h

"""
    n = sampling_dh_n(L::Integer)
    n::Int

Total number of unique collocation points for `L = lmax+1` modes that.
(If the collocation scheme places points onto the poles, then not all
collocation points are unique.)

See also: [`SSHT.sampling_dh_nphi`](@ref),
[`SSHT.sampling_dh_ntheta`](@ref).
"""
sampling_dh_n(L::Integer) = Int(ccall((:ssht_sampling_dh_n, libssht), Cint, (Cint,), L))

"""
    nphi = sampling_dh_nphi(L::Integer)
    nphi::Int

Number of collocation points in the `phi` direction for `L = lmax+1`
modes.

See also: [`SSHT.sampling_dh_n`](@ref),
[`SSHT.sampling_dh_ntheta`](@ref).
"""
sampling_dh_nphi(L::Integer) = Int(ccall((:ssht_sampling_dh_nphi, libssht), Cint, (Cint,), L))

"""
    ntheta = sampling_dh_ntheta(L::Integer)
    ntheta::Int

Number of collocation points in the `theta` direction for `L = lmax+1`
modes.

See also: [`SSHT.sampling_dh_n`](@ref),
[`SSHT.sampling_dh_nphi`](@ref).
"""
sampling_dh_ntheta(L::Integer) = Int(ccall((:ssht_sampling_dh_ntheta, libssht), Cint, (Cint,), L))

"""
    phi = sampling_dh_p2phi(p::Integer, L::Integer)
    phi::Float64

Calculate the `phi` coordinate for point `p` in the `phi` direction
(`1 ≤ p ≤ nphi`).

See also: [`SSHT.sampling_dh_t2theta`](@ref),
[`SSHT.sampling_dh_nphi`](@ref).
"""
sampling_dh_p2phi(p::Integer, L::Integer) = Float64(ccall((:ssht_sampling_dh_p2phi, libssht), Cdouble, (Cint, Cint), p - 1, L))

"""
    theta = sampling_dh_p2theta(t::Integer, L::Integer)
    theta::Float64

Calculate the `theta` coordinate for point `t` in the `theta` direction
(`1 ≤ t ≤ ntheta`).

See also: [`SSHT.sampling_dh_p2phi`](@ref),
[`SSHT.sampling_dh_ntheta`](@ref).
"""
sampling_dh_t2theta(t::Integer, L::Integer) = Float64(ccall((:ssht_sampling_dh_t2theta, libssht), Cdouble, (Cint, Cint), t - 1, L))

"""
    w = sampling_weight_dh(theta_t::Real, L::Integer)
    w::Float64

Calculate the sampling weight of a point with the given `theta`
coordinate. This value is essentially `sin(theta) * dtheta`, possibly
modified to achive higher accuracy depending on how the collocation
points are spaced out in the `theta` direction.

Integrating a function over the sphere should be done as follows:

    # Input: Choose `L` and array `f`
    s = 0.0
    nphi = SSHT.sampling_dh_nphi(L)
    ntheta = SSHT.sampling_dh_ntheta(L)
    for p in 1:nphi, t in 1:ntheta
        phi = SSHT.sampling_dh_p2phi(p, L)
        theta = SSHT.sampling_dh_t2theta(t, L)
        # dtheta = π / ntheta
        dtheta = SSHT.sampling_weight_dh(theta, L)
        dphi = 2π / nphi
        s += f[p, t] * dtheta * dphi
    end
    # Output: s

See also: [`SSHT.sampling_dh_nphi`](@ref),
[`SSHT.sampling_dh_ntheta`](@ref), [`SSHT.sampling_dh_p2phi`](@ref),
[`SSHT.sampling_dh_t2theta`](@ref).
"""
function sampling_weight_dh(theta_t::Real, L::Integer)
    return Float64(ccall((:ssht_sampling_weight_dh, libssht), Cdouble, (Cdouble, Cint), theta_t, L))
end

sampling_mw_n(L::Integer) = Int(ccall((:ssht_sampling_mw_n, libssht), Cint, (Cint,), L))
sampling_mw_nphi(L::Integer) = Int(ccall((:ssht_sampling_mw_nphi, libssht), Cint, (Cint,), L))
sampling_mw_ntheta(L::Integer) = Int(ccall((:ssht_sampling_mw_ntheta, libssht), Cint, (Cint,), L))
sampling_mw_p2phi(p::Integer, L::Integer) = Float64(ccall((:ssht_sampling_mw_p2phi, libssht), Cdouble, (Cint, Cint), p - 1, L))
sampling_mw_t2theta(t::Integer, L::Integer) = Float64(ccall((:ssht_sampling_mw_t2theta, libssht), Cdouble, (Cint, Cint), t - 1, L))
function sampling_weight_mw(theta_t::Real, L::Integer)
    return Float64(ccall((:ssht_sampling_weight_mw, libssht), Cdouble, (Cdouble, Cint), theta_t, L))
end

"""
    ind = sampling_elm2ind(el::Integer, m::Integer)
    ind::Int

Calculate the mode array index `ind` for a given mode `l`, `m`. For `L
= lmax+1` modes, there are `L^2` modes in total.

See also: [`sampling_ind2elm`](@ref).
"""
function sampling_elm2ind(el::Int, m::Int)
    0 ≤ el || throw(DomainError(el, "Need 0 ≤ el"))
    -el ≤ m ≤ el || throw(DomainError(m, "Need -el ≤ m ≤ el"))
    return el^2 + el + m + 1
end
sampling_elm2ind(el::Integer, m::Integer) = sampling_elm2ind(Int(el), Int(m))

"""
    l, m = sampling_ind2elm(ind::Integer)
    l::Int
    m::Int

Calculate the mode numbers `l` and `m` from a given mode array index
`ind`. For `L = lmax+1` modes, there are `L^2` modes in total: `1 ≤
ind ≤ L^2`.

This function needs to evaluate a square root internally. If possible,
using [`sampling_elm2ind`](@ref) instead is slightly preferred.

See also: [`sampling_elm2ind`](@ref).
"""
function sampling_ind2elm(ind::Int)
    1 ≤ ind || throw(DomainError(ind, "Need 1 ≤ ind"))
    el = isqrt(ind - 1)
    m = ind - el^2 - el - 1
    return el, m
end
sampling_ind2elm(ind::Integer) = sampling_ind2elm(Int(ind))

################################################################################

"""
    eth!(ðflm::AbstractVector, flm::AbstractVector, L::Integer, spin::Integer)

Calculate the `ð` (eth) derivative of the spin-weighted spherical
harmonic coefficients `flm` with spin weight `s`. `L = lmax+1` is the
number of modes in `flm`. The result has spin weight `s+1`.

The arrays `ðflm` and `flm` must have length `L^2`.

See
[Wikipedia](https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics)
for a definition of the `ð` (eth) operator.

See also: [`SSHT.eth`](@ref), [`SSHT.ethbar!`](@ref),
[`SSHT.core_dh_transform_sov!`](@ref),
[`SSHT.core_dh_inverse_sov!`](@ref).
"""
function eth!(ðflm::AbstractVector, flm::AbstractVector, L::Int, spin::Int)
    0 < L || throw(DomainError(L, "Need 0 < L"))
    T = typeof(abs(zero(eltype(flm))))
    # @assert T <: Real
    lmin = max(abs(spin), abs(spin + 1))
    lmax = L - 1
    for l in 0:min(lmin - 1, lmax), m in (-l):l
        # @assert (l - spin) * (l + spin + 1) ≤ 0
        ind = sampling_elm2ind(l, m)
        ðflm[ind] = 0
    end
    for l in lmin:lmax, m in (-l):l
        # @assert (l - spin) * (l + spin + 1) > 0
        ind = sampling_elm2ind(l, m)
        # <https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics>
        # ð sYlm = + sqrt((l-s)*(l+s+1)) (s+1)Ylm
        ðflm[ind] = sqrt(T((l - spin) * (l + spin + 1))) * flm[ind]
    end
    return ðflm
end
eth!(ðflm::AbstractVector, flm::AbstractVector, L::Integer, spin::Integer) = eth!(ðflm, flm, Int(L), Int(spin))

"""
    ðflm = eth(flm::AbstractVector, L::Integer, spin::Integer)
    ðflm::AbstractVector

Calculate the `ð` (eth) derivative of the spin-weighted spherical
harmonic coefficients `flm` with spin weight `s`. `L = lmax+1` is the
number of modes in `flm`. The result has spin weight `s+1`.

The array `flm` must have length `L^2`. The result `ðflm` will also
have length `L^2`.

See
[Wikipedia](https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics)
for a definition of the `ð` (eth) operator.

See also: [`SSHT.eth`](@ref), [`SSHT.ethbar!`](@ref),
[`SSHT.core_dh_transform_sov!`](@ref),
[`SSHT.core_dh_inverse_sov!`](@ref).
"""
eth(flm::AbstractVector, L::Integer, spin::Integer) = eth!(similar(flm), flm, L, spin)

"""
    ethbar!(ð̄flm::AbstractVector, flm::AbstractVector, L::Integer, spin::Integer)

Calculate the `ð̄` (eth-bar) derivative of the spin-weighted spherical
harmonic coefficients `flm` with spin weight `s`. `L = lmax+1` is the
number of modes in `flm`. The result has spin weight `s-1`.

The arrays `ð̄flm` and `flm` must have length `L^2`.

See
[Wikipedia](https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics)
for a definition of the `ð̄` (eth-bar) operator.

See also: [`SSHT.ethbar`](@ref), [`SSHT.eth!`](@ref),
[`SSHT.core_dh_transform_sov!`](@ref),
[`SSHT.core_dh_inverse_sov!`](@ref).
"""
function ethbar!(ð̄flm::AbstractVector, flm::AbstractVector, L::Int, spin::Int)
    0 < L || throw(DomainError(L, "Need 0 < L"))
    T = typeof(abs(zero(eltype(flm))))
    # @assert T <: Real
    lmin = max(abs(spin), abs(spin - 1))
    lmax = L - 1
    for l in 0:min(lmin - 1, lmax), m in (-l):l
        # @assert (l + spin) * (l - spin + 1) ≤ 0
        ind = sampling_elm2ind(l, m)
        ð̄flm[ind] = 0
    end
    for l in lmin:lmax, m in (-l):l
        # @assert (l + spin) * (l - spin + 1) > 0
        ind = sampling_elm2ind(l, m)
        # <https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics>
        # ð̄ sYlm = - sqrt((l+s)*(l-s+1)) (s-1)Ylm
        ð̄flm[ind] = -sqrt(T((l + spin) * (l - spin + 1))) * flm[ind]
    end
    return ð̄flm
end
ethbar!(ð̄flm::AbstractVector, flm::AbstractVector, L::Integer, spin::Integer) = ethbar!(ð̄flm, flm, Int(L), Int(spin))

"""
    ð̄flm = ethbar(flm::AbstractVector, L::Integer, spin::Integer)
    ð̄flm::AbstractVector

Calculate the `ð̄` (eth) derivative of the spin-weighted spherical
harmonic coefficients `flm` with spin weight `s`. `L = lmax+1` is the
number of modes in `flm`. The result has spin weight `s-1`.

The array `flm` must have length `L^2`. The result `ð̄flm` will also
have length `L^2`.

See
[Wikipedia](https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics)
for a definition of the `ð` (eth) operator.

See also: [`SSHT.ethbar`](@ref), [`SSHT.eth!`](@ref),
[`SSHT.core_dh_transform_sov!`](@ref),
[`SSHT.core_dh_inverse_sov!`](@ref).
"""
ethbar(flm::AbstractVector, L::Integer, spin::Integer) = ethbar!(similar(flm), flm, L, spin)

################################################################################

export ash_grid_size, ash_nmodes
function ash_grid_size(lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    L = Int(lmax) + 1
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    return (nphi, ntheta)::NTuple{2,Int}
end
function ash_nmodes(lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    L = Int(lmax) + 1
    return (L^2,)::NTuple{1,Int}
end

export ash_ntheta, ash_nphi, ash_thetas, ash_phis, ash_point_coord, ash_point_delta, ash_grid_as_phi_theta
ash_ntheta(lmax) = ash_grid_size(lmax)[2]::Int
ash_nphi(lmax) = ash_gri_size(lmax)[1]::Int
function ash_thetas(lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    L = Int(lmax) + 1
    ntheta = sampling_dh_ntheta(L)
    return [sampling_dh_t2theta(t, L) for t in 1:ntheta]
end
function ash_phis(lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    L = Int(lmax) + 1
    nphi = sampling_dh_nphi(L)
    return [sampling_dh_p2phi(p, L) for p in 1:nphi]
end
function ash_point_coord(ij::Union{CartesianIndex{2},NTuple{2,Int}}, lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    L = Int(lmax) + 1
    p, t = Tuple(ij)
    return sampling_dh_t2theta(t, L), sampling_dh_p2phi(p, L)
end
function ash_point_delta(ij::Union{CartesianIndex{2},NTuple{2,Int}}, lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    L = Int(lmax) + 1
    p, t = Tuple(ij)
    nphi = sampling_dh_nphi(L)
    theta = sampling_dh_t2theta(t, L)
    dtheta = sampling_weight_dh(theta, L) / sin(theta)
    dphi = 2π / nphi
    return dtheta, dphi
end
ash_grid_as_phi_theta(grid::AbstractMatrix) = grid

export ash_mode_index
function ash_mode_index(s::Integer, l::Integer, m::Integer, lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    abs(s) ≤ l ≤ lmax || throw(DomainError(l, "Need abs(s) ≤ l ≤ lmax"))
    -l ≤ m ≤ l || throw(DomainError(m, "Need -l ≤ m ≤ l"))
    return CartesianIndex(sampling_elm2ind(l, m))::CartesianIndex{1}
end
export ash_mode_numbers
function ash_mode_numbers(s::Int, ind::NTuple{1,<:Int}, lmax::Int)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    l, m = sampling_ind2elm(ind[1])
    return (l, m)::NTuple{2,Int}
end
ash_mode_numbers(s::Integer, ind::NTuple{1,<:Integer}, lmax::Integer) = ash_mode_numbers(Int(s), NTuple{1,Int}(ind), Int(lmax))

export ash_transform!, ash_transform, ash_evaluate!, ash_evaluate
function ash_transform!(flm::AbstractArray{<:Complex}, f::AbstractMatrix{<:Complex}, s::Integer, lmax::Integer)
    L = lmax + 1
    return core_dh_forward_sov!(flm, f, L, s)
end
function ash_transform(f::AbstractMatrix{<:Complex}, s::Integer, lmax::Integer)
    L = lmax + 1
    return core_dh_forward_sov(f, L, s)
end

function ash_evaluate!(f::AbstractMatrix{<:Complex}, flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    L = lmax + 1
    return core_dh_inverse_sov!(f, flm, L, s)
end
function ash_evaluate(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    L = lmax + 1
    return core_dh_inverse_sov(flm, L, s)
end

export ash_eth!, ash_eth, ash_ethbar!, ash_ethbar
function ash_eth!(ðflm::AbstractArray{<:Complex}, flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    L = lmax + 1
    return eth!(ðflm, flm, L, s)
end
function ash_eth(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    L = lmax + 1
    return eth(flm, L, s)
end

function ash_ethbar!(ðflm::AbstractArray{<:Complex}, flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    L = lmax + 1
    return ethbar!(ðflm, flm, L, s)
end
function ash_ethbar(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    L = lmax + 1
    return ethbar(flm, L, s)
end

end
