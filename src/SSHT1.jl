module SSHT

using ssht_jll

# File ssht_adjoint.h

# File ssht_core.h

# L is number of modes; 0 ≤ l < L
# verbosity: 0 ≤ v ≤ 5

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
function core_dh_inverse_sov(flm::AbstractVector{Complex{Float64}}, L::Integer, spin::Integer, verbosity::Integer=0)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    f = Array{Complex{Float64}}(undef, nphi, ntheta)
    core_dh_inverse_sov!(f, flm, L, spin, verbosity)
    return f
end

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
function core_dh_inverse_sov_real(flm::AbstractVector{Complex{Float64}}, L::Integer, verbosity::Integer=0)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    f = Array{Float64}(undef, nphi, ntheta)
    core_dh_inverse_sov_real!(f, flm, L, verbosity)
    return f
end

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
function core_dh_forward_sov(f::AbstractArray{Complex{Float64},2}, L::Integer, spin::Integer, verbosity::Integer=0)
    flm = Array{Complex{Float64}}(undef, L^2)
    core_dh_forward_sov!(flm, f, L, spin, verbosity)
    return flm
end

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

sampling_dh_n(L::Integer) = Int(ccall((:ssht_sampling_dh_n, libssht), Cint, (Cint,), L))
sampling_dh_nphi(L::Integer) = Int(ccall((:ssht_sampling_dh_nphi, libssht), Cint, (Cint,), L))
sampling_dh_ntheta(L::Integer) = Int(ccall((:ssht_sampling_dh_ntheta, libssht), Cint, (Cint,), L))
sampling_dh_p2phi(p::Integer, L::Integer) = Float64(ccall((:ssht_sampling_dh_p2phi, libssht), Cdouble, (Cint, Cint), p - 1, L))
sampling_dh_t2theta(t::Integer, L::Integer) = Float64(ccall((:ssht_sampling_dh_t2theta, libssht), Cdouble, (Cint, Cint), t - 1, L))
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

# There are L^2 coefficients
function sampling_elm2ind(el::Int, m::Int)
    0 ≤ el || throw(DomainError(el, "Need 0 ≤ el"))
    -el ≤ m ≤ el || throw(DomainError(m, "Need -el ≤ m ≤ el"))
    return el^2 + el + m + 1
end
sampling_elm2ind(el::Integer, m::Integer) = sampling_elm2ind(Int(el), Int(m))
function sampling_ind2elm(ind::Int)
    1 ≤ ind || throw(DomainError(ind, "Need 1 ≤ ind"))
    el = isqrt(ind - 1)
    m = ind - el^2 - el - 1
    return el, m
end
sampling_ind2elm(ind::Integer) = sampling_ind2elm(Int(ind))

################################################################################

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
eth(flm::AbstractVector, L::Integer, spin::Integer) = eth!(similar(flm), flm, L, spin)

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
ethbar(flm::AbstractVector, L::Integer, spin::Integer) = ethbar!(similar(flm), flm, L, spin)

################################################################################

export ash_grid_size, ash_nmodes
function ash_grid_size(lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    L = Int(lmax) + 1
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    return nphi, ntheta
end
function ash_nmodes(lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    L = Int(lmax) + 1
    return (L^2,)
end

export ash_ntheta, ash_nphi, ash_thetas, ash_phis, ash_point_coord, ash_point_delta, ash_grid_as_phi_theta
ash_ntheta(lmax) = ash_grid_size(lmax)[2]
ash_nphi(lmax) = ash_gri_size(lmax)[1]
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
function ash_point_coord(ij::CartesianIndex{2}, lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    L = Int(lmax) + 1
    p, t = Tuple(ij)
    return sampling_dh_t2theta(t, L), sampling_dh_p2phi(p, L)
end
function ash_point_delta(ij::CartesianIndex{2}, lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    L = Int(lmax) + 1
    p, t = Tuple(ij)
    nphi = ssht.sampling_dh_nphi(L)
    theta = ssht.sampling_dh_t2theta(t, L)
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
    return sampling_elm2ind(l, m)
end

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
