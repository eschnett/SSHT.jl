module ssht

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
function core_dh_inverse_sov_real(Lflm::AbstractVector{Complex{Float64}}, spin::Integer, verbosity::Integer=0)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    f = Array{Float64}(undef, nphi, ntheta)
    core_dh_inverse_sov_real!(f, flm, L, spin, verbosity)
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
    0 ≤ el || throw(DomainError())
    -el ≤ m ≤ el || throw(DomainError())
    return el^2 + el + m + 1
end
sampling_elm2ind(el::Integer, m::Integer) = sampling_elm2ind(Int(el), Int(m))
function sampling_ind2elm(ind::Int)
    1 ≤ ind || throw(DomainError())
    el = isqrt(ind - 1)
    m = ind - el^2 - el - 1
    return el, m
end
sampling_ind2elm(ind::Integer) = sampling_ind2elm(Int(ind))

################################################################################

function eth!(ðflm::AbstractVector, flm::AbstractVector, L::Int, spin::Int)
    0 < L || throw(DomainError())
    for l in 0:min(abs(spin), L - 1), m in (-l):l
        ind = sampling_elm2ind(l, m)
        ðflm[ind] = 0
    end
    for l in (abs(spin) + 1):(L - 1), m in (-l):l
        ind = sampling_elm2ind(l, m)
        # <https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics>
        # ð sYlm = + sqrt((l-s)*(l+s+1)) (s+1)Ylm
        ðflm[ind] = sqrt((l - spin) * (l + spin + 1)) * flm[ind]
    end
    return ðflm
end
eth!(ðflm::AbstractVector, flm::AbstractVector, L::Integer, spin::Integer) = eth!(ðflm, flm, Int(L), Int(spin))
eth(flm::AbstractVector, L::Integer, spin::Integer) = eth!(similar(flm), flm, L, spin)

function ethbar!(ð̄flm::AbstractVector, flm::AbstractVector, L::Int, spin::Int)
    0 < L || throw(DomainError())
    for l in 0:min(abs(spin), L - 1), m in (-l):l
        ind = sampling_elm2ind(l, m)
        ð̄flm[ind] = 0
    end
    for l in (abs(spin) + 1):(L - 1), m in (-l):l
        ind = sampling_elm2ind(l, m)
        # <https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics>
        # ð̄ sYlm = - sqrt((l+s)*(l-s+1)) (s-1)Ylm
        ð̄flm[ind] = -sqrt((l + spin) * (l - spin + 1)) * flm[ind]
    end
    return ð̄flm
end
ethbar!(ð̄flm::AbstractVector, flm::AbstractVector, L::Integer, spin::Integer) = ethbar!(ð̄flm, flm, Int(L), Int(spin))
ethbar(flm::AbstractVector, L::Integer, spin::Integer) = ethbar!(similar(flm), flm, L, spin)

end
