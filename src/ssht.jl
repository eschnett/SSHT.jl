module ssht

using ssht_jll

# File ssht_adjoint.h

# File ssht_core.h

function core_dh_inverse_sov!(f::Array{Complex{Float64},2}, flm::Vector{Complex{Float64}}, L::Integer, spin::Integer,
                              verbosity::Integer)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    any(size(f) .< (nphi, ntheta)) && throw(DimensionMismatch())
    length(flm) < L^2 && throw(DimensionMismatch())
    ccall((:ssht_core_dh_inverse_sov, libssht), Cvoid, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Cint, Cint, Cint), f, flm, L,
          spin, verbosity)
    return nothing
end

function core_dh_inverse_sov_real!(f::Array{Float64,2}, flm::Vector{Complex{Float64}}, L::Integer, verbosity::Integer)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    any(size(f) .< (nphi, ntheta)) && throw(DimensionMismatch())
    length(flm) < L^2 && throw(DimensionMismatch())
    ccall((:ssht_core_dh_inverse_sov_real, libssht), Cvoid, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Cint, Cint), f, flm, L,
          verbosity)
    return nothing
end

function core_dh_forward_sov!(flm::Vector{Complex{Float64}}, f::Array{Complex{Float64},2}, L::Integer, spin::Integer,
                              verbosity::Integer)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    length(flm) < L^2 && throw(DimensionMismatch())
    any(size(f) .< (nphi, ntheta)) && throw(DimensionMismatch())
    ccall((:ssht_core_dh_forward_sov, libssht), Cvoid, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Cint, Cint, Cint), flm, f, L,
          spin, verbosity)
    return nothing
end

function core_dh_forward_sov_real!(flm::Vector{Complex{Float64}}, f::Array{Float64,2}, L::Integer, verbosity::Integer)
    nphi = sampling_dh_nphi(L)
    ntheta = sampling_dh_ntheta(L)
    length(flm) < L^2 && throw(DimensionMismatch())
    any(size(f) .< (nphi, ntheta)) && throw(DimensionMismatch())
    ccall((:ssht_core_dh_forward_sov_real, libssht), Cvoid, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Cint, Cint), flm, f, L,
          verbosity)
    return nothing
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

"L is number of modes; 0 ≤ l < L"
sampling_dh_n(L::Integer) = Int(ccall((:ssht_sampling_dh_n, libssht), Cint, (Cint,), L))
sampling_dh_nphi(L::Integer) = Int(ccall((:ssht_sampling_dh_nphi, libssht), Cint, (Cint,), L))
sampling_dh_ntheta(L::Integer) = Int(ccall((:ssht_sampling_dh_ntheta, libssht), Cint, (Cint,), L))
sampling_dh_p2phi(p::Integer, L::Integer) = Float64(ccall((:ssht_sampling_dh_p2phi, libssht), Cdouble, (Cint, Cint), p - 1, L))
sampling_dh_t2theta(t::Integer, L::Integer) = Float64(ccall((:ssht_sampling_dh_t2theta, libssht), Cdouble, (Cint, Cint), t - 1, L))
function sampling_weight_dh(theta_t::Real, L::Integer)
    return Float64(ccall((:ssht_sampling_weight_dh, libssht), Cdouble, (Cdouble, Cint), theta_t, L))
end

# There are L^2 coefficients
function sampling_elm2ind(el::Int, m::Int)
    0 ≤ el || throw(DomainError())
    -el ≤ m ≤ el || throw(DomainError())
    return el * el + el + m + 1
end
sampling_elm2ind(el::Integer, m::Integer) = sampling_elm2ind(Int(el), Int(m))
function sampling_ind2elm(ind::Int)
    1 ≤ ind || throw(DomainError())
    el = isqrt(ind - 1)
    m = ind - 1 - el^2 - el
    return (el, m)
end
sampling_ind2elm(ind::Integer) = sampling_ind2elm(Int(ind))

end
