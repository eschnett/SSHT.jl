# SSHT.jl: Fast and exact spin spherical harmonic transforms

This package is a Julia wrapper for the
[SSHT](https://astro-informatics.github.io/ssht/) library. It
calculates spin spherical harmonic transforms.

Most functions come in two versions, one that mutates its arguments
and one that allocates its output.

## Spin-weighted Spherical Harmonic transforms

The functions with `dh` in their names use the Driscoll & Healy
quadrature points. With `L = lmax+1` modes there are `nphi * ntheta`
quadrature points on the sphere, with `nphi = 2L-1` and `ntheta = 2L`.
The points are equispaced in the angles `theta` (latitude) and `phi`
(longitude), and they straddle (avoid) the poles.

```@docs
SSHT.core_dh_inverse_sov!
SSHT.core_dh_inverse_sov
SSHT.core_dh_inverse_sov_real!
SSHT.core_dh_inverse_sov_real
SSHT.core_dh_forward_sov!
SSHT.core_dh_forward_sov
SSHT.core_dh_forward_sov_real!
SSHT.core_dh_forward_sov_real
```

## Derivatives ð (eth) and ð̄ (eth-bar)

See
[Wikipedia](https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics)
for the definition of these operators.

```@docs
SSHT.eth!
SSHT.eth
SSHT.ethbar!
SSHT.ethbar
```

## Helper functions

```@docs
SSHT.sampling_dh_n
SSHT.sampling_dh_nphi
SSHT.sampling_dh_ntheta
SSHT.sampling_dh_p2phi
SSHT.sampling_dh_t2theta
SSHT.sampling_weight_dh
SSHT.sampling_elm2ind
SSHT.sampling_ind2elm
```
