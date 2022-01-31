# SSHT.jl: Fast and exact spin-weighted spherical harmonic transforms

This package is a Julia wrapper for the
[SSHT](https://astro-informatics.github.io/ssht/) library. It
calculates spin-weighted wpherical harmonic transforms.

* [![Documenter](https://img.shields.io/badge/docs-dev-blue.svg)](https://eschnett.github.io/SSHT.jl/dev)
* [![GitHub
  CI](https://github.com/eschnett/SSHT.jl/workflows/CI/badge.svg)](https://github.com/eschnett/SSHT.jl/actions)
* [![Codecov](https://codecov.io/gh/eschnett/SSHT.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/eschnett/SSHT.jl)

## Example

Plot the `s=+1, l=1, m=+1` spin-weighted spherical harmonic function:
```Julia
using SSHT

lmax = 40

L = lmax + 1
nphi = SSHT.sampling_dh_nphi(L)
ntheta = SSHT.sampling_dh_ntheta(L)

phis = [SSHT.sampling_dh_p2phi(p, L) for p in 1:nphi];
thetas = [SSHT.sampling_dh_t2theta(t, L) for t in 1:ntheta];

flm = zeros(Complex{Float64}, L^2);
flm[SSHT.sampling_elm2ind(1, +1)] = 1;

f = SSHT.core_dh_inverse_sov(flm, L, +1);

################################################################################

using GLMakie

fig = Figure(; resolution=(1000, 300));

Axis(fig[1, 1]; title="real(f)")
Axis(fig[1, 3]; title="imag(f)")
hm = heatmap!(fig[1, 1], phis, thetas, real.(f); colormap=:magma)
Colorbar(fig[1, 2], hm)
hm = heatmap!(fig[1, 3], phis, thetas, imag.(f); colormap=:magma)
Colorbar(fig[1, 4], hm)
rowsize!(fig.layout, 1, Aspect(1, 1 / 2))

display(fig)
```

![s=_1, l=1, m=+1 mode](https://github.com/eschnett/SSHT.jl/blob/main/figures/sYlm.png)


## Related packages

- [FFTW.jl](https://github.com/JuliaMath/FFTW.jl)
- [FastSphericalHarmonics.jl](https://github.com/eschnett/FastSphericalHarmonics.jl)
- [FastTransforms.jl](https://github.com/JuliaApproximation/FastTransforms.jl)
- [SHTOOLS.jl](https://github.com/eschnett/SHTOOLS.jl)
- [SphericalFunctions.jl](https://github.com/moble/SphericalFunctions.jl)
