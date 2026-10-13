# CLAUDE.md

Guidance for coding agents working on SSHT.jl.

## What this package is

SSHT.jl is a thin Julia wrapper around the C library
[ssht](https://astro-informatics.github.io/ssht/) (spin-weighted
spherical harmonic transforms). The C library comes from the binary
package `ssht_jll`; it is called directly via `ccall((:ssht_..., libssht), ...)`.
There are no other runtime dependencies.

## Layout

- `src/SSHT.jl` — the whole package, in a single module. Sections follow
  the C headers (`# File ssht_core.h`, `# File ssht_sampling.h`, ...),
  followed by pure-Julia code (`eth`, `ethbar`).
- `test/runtests.jl` — the test suite; test-only deps (`Random`, `Test`)
  are in `test/Project.toml`, not in the main `Project.toml`.
- `docs/` — Documenter.jl setup; `docs/src/index.md` lists every
  documented function in `@docs` blocks.
- `figures/` — images used by `README.md`.

## Commands

```sh
# Run the tests (takes about 10 seconds)
julia --project=. -e 'using Pkg; Pkg.test()'

# Build the documentation locally (output in docs/build/).
# Pkg.develop adds SSHT (with a local path) to docs/Project.toml; revert that afterwards.
julia --project=docs -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=docs docs/make.jl
git checkout docs/Project.toml

# Format the code (settings in .JuliaFormatter.toml: YAS style, 4-space indent, 132-char margin)
julia -e 'using JuliaFormatter; format(".")'
```

## Conventions

- **Function pairs.** Most functions come as a mutating `foo!(out, in, ...)`,
  which checks array sizes, calls the C function, and returns `out`, plus
  an allocating `foo(in, ...)` that allocates the output and calls `foo!`.
  Keep this pattern for new wrappers.
- **Names** mirror the C names without the `ssht_` prefix
  (`ssht_core_dh_forward_sov` → `SSHT.core_dh_forward_sov`). Nothing is
  exported; users call `SSHT.foo`.
- **Indexing.** The Julia API is 1-based; the C API is 0-based. Wrappers
  convert (e.g. `sampling_dh_p2phi(p, L)` passes `p - 1`). Mode indices
  use `ind = l^2 + l + m + 1` (`sampling_elm2ind` / `sampling_ind2elm`,
  implemented in Julia, not via `ccall`).
- **Array layout.** Grid arrays have size `(nphi, ntheta)` (phi varies
  fastest, matching ssht's row-major `[theta][phi]` layout); mode arrays
  have length `L^2`, where `L = lmax + 1`. Arrays are passed to C as raw
  pointers, so they must be dense, column-major `Float64` /
  `Complex{Float64}` arrays.
- **Spin.** The complex transforms take a `spin` argument; the
  real-valued `*_real` transforms do not (they are spin 0).
- **Integer types.** C `int` arguments are declared `Cint` in the `ccall`
  signature. Public methods accept any `Integer` and convert to `Int`
  internally. CI also runs on 32-bit x86, where `Int == Int32`, so do
  not assume `Int == Int64`.
- **Docstrings.** Every public function has a docstring with a
  hand-aligned signature block, size requirements, and a "See also" list
  using `[`SSHT.foo`](@ref)` links. When you add a docstring, also add the
  function to the appropriate `@docs` block in `docs/src/index.md`;
  Documenter fails the build on documented functions that are not
  included in the manual. Update `README.md` if user-visible usage changes.
- **Julia compatibility.** `[compat]` allows Julia 1.6, and CI tests 1.6
  through the current release on Linux, macOS, and Windows (x64, arm64,
  x86). Do not use language or Base features newer than Julia 1.6
  (e.g. `public`, package extensions, `stack`, `Returns`) without
  raising the lower bound in `Project.toml` and the CI matrix in
  `.github/workflows/CI.yml`.

## Testing

- Every wrapped function should have tests. The test suite checks
  transforms against hard-coded closed forms of the spin-weighted
  spherical harmonics `sYlm` (and their `ð` / `ð̄` derivatives) for
  `l ≤ 2`, `|s| ≤ 2`, plus randomized checks (linearity, orthonormality)
  with fixed `Random.seed!` values. Follow the same style for new tests.
- Use `isapprox` with explicit `atol` (multiples of `eps()`) when
  comparing results that should be zero.

## Commits and releases

- Run the tests (and the formatter) before committing.
- Do not add `Co-Authored-By` trailers for AI agents to commit messages or PR
  descriptions.
- Releases: bump `version` in `Project.toml`, register via the Julia
  General registry (JuliaRegistrator); TagBot then creates the git tag
  and GitHub release. Do not create release tags by hand.
