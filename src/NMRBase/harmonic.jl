# Parametric signal analysis ######################################################################
#
# Two distinct methods on 1-D complex signals, kept separate because they are used for different
# things (and the earlier code wrongly used the first for both):
#
#   Linear prediction (`linearpredict`) — the standard autoregressive method (as in nmrPipe /
#   nmrglue). Solve for LP coefficients `aⱼ` such that `x[n] ≈ Σⱼ aⱼ x[n∓j]`, optionally reflect the
#   roots of the prediction polynomial for stability, then EXTRAPOLATE by recurrence. Cheap per
#   vector (small order), no pole/amplitude fitting.
#
#   Harmonic inversion / HSVD (`solventfilter(method=:hsvd)`) — the state-space (ESPRIT) method:
#   model `x[n] ≈ Σₖ cₖ zₖⁿ` as damped exponentials via SVD of a Hankel matrix, giving explicit
#   frequencies used to subtract only the components in a chosen band.
#
# Kernels operate on plain complex vectors (testable in isolation); multicomplex data is fed one
# imₖ-plane vector at a time by `_maptransform`.

# ---------------------------------------------------------------------------------------------------
# Linear prediction (autoregressive coefficient method)
# ---------------------------------------------------------------------------------------------------

"""
    _lpcoefficients(x, order, direction) -> Vector{Complex}

Least-squares LP coefficients `a` of the given `order` for complex vector `x`. For `:forward`,
`x[n] ≈ Σⱼ aⱼ x[n-j]`; for `:backward`, `x[n] ≈ Σⱼ aⱼ x[n+j]`. Solved by QR (`D \\ d`).
"""
function _lpcoefficients(x::AbstractVector{<:Complex}, order::Integer, direction::Symbol)
    N = length(x)
    m = clamp(order, 1, N - 1)
    K = N - m
    D = Matrix{eltype(x)}(undef, K, m)
    d = Vector{eltype(x)}(undef, K)
    if direction === :forward
        @inbounds for i in 1:K
            for j in 1:m
                D[i, j] = x[m + i - j]     # x[n-1 … n-m], n = m+i
            end
            d[i] = x[m + i]
        end
    else                                    # backward
        @inbounds for i in 1:K
            for j in 1:m
                D[i, j] = x[i + j]          # x[n+1 … n+m], n = i
            end
            d[i] = x[i]
        end
    end
    return D \ d
end

# polynomial helpers for root reflection --------------------------------------------------------

"coefficients (descending, leading first) of ∏(z − rₖ)"
function _polyfromroots(r::AbstractVector{<:Complex})
    p = Complex{Float64}[1]
    for rk in r
        q = similar(p, length(p) + 1)
        fill!(q, 0)
        for i in eachindex(p)
            q[i] += p[i]
            q[i + 1] -= p[i] * rk
        end
        p = q
    end
    return p
end

"roots of the monic-normalised polynomial with descending coefficients `pc`"
function _polyroots(pc::AbstractVector{<:Complex})
    a = pc ./ pc[1]
    m = length(a) - 1
    m == 0 && return Complex{Float64}[]
    Cm = zeros(Complex{Float64}, m, m)
    Cm[1, :] .= -a[2:end]
    for i in 2:m
        Cm[i, i - 1] = 1
    end
    return eigvals(Cm)
end

"""
    _reflectroots(a, direction) -> Vector{Complex}

Stabilise LP coefficients `a` by reflecting the roots of the prediction polynomial across the unit
circle: for `:forward` extension, roots with `|z|>1` (growing) are mapped to `1/conj(z)` (inside);
for `:backward`, roots with `|z|<1` are pushed outside. Coefficients are rebuilt from the reflected
roots.
"""
function _reflectroots(a::AbstractVector{<:Complex}, direction::Symbol)
    # recurrence characteristic polynomial: zᵐ − a₁zᵐ⁻¹ − … − aₘ  →  coeffs [1, −a₁, …, −aₘ]
    r = _polyroots(vcat(one(eltype(a)), -a))
    for k in eachindex(r)
        ρ = abs(r[k])
        if direction === :forward && ρ > 1
            r[k] = 1 / conj(r[k])
        elseif direction === :backward && 0 < ρ < 1
            r[k] = 1 / conj(r[k])
        end
    end
    pc = _polyfromroots(r)          # [1, −a₁′, …]
    return -pc[2:end]
end

"""
    _lpextrapolate(x, a, npred, direction) -> Vector{Complex}

Extend `x` by `npred` points by iterating the LP recurrence with coefficients `a`: `:forward` appends
`Σⱼ aⱼ x[end-j+1]`, `:backward` prepends `Σⱼ aⱼ x[j]`.
"""
function _lpextrapolate(x::AbstractVector{<:Complex}, a::AbstractVector{<:Complex}, npred::Integer,
                        direction::Symbol)
    m = length(a)
    y = collect(x)
    if direction === :forward
        for _ in 1:npred
            v = zero(eltype(y))
            @inbounds for j in 1:m
                v += a[j] * y[end - j + 1]
            end
            push!(y, v)
        end
    else
        for _ in 1:npred
            v = zero(eltype(y))
            @inbounds for j in 1:m
                v += a[j] * y[j]
            end
            pushfirst!(y, v)
        end
    end
    return y
end

"""
    _zhubaxforward(x, order, rootfix) -> Vector{Complex}

Forward-extension LP coefficients by the Zhu–Bax scheme: fit *backward* LP, treat the roots of its
prediction polynomial (`Σⱼ bⱼwʲ = 1`) as the candidate poles, reflect any with `|z|>1` to `1/conj(z)`
so the forward extrapolation cannot grow, and rebuild forward coefficients. This separates signal
poles (inside the unit circle) from noise poles more cleanly than a direct forward fit, and is the
variant NMR usually prefers. With `rootfix=false`, falls back to a plain forward least-squares fit.
"""
function _zhubaxforward(x::AbstractVector{<:Complex}, order::Integer, rootfix::Bool)
    rootfix || return _lpcoefficients(x, order, :forward)
    b = _lpcoefficients(x, order, :backward)
    poles = _polyroots(vcat(reverse(b), -one(eltype(b))))     # roots of Σ bⱼ wʲ − 1
    poles = [abs(z) > 1 ? 1 / conj(z) : z for z in poles]     # keep inside for stable forward
    return -_polyfromroots(poles)[2:end]                      # rebuild forward coefficients
end

"""
    _lppredict(x; npred, order, direction, rootfix) -> Vector{Complex}

Linear-prediction extension of `x` by `npred` points, extrapolating by recurrence. Forward extension
uses the Zhu–Bax backward-fit-and-reflect coefficients (stable); backward extension uses a direct
backward fit with root reflection.
"""
function _lppredict(x::AbstractVector{<:Complex}; npred::Integer, order::Integer,
                    direction::Symbol=:forward, rootfix::Bool=true)
    if direction === :forward
        a = _zhubaxforward(x, order, rootfix)
    else
        a = _lpcoefficients(x, order, :backward)
        rootfix && (a = _reflectroots(a, :backward))
    end
    return _lpextrapolate(x, a, npred, direction)
end

"""
    _lprepair(x; nrepair, order, rootfix) -> Vector{Complex}

Replace the first `nrepair` points of `x` by backward prediction from the remaining points (length
unchanged).
"""
function _lprepair(x::AbstractVector{<:Complex}; nrepair::Integer, order::Integer,
                   rootfix::Bool=true)
    good = collect(@view x[(nrepair + 1):end])
    a = _lpcoefficients(good, order, :backward)
    rootfix && (a = _reflectroots(a, :backward))
    return _lpextrapolate(good, a, nrepair, :backward)   # prepends nrepair points
end

# Pooled (shared-coefficient) LP ------------------------------------------------------------------
#
# Vectors along a dimension share the same poles ⇒ the same LP recurrence, so the coefficients are
# estimated ONCE from a stacked fit of the highest-energy vectors and the (cheap, LAPACK-free)
# recurrence is applied to every column. One QR solve + one eigvals instead of one per column.

"indices of the `nselect` highest-energy columns of `mat` (or all, if fewer)"
function _selectcols(mat::AbstractMatrix, nselect::Integer)
    size(mat, 2) <= nselect && return collect(axes(mat, 2))
    return partialsortperm(vec(sum(abs2, mat; dims=1)), 1:nselect; rev=true)
end

"stacked least-squares LP coefficients (order `m`, direction `fitdir`) from columns `sel` of `mat`"
function _pooledfit(mat::AbstractMatrix{<:Complex}, sel, m::Integer, fitdir::Symbol)
    N = size(mat, 1)
    K = N - m
    D = Matrix{ComplexF64}(undef, length(sel) * K, m)
    d = Vector{ComplexF64}(undef, length(sel) * K)
    row = 0
    @inbounds for jj in sel
        x = view(mat, :, jj)
        if fitdir === :backward
            for i in 1:K, j in 1:m
                D[row + i, j] = x[i + j]
            end
            for i in 1:K
                d[row + i] = x[i]
            end
        else
            for i in 1:K, j in 1:m
                D[row + i, j] = x[m + i - j]
            end
            for i in 1:K
                d[row + i] = x[m + i]
            end
        end
        row += K
    end
    return D \ d
end

"shared extension coefficients for `direction` (Zhu–Bax for forward), pooled over columns"
function _pooledcoeffs(mat, m, direction, rootfix, nselect)
    sel = _selectcols(mat, nselect)
    if direction === :forward
        rootfix || return _pooledfit(mat, sel, m, :forward)
        b = _pooledfit(mat, sel, m, :backward)
        poles = _polyroots(vcat(reverse(b), -one(eltype(b))))
        poles = [abs(z) > 1 ? 1 / conj(z) : z for z in poles]
        return -_polyfromroots(poles)[2:end]
    else
        a = _pooledfit(mat, sel, m, :backward)
        rootfix && (a = _reflectroots(a, :backward))
        return a
    end
end

"""
    _lpmatrix(mat; npred, order, direction, rootfix, pool, nselect) -> Matrix

Extend every column of `mat` by `npred` points. With `pool` (default), one shared set of
coefficients is estimated from the `nselect` highest-energy columns and applied to all columns by
recurrence; otherwise each column is fitted independently (matching nmrglue).
"""
function _lpmatrix(mat::AbstractMatrix{<:Complex}; npred::Integer, order::Integer,
                   direction::Symbol=:forward, rootfix::Bool=true, pool::Bool=true,
                   nselect::Integer=128)
    N, ncols = size(mat)
    m = clamp(order, 1, N - 1)
    pool ||
        return _percolumn(x -> _lppredict(x; npred=npred, order=m, direction=direction,
                                          rootfix=rootfix))(mat)

    a = _pooledcoeffs(mat, m, direction, rootfix, nselect)
    out = Matrix{ComplexF64}(undef, N + npred, ncols)
    @inbounds for jc in 1:ncols
        if direction === :forward
            for i in 1:N
                out[i, jc] = mat[i, jc]
            end
            for k in 1:npred
                v = zero(ComplexF64)
                for j in 1:m
                    v += a[j] * out[N + k - j, jc]
                end
                out[N + k, jc] = v
            end
        else
            for i in 1:N
                out[npred + i, jc] = mat[i, jc]
            end
            for k in npred:-1:1
                v = zero(ComplexF64)
                for j in 1:m
                    v += a[j] * out[k + j, jc]
                end
                out[k, jc] = v
            end
        end
    end
    return out
end

# ---------------------------------------------------------------------------------------------------
# Harmonic inversion / HSVD (used by solvent suppression)
# ---------------------------------------------------------------------------------------------------

# cap the Hankel row count so the per-vector SVD stays cheap on long FIDs.
const _HANKEL_LMAX = 256

"""
    _hankelpoles(x, K; L=0) -> Vector{Complex}

Estimate up to `K` signal poles `zₖ` of the complex vector `x` by ESPRIT: SVD of the Hankel matrix,
then the eigenvalues of the shift-invariance operator on the truncated left singular subspace. `L`
(Hankel rows) is capped at `_HANKEL_LMAX`; the model order `K` sets the eigenproblem size.
"""
function _hankelpoles(x::AbstractVector{<:Complex}, K::Integer; L::Integer=0)
    N = length(x)
    K = clamp(K, 1, (N - 1) ÷ 2)
    L = L > 0 ? L : min((N + 1) ÷ 2, _HANKEL_LMAX)
    L = clamp(L, K + 1, N - K)
    M = N - L + 1
    H = Matrix{eltype(x)}(undef, L, M)
    for j in 1:M, i in 1:L
        H[i, j] = x[i + j - 1]
    end
    F = svd(H)
    Uk = @view F.U[:, 1:K]
    Z = Uk[1:(end - 1), :] \ Uk[2:end, :]
    return eigvals(Z)
end

"least-squares complex amplitudes `cₖ` such that `x[n] ≈ Σ cₖ·polesₖ^(n-1)`"
function _amplitudes(x::AbstractVector{<:Complex}, poles::AbstractVector{<:Complex})
    N = length(x)
    V = Matrix{eltype(x)}(undef, N, length(poles))
    for k in eachindex(poles)
        p = poles[k]
        acc = one(p)
        for n in 1:N
            V[n, k] = acc
            acc *= p
        end
    end
    return V \ collect(x)
end

"evaluate `Σ cₖ·polesₖ^(n-1)` over the index range `ns`, iterating the powers"
function _evalmodel(poles, c, ns::AbstractRange)
    K = length(poles)
    T = eltype(c)
    cur = T[poles[k]^(first(ns) - 1) for k in 1:K]
    stepfac = step(ns) == 1 ? poles : T[poles[k]^step(ns) for k in 1:K]
    out = Vector{T}(undef, length(ns))
    @inbounds for i in eachindex(out)
        s = zero(T)
        for k in 1:K
            s += c[k] * cur[k]
            cur[k] *= stepfac[k]
        end
        out[i] = s
    end
    return out
end

"""
    _hsvdsubtract(x; dt, targethz, bandhz, ncoef) -> Vector{Complex}

Fit `x` with `ncoef` damped exponentials and subtract those whose frequency lies within `bandhz` of
`targethz` (both in Hz). Frequency-selective solvent/peak removal.
"""
function _hsvdsubtract(x::AbstractVector{<:Complex}; dt::Real, targethz::Real, bandhz::Real,
                       ncoef::Integer)
    poles = _hankelpoles(x, ncoef)
    c = _amplitudes(x, poles)
    freqs = angle.(poles) ./ (2π * dt)
    keep = abs.(freqs .- targethz) .<= (bandhz / 2)
    any(keep) || return collect(x)
    solvent = _evalmodel(poles[keep], c[keep], 1:length(x))
    return collect(x) .- solvent
end

# ---------------------------------------------------------------------------------------------------
# Multicomplex glue
# ---------------------------------------------------------------------------------------------------

"""
    _maptransform(A, di, mcindex, matkernel) -> Array{<:Multicomplex}

Apply the matrix kernel `matkernel` — mapping the `(L, ncols)` complex matrix of imₖ-plane vectors
(one per column) to a `(L2, ncols)` matrix, possibly with `L2 ≠ L` — to the multicomplex array `A`
along array dimension `di`, where `mcindex` is that dimension's imaginary unit. Uses the same
zero-copy reinterpretation as `_mcfft!`; a fresh output array is built so the kernel is free to change
the length along `di`.
"""
function _maptransform(A::AbstractArray{<:Multicomplex{T,N,C}}, di::Integer, mcindex::Integer,
                       matkernel) where {T,N,C}
    D = ndims(A)
    Ac = copy(A)
    tax = di + (N - 1)
    # the complex view wraps Ac's memory without referencing Ac: keep Ac alive while it is read
    wp, perm = GC.@preserve Ac begin
        w = MulticomplexNumbers.unsafe_ascomplex!(Ac, mcindex)   # transform axis at di+(N-1)
        perm = pushfirst!(deleteat!(collect(1:ndims(w)), tax), tax)  # transform axis to the front
        permutedims(w, perm), perm
    end
    L = size(wp, 1)
    rest = size(wp)[2:end]
    mat = reshape(wp, L, :)

    outmat = matkernel(mat)::AbstractMatrix
    L2 = size(outmat, 1)
    wnew = permutedims(reshape(outmat, L2, rest...), invperm(perm))

    outsize = ntuple(i -> i == di ? L2 : size(A, i), D)
    Aout = zeros(Multicomplex{T,N,C}, outsize)
    GC.@preserve Aout begin
        wov = MulticomplexNumbers.unsafe_ascomplex!(Aout, mcindex)
        wov .= wnew
        MulticomplexNumbers.unsafe_fromcomplex!(Aout, mcindex)
    end
    return Aout
end

"wrap a per-column complex-vector kernel `f` as a matrix kernel"
function _percolumn(f)
    return function (mat)
        cols = [f(view(mat, :, j)) for j in axes(mat, 2)]
        L2 = length(first(cols))
        out = Matrix{eltype(mat)}(undef, L2, size(mat, 2))
        for j in axes(mat, 2)
            @views out[:, j] .= cols[j]
        end
        return out
    end
end

# ---------------------------------------------------------------------------------------------------
# linearpredict
# ---------------------------------------------------------------------------------------------------

"""
    linearpredict(spec, [dim]; npoints=nothing, factor=2, order=:auto, direction=:forward,
                  replace=0, rootfix=true, pool=true)

Extend a time dimension by linear prediction, or rebuild its first points.

Linear prediction models each FID as a sum of decaying oscillations, using the autoregressive
coefficient method of nmrPipe and nmrglue, and extrapolates it by recurrence. Forward extension uses
the Zhu–Bax scheme: it fits backward coefficients, reflects the roots of the prediction polynomial
inside the unit circle, and builds forward coefficients from them, which keeps the extrapolation
stable.

With no dimension argument, every indirect Fourier dimension is extended. Prediction uses the FID
from t = 0; digital-filter build-up points before t = 0 are left unchanged.

# Keywords
- `npoints`: length after prediction. Defaults to `factor` times the current length.
- `factor`: length multiplier used when `npoints` is not given.
- `order`: number of prediction coefficients. `:auto` uses `clamp(td ÷ 2, 8, 32)`.
- `direction`: `:forward` appends points; `:backward` prepends them and moves the start of the time
  axis back. Backward extension rebuilds points missing at the start of a FID that begins after
  t = 0 (see [`loadfid`](@ref NMRTools.NMRIO.loadfid)).
- `replace`: if greater than 0, replace that many points from t = 0 by backward prediction, keeping
  the length unchanged. Overrides `npoints` and `direction`.
- `rootfix`: reflect unstable roots so that the extrapolation cannot grow.
- `pool`: estimate one set of coefficients from the strongest vectors and apply it to all of them.
  Vectors along one dimension share the same frequencies, so this is valid and much faster than
  fitting each vector separately (`pool=false`, as nmrglue does).

Use on time-domain data, before [`apodize`](@ref). The `:td` metadata is updated.

# Examples
```julia
fid |> linearpredict(T2Dim) |> apodize |> zerofill |> ft     # double the indirect dimension
linearpredict(fid, T1Dim; replace=2)                          # rebuild the first two points
```
"""
function linearpredict(spec::NMRData, dim; npoints=nothing, factor::Integer=2, order=:auto,
                       direction::Symbol=:forward, replace::Integer=0, rootfix::Bool=true,
                       pool::Bool=true)
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    d isa TimeDimension ||
        throw(ArgumentError("linearpredict: dimension $di ($(typeof(d))) is not a TimeDimension"))
    mcindex = metadata(d, :mcindex)
    isnothing(mcindex) &&
        throw(ArgumentError("linearpredict: dimension $di has no :mcindex"))

    n = size(spec, di)
    m = order === :auto ? clamp(n ÷ 2, 8, 32) : order
    # Digital-filter build-up points before t = 0 (direct dimension of `loadfid` data) are left
    # alone: prediction works on the FID from t = 0.
    i0 = max(round(Int, _groupdelaypoints(d)), 0) + 1
    if i0 > 1 && replace == 0 && direction === :backward
        throw(ArgumentError("linearpredict: backward extension is not possible while the " *
                            "dimension holds digital-filter build-up points before t = 0; use " *
                            "`replace` to rebuild the first points instead"))
    end
    fromzero(f) = mat -> vcat(mat[1:(i0 - 1), :], f(mat[i0:end, :]))

    if replace > 0
        kernel = fromzero(_percolumn(x -> _lprepair(x; nrepair=replace, order=m,
                                                     rootfix=rootfix)))
        newlen = n
    else
        target = isnothing(npoints) ? factor * n : npoints
        npred = target - n
        npred > 0 ||
            throw(ArgumentError("linearpredict: target length $target ≤ current length $n"))
        kernel = fromzero(mat -> _lpmatrix(mat; npred=npred, order=m, direction=direction,
                                           rootfix=rootfix, pool=pool))
        newlen = target
    end

    newdata = _maptransform(parent(spec), di, mcindex, kernel)
    prepended = (replace == 0 && direction === :backward) ? newlen - n : 0
    newdim = _extendtimedim(d, di, newlen; prepend=prepended)
    newdims = _replacedimtuple(dims(spec), di, newdim)
    newspec = _rebuilddims(spec, newdata, newdims)
    metadata(newspec, di)[:td] = newlen          # LP grows the *acquired* points
    return newspec
end

# no-argument form: extend every indirect Fourier dimension (skip the direct dimension)
function linearpredict(spec::NMRData; kw...)
    out = spec
    for di in _fourierdims(spec)
        di == 1 && continue
        out = linearpredict(out, di; kw...)
    end
    return out
end

linearpredict(dim; kw...) = spec -> linearpredict(spec, dim; kw...)
linearpredict(; kw...) = spec -> linearpredict(spec; kw...)
function linearpredict(d::AbstractDict{Symbol,<:NMRData}, args...; kw...)
    return Dict(k => linearpredict(v, args...; kw...) for (k, v) in d)
end
