# Harmonic inversion (Hankel-SVD) core ############################################################
#
# Models a 1-D complex signal as a sum of K damped complex exponentials
#   x[n] ≈ Σ_k cₖ · zₖⁿ ,   zₖ = exp((-λₖ + i·2π·fₖ)·dt)
# via the state-space / matrix-pencil (ESPRIT/Kung) method: SVD of a Hankel matrix gives the signal
# subspace, whose shift-invariance yields the poles zₖ; amplitudes cₖ follow from a Vandermonde fit.
#
# This single core backs both linear prediction (extrapolate with the fitted model) and HSVD solvent
# suppression (reconstruct and subtract only the components whose frequency falls in a chosen band).
# The kernels operate on plain complex vectors so they can be tested in isolation; multicomplex data
# is fed in one imₖ-plane vector at a time by `_maptransform` (see below).

# cap the Hankel row count so the per-vector SVD stays cheap on long FIDs. `L` bounds both the
# frequency resolution (~sw/L) and the cost (SVD is O(L²·M)); the ESPRIT eigenproblem is O(K³).
const _HANKEL_LMAX = 256

"""
    _hankelpoles(x, K; L=0) -> Vector{Complex}

Estimate up to `K` signal poles `zₖ` of the complex vector `x` by ESPRIT: SVD of the Hankel matrix,
then the eigenvalues of the shift-invariance operator on the truncated left singular subspace. `L`
(Hankel rows) is capped at `_HANKEL_LMAX` by default to bound cost; the model order `K` is the main
speed lever (eigenproblem is O(K³)).
"""
function _hankelpoles(x::AbstractVector{<:Complex}, K::Integer; L::Integer=0)
    N = length(x)
    K = clamp(K, 1, (N - 1) ÷ 2)
    L = L > 0 ? L : min((N + 1) ÷ 2, _HANKEL_LMAX)
    L = clamp(L, K + 1, N - K)      # need L>K and M=N-L+1>K for the pencil
    M = N - L + 1
    H = Matrix{eltype(x)}(undef, L, M)
    for j in 1:M, i in 1:L
        H[i, j] = x[i + j - 1]
    end
    F = svd(H)
    Uk = @view F.U[:, 1:K]
    Z = Uk[1:(end - 1), :] \ Uk[2:end, :]   # least-squares shift invariance
    return eigvals(Z)
end

"""
    _amplitudes(x, poles) -> Vector{Complex}

Least-squares complex amplitudes `cₖ` such that `x[n] ≈ Σ cₖ·polesₖ^(n-1)` (Vandermonde solve).
"""
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

"reflect poles outside the unit circle to inside (keep angle) so extrapolation cannot diverge"
_stabilize(poles, ::Val{:decreasing}) = [abs(z) > 1 ? 1 / conj(z) : z for z in poles]
"reflect poles inside the unit circle to outside (for backward prediction)"
_stabilize(poles, ::Val{:increasing}) = [abs(z) < 1 ? 1 / conj(z) : z for z in poles]
_stabilize(poles, ::Val{:none}) = poles

"""
    _evalmodel(poles, c, ns) -> Vector

Evaluate the fitted model `Σ cₖ·polesₖ^(n-1)` at the sample indices in the range `ns`, iterating the
powers (`cur .*= poles`) rather than recomputing `polesₖ^(n-1)` each time.
"""
function _evalmodel(poles, c, ns::AbstractRange)
    K = length(poles)
    T = eltype(c)
    cur = T[poles[k]^(first(ns) - 1) for k in 1:K]     # powers at the first index
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
    _lpextend(x; npred, ncoef, direction, rootfix) -> Vector{Complex}

Linear-prediction extension of complex vector `x` by `npred` points, using the Hankel-SVD model with
`ncoef` poles. `direction` is `:forward` (append) or `:backward` (prepend); `rootfix` stabilises the
poles (`:decreasing` for forward, `:increasing` for backward, or `:none`).
"""
function _lpextend(x::AbstractVector{<:Complex}; npred::Integer, ncoef::Integer,
                   direction::Symbol=:forward, rootfix::Symbol=:decreasing)
    poles = _stabilize(_hankelpoles(x, ncoef), Val(rootfix))
    c = _amplitudes(x, poles)
    N = length(x)
    if direction === :forward
        ext = _evalmodel(poles, c, (N + 1):(N + npred))
        return vcat(collect(x), ext)
    elseif direction === :backward
        ext = _evalmodel(poles, c, (1 - npred):0)
        return vcat(ext, collect(x))
    else
        throw(ArgumentError("linearpredict: direction must be :forward or :backward, got $direction"))
    end
end

"""
    _lprepair(x; nrepair, ncoef, rootfix) -> Vector{Complex}

Replace the first `nrepair` points of `x` by backward prediction from the remaining points (repairs
distorted initial points without changing the length).
"""
function _lprepair(x::AbstractVector{<:Complex}; nrepair::Integer, ncoef::Integer,
                   rootfix::Symbol=:increasing)
    good = @view x[(nrepair + 1):end]
    poles = _stabilize(_hankelpoles(good, ncoef), Val(rootfix))
    c = _amplitudes(good, poles)
    repaired = _evalmodel(poles, c, (1 - nrepair):0)   # extrapolate backwards from `good`
    return vcat(repaired, collect(good))
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
    isempty(findall(keep)) && return collect(x)
    solvent = _evalmodel(poles[keep], c[keep], 1:length(x))
    return collect(x) .- solvent
end

# Multicomplex glue ###############################################################################

"""
    _maptransform(A, di, mcindex, f) -> Array{<:Multicomplex}

Apply the complex-vector kernel `f` (which may change length) to every imₖ-plane vector of the
multicomplex array `A` along array dimension `di`, where `mcindex` is that dimension's imaginary
unit. Uses the same zero-copy reinterpretation as `_mcfft!`; a fresh output array is built so `f` is
free to change the length along `di`. The other imaginary units and the other array axes ride along
as independent vectors (so, e.g., LP of an indirect dimension carries the direct-dimension quadrature
as extra realisations of the same poles).
"""
function _maptransform(A::AbstractArray{<:Multicomplex{T,N,C}}, di::Integer, mcindex::Integer,
                       f) where {T,N,C}
    D = ndims(A)
    Ac = copy(A)
    w = MulticomplexNumbers.unsafe_ascomplex!(Ac, mcindex)   # transform axis at di+(N-1)
    tax = di + (N - 1)
    nd = ndims(w)
    perm = pushfirst!(deleteat!(collect(1:nd), tax), tax)    # move transform axis to the front
    wp = permutedims(w, perm)
    L = size(wp, 1)
    rest = size(wp)[2:end]
    mat = reshape(wp, L, :)

    outcols = [f(view(mat, :, j)) for j in axes(mat, 2)]
    L2 = length(first(outcols))
    all(c -> length(c) == L2, outcols) ||
        throw(DimensionMismatch("_maptransform: kernel returned inconsistent lengths"))
    outmat = Matrix{Complex{T}}(undef, L2, size(mat, 2))
    for j in axes(mat, 2)
        @views outmat[:, j] .= outcols[j]
    end
    wnew = permutedims(reshape(outmat, L2, rest...), invperm(perm))

    outsize = ntuple(i -> i == di ? L2 : size(A, i), D)
    Aout = zeros(Multicomplex{T,N,C}, outsize)
    wov = MulticomplexNumbers.unsafe_ascomplex!(Aout, mcindex)
    wov .= wnew
    MulticomplexNumbers.unsafe_fromcomplex!(Aout, mcindex)
    return Aout
end

# linearpredict ###################################################################################

"""
    linearpredict(spec, [dim]; npoints=nothing, factor=2, ncoef=:auto,
                  direction=:forward, replace=0, rootfix=:auto)

Linear-prediction extension (or first-point repair) of the specified `TimeDimension`(s), using the
Hankel-SVD model. With no dimension argument, every *indirect* Fourier dimension is extended (forward
LP of the directly detected dimension is rarely wanted).

# Keywords
- `npoints`: target length after prediction. Defaults to `factor × td`.
- `factor`: length multiplier when `npoints` is not given (default `2`).
- `ncoef`: number of LP coefficients (model order); `:auto` ⇒ `clamp(td ÷ 4, 6, 24)`. This is the
  dominant cost lever (the ESPRIT eigenproblem is O(ncoef³)); raise it for crowded dimensions.
- `direction`: `:forward` (append points) or `:backward` (prepend).
- `replace`: if `> 0`, repair the first `replace` points by backward prediction (length unchanged;
  overrides `npoints`/`direction`).
- `rootfix`: pole stabilisation, `:auto` (`:decreasing` for forward, `:increasing` for backward),
  `:decreasing`, `:increasing`, or `:none`.

Runs on time-domain data before `apodize`; updates `:td`/`:tdzf`. Element type stays `Multicomplex`.
"""
function linearpredict(spec::NMRData, dim; npoints=nothing, factor::Integer=2, ncoef=:auto,
                       direction::Symbol=:forward, replace::Integer=0, rootfix::Symbol=:auto)
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    d isa TimeDimension ||
        throw(ArgumentError("linearpredict: dimension $di ($(typeof(d))) is not a TimeDimension"))
    mcindex = metadata(d, :mcindex)
    isnothing(mcindex) &&
        throw(ArgumentError("linearpredict: dimension $di has no :mcindex"))

    n = size(spec, di)
    nc = ncoef === :auto ? clamp(n ÷ 4, 6, 24) : ncoef

    if replace > 0
        kernel = x -> _lprepair(x; nrepair=replace, ncoef=nc,
                                rootfix=(rootfix === :auto ? :increasing : rootfix))
        newlen = n
    else
        target = isnothing(npoints) ? factor * n : npoints
        npred = target - n
        npred > 0 ||
            throw(ArgumentError("linearpredict: target length $target ≤ current length $n"))
        rf = rootfix === :auto ? (direction === :backward ? :increasing : :decreasing) : rootfix
        kernel = x -> _lpextend(x; npred=npred, ncoef=nc, direction=direction, rootfix=rf)
        newlen = target
    end

    newdata = _maptransform(parent(spec), di, mcindex, kernel)
    newdim = _extendtimedim(d, di, newlen)
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
