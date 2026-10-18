"""
    WindowFunction

Abstract type to represent apodization functions.

Window functions are represented by subtypes of the abstract type `WindowFunction`,
each of which contain appropriate parameters to specify the particular function
applied. In addition, the acquisition time `tmax` is also stored (calculated at
the point the window function is applied, i.e. after linear prediction but before
zero filling).
"""
abstract type WindowFunction end

Base.Broadcast.broadcastable(w::WindowFunction) = Ref(w)

"""
    NullWindow(tmax)

No apodization applied. Acquisition time is `tmax`.
"""
struct NullWindow <: WindowFunction
    tmax::Float64
    NullWindow(tmax=Inf) = new(tmax)
end

"""
    UnknownWindow(tmax)

Unknown apodization applied. Acquisition time is `tmax`.
"""
struct UnknownWindow <: WindowFunction
    tmax::Float64
    UnknownWindow(tmax=Inf) = new(tmax)
end

"""
    ExponentialWindow(lb, tmax)

Exponential window function, with a line broadening of `lb` Hz. Acquisition time is `tmax`.
"""
struct ExponentialWindow <: WindowFunction
    lb::Float64
    tmax::Float64
    ExponentialWindow(lb=0.0, tmax=Inf) = new(lb, tmax)
end

"""
    SineWindow(offset, endpoint, power, tmax)

Abstract window function representing multiplication by sine/cosine functions.
Acquisition time is `tmax`.
```math
\\left[\\sin\\left(
    \\pi\\cdot\\mathrm{offset} +
    \\frac{\\left(\\mathrm{end} - \\mathrm{offset}\\right)\\pi t}{\\mathrm{tmax}}
    \\right)\\right]^\\mathrm{power}
```

Specialises to `CosWindow`, `Cos2Window` or `GeneralSineWindow`.

# Arguments
- `offset`: initial value is ``\\sin(\\mathrm{offset}\\cdot\\pi)`` (0 to 1)
- `endpoint`: initial value is ``\\sin(\\mathrm{endpoint}\\cdot\\pi)`` (0 to 1)
- `pow`: sine exponent
"""
abstract type SineWindow <: WindowFunction end

function SineWindow(offset=0.0, endpoint=1.0, power=1.0, tmax=Inf)
    if power ≈ 1.0 && offset ≈ 0.5
        return CosWindow(endpoint * tmax)
    elseif power ≈ 2.0 && offset ≈ 0.5
        return Cos2Window(endpoint * tmax)
    else
        return GeneralSineWindow(offset, endpoint, power, tmax)
    end
end

struct GeneralSineWindow <: SineWindow
    offset::Float64
    endpoint::Float64
    power::Float64
    tmax::Float64
    function GeneralSineWindow(offset=0.0, endpoint=1.0, power=1.0, tmax=Inf)
        return new(offset, endpoint, power, tmax)
    end
end

"""
    CosWindow(tmax)

Apodization by a pure cosine function. Acquisition time is `tmax`.

See also [`Cos2Window`](@ref), [`SineWindow`](@ref).
"""
struct CosWindow <: SineWindow
    tmax::Float64
    CosWindow(tmax=Inf) = new(tmax)
end

"""
    Cos2Window(tmax=Inf)

Apodization by a cosine-squared function, `cos²(πt / 2tmax)`, falling from 1 at t = 0 to 0 at
`tmax`. When `tmax` is left at `Inf`, [`apodize`](@ref) sets it to the acquisition time.
Equivalent to TopSpin's `QSINE` with `SSB = 2`.

`Cos²Window` is an alias, kept for compatibility.

See also [`CosWindow`](@ref), [`SineWindow`](@ref).
"""
struct Cos2Window <: SineWindow
    tmax::Float64
    Cos2Window(tmax=Inf) = new(tmax)
end

"Alias of [`Cos2Window`](@ref), kept for compatibility."
const Cos²Window = Cos2Window

"""
    GaussWindow(expHz, gaussHz, center, tmax)

Abstract representation of Lorentz-to-Gauss window functions, applying an inverse
exponential of `expHz` Hz, and a gaussian broadening of `gaussHz` Hz, with maximum
at `center` (between 0 and 1). Acquisition time is `tmax`.

Specialises to `LorentzToGaussWindow` when `center` is zero, otherwise `GeneralGaussWindow`.
"""
abstract type GaussWindow <: WindowFunction end

function GaussWindow(expHz=0.0, gaussHz=0.0, center=0.0, tmax=Inf)
    if center ≈ 0.0
        return LorentzToGaussWindow(expHz, gaussHz, tmax)
    else
        return GeneralGaussWindow(expHz, gaussHz, center, tmax)
    end
end

struct GeneralGaussWindow <: GaussWindow
    expHz::Float64
    gaussHz::Float64
    center::Float64
    tmax::Float64
    function GeneralGaussWindow(expHz=0.0, gaussHz=0.0, center=0.0, tmax=Inf)
        return new(expHz, gaussHz, center, tmax)
    end
end

struct LorentzToGaussWindow <: GaussWindow
    expHz::Float64
    gaussHz::Float64
    tmax::Float64
    LorentzToGaussWindow(expHz=0.0, gaussHz=0.0, tmax=Inf) = new(expHz, gaussHz, tmax)
end

"""
    LineshapeComplexity

Abstract type to specify calculation of `RealLineshape` or `ComplexLineshape` in function calls.
"""
abstract type LineshapeComplexity end

"""
    RealLineshape

Return a real-valued lineshape when used in calculations
"""
struct RealLineshape <: LineshapeComplexity end

"""
    ComplexLineshape

Return a complex-valued lineshape when used in calculations
"""
struct ComplexLineshape <: LineshapeComplexity end

"""
    lineshape(axis, δ, R2, complexity=RealLineshape())

Return a simulated real- or complex-valued spectrum for a resonance with chemical
shift `δ` and relaxation rate `R2`, using the parameters and window function associated
with the specified axis.
"""
function lineshape(ax, δ, R2, complexity::LineshapeComplexity)
    return _lineshape(2π * hz(δ, ax), R2, 2π * hz(ax), ax[:window], complexity)
end
# default to a real return type
lineshape(ax, δ, R2) = lineshape(ax, δ, R2, RealLineshape())

"""
    _lineshape(ω, R2, ωaxis, window, complexity)

Internal function to calculate a resonance lineshape with frequency `ω` and
relaxation rate `R2`, calculated at frequencies `ωaxis` and with apodization according
to the specified window function.
"""
function _lineshape end

function _lineshape(ω, R, ωax, ::WindowFunction, ::RealLineshape)
    # generic case - return a default real-valued Lorentzian
    @. R / ((ωax - ω)^2 + R^2)
end

function _lineshape(ω, R, ωax, ::WindowFunction, ::ComplexLineshape)
    # generic case - return a default (complex-valued) Lorentzian
    @. 1 / (R + 1im * (ω - ωax))
end

# exponential functions

function _lineshape(ω, R, ωax, w::ExponentialWindow, ::ComplexLineshape)
    x = @. R + 1im * (ω - ωax) + π * w.lb
    T = w.tmax

    return @. (1 - exp(-T * x)) / x
end

function _lineshape(ω, R, ωax, w::ExponentialWindow, ::RealLineshape)
    return real(_lineshape(ω, R, ωax, w, ComplexLineshape()))
end

# cosine

function _lineshape(ω, R, ωax, w::CosWindow, ::ComplexLineshape)
    x = @. R + 1im * (ω - ωax)
    T = w.tmax
    Tx = T * x
    return @. 2 * T * (π * exp(-Tx) + 2 * Tx) / (π^2 + 4 * Tx^2)
end

function _lineshape(ω, R, ωax, w::CosWindow, ::RealLineshape)
    return real(_lineshape(ω, R, ωax, w, ComplexLineshape()))
end

function _lineshape(ω, R, ωax, w::Cos2Window, ::ComplexLineshape)
    x = @. R + 1im * (ω - ωax)
    Tx = w.tmax * x

    return @. (π^2 * (1 - exp(-Tx)) + 2 * Tx^2) / (2 * (π^2 + Tx^2) * x)
end

function _lineshape(ω, R, ωax, w::Cos2Window, ::RealLineshape)
    return real(_lineshape(ω, R, ωax, w, ComplexLineshape()))
end

"""
    apod(spec::NMRData, dimension, zerofill=true)

Return the time-domain apodization function for the specified axis,
as a vector of values.
"""
function apod(spec::AbstractNMRData, dimension, zerofill=true)
    return apod(dims(spec, dimension), zerofill)
end

function apod(ax::FrequencyDimension, zerofill=true)
    td = ax[:td]
    sw = ax[:swhz]
    window = ax[:window]

    dt = 1 / sw
    t = dt * (0:(td - 1))

    if zerofill
        tdzf = ax[:tdzf]
        w = zeros(tdzf)
    else
        w = zeros(td)
    end
    w[1:td] = apod(t, window)

    return w
end

# function apod(t, w::WindowFunction)
#     @warn "time-domain apodization function not yet defined for window $w - treating as no apodization"
#     ones(length(t))
# end

apod(t, ::NullWindow) = ones(length(t))
apod(t, w::ExponentialWindow) = exp.(-π * w.lb * t)
function apod(t, w::GeneralSineWindow)
    @. sin(π * w.offset + π * ((w.endpoint - w.offset) * t / w.tmax))^w.power
end
apod(t, w::CosWindow) = @. cos(π / 2 * t / w.tmax)
apod(t, w::Cos2Window) = @. cos(π / 2 * t / w.tmax)^2
