"""
    Power{T}

A type to represent NMR pulse powers that can be described either in Watts (W)
or as dB attenuation (dB). Internally stores the dB value as type `T`.

# Constructors
- `Power(value, :dB)` - Create from dB attenuation value
- `Power(value, :W)` - Create from Watts value

# Getters
- `db(p::Power)` - Get dB attenuation value
- `watts(p::Power)` - Get Watts value

# Examples
```jldoctest
julia> p1 = Power(30.0, :dB)
Power(30.0 dB, 0.001 W)

julia> db(p1), watts(p1)
(30.0, 0.001)

julia> p2 = Power(0.5, :W);

julia> round(db(p2); digits=2)
3.01

julia> db(Power(10, :W))       # integers and other numeric types work
-10.0

julia> db(Power(0.0, :W))      # zero watts is 120 dB attenuation
120.0
```

# Conversion
Converting W to dB: -10 * log10(plW)
Converting dB to W: 10^(-dB/10)
"""
struct Power{T}
    db::T

    function Power(value::T, unit::Symbol) where {T}
        x = float(value)
        FT = typeof(x)
        if unit == :dB
            new{FT}(x)
        elseif unit == :W
            if iszero(x)
                new{FT}(FT(120))  # represent 0 W as 120 dB attenuation
            else
                db_value = -10 * log10(x)
                new{FT}(FT(db_value))
            end
        else
            throw(ArgumentError("Unit must be :dB or :W"))
        end
    end
end

# Error for Power(x) without unit
function Power(value)
    throw(ArgumentError("Power constructor requires unit specification. Use Power(value, :dB) or Power(value, :W)"))
end

# Getters
"""
    db(p::Power)

Get the power level in dB attenuation.

# Example
```julia
p = Power(20.0, :dB)
db(p)  # 20.0
```
"""
db(p::Power) = p.db

"""
    watts(p::Power)

Get the power level in Watts.

# Example
```julia
p = Power(20.0, :dB)
watts(p)  # 0.01
```
"""
watts(p::Power) = 10^(-p.db / 10)

# Pretty printing
Base.show(io::IO, p::Power) = print(io, "Power($(p.db) dB, $(watts(p)) W)")

# Array interface
Base.length(p::Power) = 1
Base.getindex(p::Power, i::Int) = (i == 1) ? p : throw(BoundsError(p, i))
Base.iterate(p::Power, state=1) = (state > 1) ? nothing : (p, state + 1)

"""
    hz(p::Power, ref_p::Power, ref_Hz)

Convert power to radiofrequency strength in Hz using a reference power and known Hz value.

Uses the relationship: Hz = ref_Hz * 10^(-ΔdB/20) where ΔdB is the power difference.

# Arguments
- `p::Power`: Power to convert
- `ref_p::Power`: Reference power with known Hz value
- `ref_Hz`: Radiofrequency strength in Hz at the reference power

# Returns
Radiofrequency strength in Hz for power `p`
"""
function hz(p::Power, ref_p::Power, ref_Hz)
    ΔdB = db(p) .- db(ref_p)
    return @. ref_Hz * 10^(-ΔdB / 20)
end

"""
    hz(p::Power, ref_p::Power, ref_pulselength, ref_pulseangle_deg)

Convert power to radiofrequency strength in Hz using reference pulse parameters.

# Arguments
- `p::Power`: Power to convert
- `ref_p::Power`: Reference power
- `ref_pulselength`: Reference pulse length in microseconds
- `ref_pulseangle_deg`: Reference pulse flip angle in degrees

# Returns
Radiofrequency strength in Hz for power `p`
"""
function hz(p::Power, ref_p::Power, ref_pulselength, ref_pulseangle_deg)
    ref_hz = ref_pulseangle_deg / (360 * ref_pulselength) # pulse length in s
    return hz(p, ref_p, ref_hz)
end