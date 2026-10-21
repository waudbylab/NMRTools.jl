# Power levels

Power levels from acquisition parameters are represented using the [`Power`](@ref) type, which handles both Watts (W) and dB attenuation units with automatic conversion.

## Creating Power objects

Power values must be created with an explicit unit specification:

```julia
# From dB attenuation
p1 = Power(30.0, :dB)

# From Watts
p2 = Power(0.001, :W)

# Works with any numeric type
p3 = Power(20, :dB)
```

!!! note
    Calling `Power(value)` without a unit will throw an error.

## Accessing from acqus data

Power levels are automatically parsed when loading data:

```julia
spec = loadnmr("experiment")

# Get power level (returns Power object)
p = acqus(spec, :pl, 1)  # First power level
```

## Accessing values

Use `db()` and `watts()` to retrieve values in either unit:

```julia
p = Power(20.0, :dB)

db(p)     # 20.0 (dB)
watts(p)  # 0.01 (W)
```

**Conversion formulas:**

```math
\mathrm{dB} = -10 \log_{10}(\mathrm{watts})
```

```math
\mathrm{watts} = 10^{-\mathrm{dB}/10}
```

## Zero power handling

Zero watts is handled specially to avoid logarithm issues:

```julia
p_zero = Power(0.0, :W)
db(p_zero)     # 120.0 (very high attenuation)
watts(p_zero)  # ≈ 1e-12
```

This is relevant when parsing acqus files where unused power levels may be zero.

## Converting to RF field strength

The `hz()` function converts power to radiofrequency field strength in Hz using calibration data.

### Using a reference RF strength

```julia
# Known calibration: 25 kHz at 30 dB
ref_power = Power(30.0, :dB)
ref_hz = 25000.0

# Calculate RF at different power
test_power = Power(36.0, :dB)
rf = hz(test_power, ref_power, ref_hz)
```

**Conversion formula:**

```math
\mathrm{Hz} = \mathrm{Hz_{ref}} \times 10^{-\Delta\mathrm{dB}/20}
```

where $\Delta\mathrm{dB} = \mathrm{dB} - \mathrm{dB_{ref}}$.

### Using pulse calibration

```julia
# Known: 10 μs pulse gives 90° at 30 dB
ref_power = Power(30.0, :dB)
pulse_length = 10.0    # μs
flip_angle = 90.0      # degrees

# Calculate RF at different power
test_power = Power(33.0, :dB)
rf = hz(test_power, ref_power, pulse_length, flip_angle)
```

This first calculates the reference RF strength:

```math
\mathrm{Hz_{ref}} = \frac{\theta}{360 \times t}
```

where $\theta$ is the flip angle (degrees) and $t$ is the pulse length (seconds), then applies the power scaling formula above.
