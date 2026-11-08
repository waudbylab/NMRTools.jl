# Pulse program annotation schema

## Syntax – programmatic arrays

- parameter: {type, start, step}
- parameter: {type, start, end}
- type = linear | logarithmic
- ‘parameter’ must occur in ‘dimensions’ to determine the array length

### Converting from vclists to durations

- duration: {counter: array, scale: X}
  durations = X * array

## Basic information

- schema_version
- sequence_version
- title
- description
- authors
- citation
- doi
- created
- last_modified
- repository
- status = alpha | beta | stable | deprecated
- experiment_type
- features
- dimensions
- typical_nuclei
  - standard isotope notation: 1H, 13C etc.
  - array mapping to channels
  - ‘nothing’ if channel not used
- acquisition_order
- reference_pulse
  - channel
  - duration
  - power
  - angle = 90

Channel references (`f1`…`f8`) resolve against the spectrometer channel model
(`channels(spec, ...)`), which holds the nucleus, base frequency, and carrier for every
channel — not just the detected ones. As a result, an `offset` frequency list converts to
absolute ppm even when its channel's nucleus is never detected (e.g. a 13C CEST profile
read out on a 1H axis): if no detected axis matches the channel, the channel's own bf/sfo
are used as the reference.

## Experiment types and associated features

- 1d, 2d, 3d
  - watergate
- hsqc, hmqc, trosy
  - gradient_selected
  - sensitivity_enhanced
- relaxation
  - R1
  - R2
  - hahn_echo
  - perfect_echo
  - cpmg
- diffusion
  - ste
  - xste
  - led
- r1rho
  - on_resonance
  - off_resonance
- cest
  - direct_saturation
  - indirect_saturation
  - zz_saturation
- hetnoe

## Parameter blocks

- relaxation
  - model: exponential_decay | inversion_recovery | saturation_recovery
  - channel
  - duration
- diffusion
  - type = bipolar
  - coherence: {channel, coherence_order}
  - big_delta
  - little_delta
  - tau (if bipolar)
  - gmax
  - g = array
  - shape
- r1rho
  - channel
  - power
  - duration
  - offset
- cest
  - channel
  - power
  - duration
  - offset
  - reference = first | X
    X = value to indicate reference plane
- hetnoe
  - channel: the saturated channel, e.g. f1 for ¹H
  - saturated: true | false for an experiment recording one condition, or one value per
    plane, e.g. [false, true], for an interleaved series
- calibration.nutation
  - channel
  - power
  - duration
  - model = sine_modulated | cosine_modulated
  - offset = 0

