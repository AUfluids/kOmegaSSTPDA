# calcPDAFields Function Object

This function object calculates PDA (Progressive Data Augmentation) turbulence model invariants and base tensors from the velocity field. It is designed for LES simulations where the kOmegaSSTPDA library is not available.

## Compilation

To compile the function object:

```bash
cd calcPDAFields
wmake libso
```

This will create `libcalcPDAFields.so` in `$FOAM_USER_LIBBIN`.

## Usage

Add the following to your `controlDict` functions section:

```cpp
calcPDAFields
{
    type            calcPDAFields;
    libs            (calcPDAFields);
    region          region0;
    enabled         true;
    log             true;
    timeStart       0;
    timeEnd         1000;
    executeControl  timeStep;
    executeInterval 1;
    writeControl    adjustableRunTime;
    writeInterval   0.01;
    
    // Field selection
    useMean         true;              // Use mean fields (UMean, kMean, omegaMean) 
                                      // or instantaneous fields (U, k, omega)
                                      // Default: false (uses U, k, omega)
    
    // Normalisation control
    normalise       true;              // Enable/disable normalisation
                                      // When enabled, uses omega or omegaMean based on useMean flag
    
    epsilon         1e-10;            // Small value to prevent division by zero (used in omegaMin)
}
```

### Field Selection

- **`useMean`**: Set to `true` to use mean fields, `false` to use instantaneous fields:
  - `true`: Uses `UMean`, `kMean`, `omegaMean`
  - `false`: Uses `U`, `k`, `omega` (default)
  - The `U` field name can be overridden using the `U` parameter if needed

### Normalisation Options

- **`normalise`**: Set to `true` to normalise fields using a time scale, `false` to output non-normalised (dimensional) fields
- When `normalise` is `true`, the function object automatically uses:
  - `omegaMean` if `useMean` is `true`
  - `omega` if `useMean` is `false`
- If the omega field is not found, normalisation is disabled and a warning is issued

## Calculated Fields

The function object calculates and writes the following fields:

### Invariants (scalar fields):
- `I1`: First invariant
- `I2`: Second invariant  
- `I3`: Third invariant
- `I4`: Fourth invariant
- `I5`: Fifth invariant

### Base Tensors (symmetric tensor fields):
- `Tij2` through `Tij10`: Base tensors for Reynolds stress modeling

### Additional Fields:
- `Sij`: Strain rate tensor
- `Omegaij`: Rotation rate tensor
- `Rij`: Reynolds stress tensor
- `bij`: Normalised Reynolds stress anisotropy tensor

## Time Scale Calculation

The time scale calculation uses the same formula as `kOmegaSSTPDABase` for consistency:

**If `normalise` is `true` and omega field is found:**
- Uses the kOmegaSSTPDA formula: `tau = 1.0 / max(S/a1 + omegaMin, omega + omegaMin)`
  - `S` is the strain rate magnitude
  - `a1 = 0.31` (kOmegaSST model coefficient)
  - `omegaMin = 1.0e-15` (small value to prevent division by zero)
- This ensures the invariants and tensors are dimensionless and match the turbulence model implementation

**If `normalise` is `false`:**
- `tau = 1.0` (dimensionless, no normalisation applied)
- Fields are output with their original dimensions:
  - Invariants have dimensions `[0 0 -2 0 0 0 0]` (I1, I2) or `[0 0 -3 0 0 0 0]` (I3, I4) or `[0 0 -4 0 0 0 0]` (I5)
  - Tensors have dimensions `[0 0 -1 0 0 0 0]`

## Notes

- All fields are written with `AUTO_WRITE` enabled, so they will be written automatically when the mesh/time directory is written.
- The calculations match the implementation in `kOmegaSSTPDABase.C` for consistency.
- When normalisation is enabled, all invariants and tensors are dimensionless, matching the turbulence model behaviour.
- The function object is designed specifically for use with omega-based turbulence models (k-omega SST, kOmegaSSTPDA, etc.).

