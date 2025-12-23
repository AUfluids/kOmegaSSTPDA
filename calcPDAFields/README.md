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
    
    // Normalisation control
    normalise       true;              // Enable/disable normalisation
    normalisationField omega;          // Field name for normalisation
                                      // Options: "omega", "epsilon", "S" or "strain", 
                                      //          or any scalar field name available in the case
    
    // Optional: use k and epsilon for time scale calculation (if normalisationField not found)
    // useKEpsilon  false;
    // k            k;
    // epsilon      epsilon;
    epsilon         1e-10;            // Small value to prevent division by zero
}
```

### Normalisation Options

- **`normalise`**: Set to `true` to normalise fields using a time scale, `false` to output non-normalised (dimensional) fields
- **`normalisationField`**: Name of the scalar field to use for normalisation:
  - `"omega"`: Use omega field (default, falls back to strain rate if not available)
  - `"epsilon"`: Use epsilon field (falls back to strain rate if not available)
  - `"S"` or `"strain"`: Use strain rate magnitude
  - Any other scalar field name: Will use that field if available

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

## Time Scale Calculation

The time scale calculation depends on the `normalise` and `normalisationField` settings:

**If `normalise` is `true`:**
- If `normalisationField` is `"omega"` and omega field exists: `tau = 1 / (omega + epsilon)`
- If `normalisationField` is `"epsilon"` and epsilon field exists: `tau = 1 / (epsilon + epsilon)`
- If `normalisationField` is `"S"` or `"strain"`: `tau = 1 / (S + epsilon)` where S is strain rate magnitude
- If `normalisationField` is any other field name and exists: `tau = 1 / (field + epsilon)`
- If the specified field is not found, falls back to: `tau = 1 / (S + epsilon)`
- If `useKEpsilon` is `true` and both `k` and `epsilon` are available: `tau = k / (epsilon + epsilon)`

**If `normalise` is `false`:**
- `tau = 1.0` (dimensionless, no normalisation applied)
- Fields are output with their original dimensions

## Notes

- All fields are written with `AUTO_WRITE` enabled, so they will be written automatically when the mesh/time directory is written.
- The calculations match the implementation in `kOmegaSSTPDABase.C` for consistency.

