# Compilation Guide for kOmegaSSTPDA

## Quick Start

To compile the entire package:

```bash
./Allwclean
./Allwmake
```

This will compile:
- `kOmegaSSTPDA_symbolic_regression` (includes all incompressible features + symbolic regression)
- `kOmegaSSTPDA_compressible` (compressible version)
- `calcPDAFields` (utility for post-processing)

## Directory Structure

- **`kOmegaSSTPDA_symbolic_regression/`**: Main incompressible version with symbolic regression support
  - **Replaces** `kOmegaSSTPDA_incompressible/` (which is kept for reference)
  - Includes ExprTK support for runtime expression evaluation
  - Produces library: `libPDAIncompressibleTurbulenceModels.so`

- **`kOmegaSSTPDA_incompressible/`**: Legacy version (kept for reference, not compiled by default)
  - Use `kOmegaSSTPDA_symbolic_regression` instead

- **`kOmegaSSTPDA_compressible/`**: Compressible version
  - Produces library: `libPDACompressibleTurbulenceModels.so`

- **`calcPDAFields/`**: Utility for calculating and writing PDA-related fields
  - Produces executable: `calcPDAFields`

## ExprTK Setup

The `Allwmake` script automatically downloads ExprTK if needed. For manual setup:

```bash
cd kOmegaSSTPDA_symbolic_regression
mkdir -p external
cd external
git clone https://github.com/ArashPartow/exprtk.git
cd ../..
```

Then ensure `Make/options` includes ExprTK configuration (copy from `Make/options.exprtk`).

## Compilation Options

### Standard Compilation (with ExprTK)

```bash
./Allwclean
./Allwmake
```

### Compilation Without ExprTK

If you don't need expression evaluation:

```bash
cd kOmegaSSTPDA_symbolic_regression
# Use Make/options (without ExprTK) instead of Make/options.exprtk
wmake libso
```

### HPC Compilation

For HPC systems, use the provided script:

```bash
cd kOmegaSSTPDA_symbolic_regression
./COMPILE_HPC.sh
```

## Output Libraries

After compilation, you should find:

- `$FOAM_USER_LIBBIN/libPDAIncompressibleTurbulenceModels.so` (or system equivalent)
- `$FOAM_USER_LIBBIN/libPDACompressibleTurbulenceModels.so` (or system equivalent)
- `$FOAM_USER_APPBIN/calcPDAFields` (or system equivalent)

## Verification

To verify compilation:

```bash
# Check libraries
ls -lh $FOAM_USER_LIBBIN/libPDA*.so

# Check utility
which calcPDAFields
```

## Troubleshooting

### ExprTK Not Found

If you see errors about ExprTK:
1. Ensure ExprTK is downloaded: `ls kOmegaSSTPDA_symbolic_regression/external/exprtk/exprtk.hpp`
2. Check `Make/options` includes ExprTK configuration
3. Or compile without ExprTK (expressions won't be evaluated, will use hardcoded coefficients)

### Library Not Found

If the library isn't found:
1. Check compilation completed successfully
2. Verify `$FOAM_USER_LIBBIN` is in your library path
3. Run `wmake libso` manually in the directory

### Symbolic Regression Not Working

1. Ensure `useSymbolicRegression true;` in `turbulenceProperties`
2. Check expression dictionaries exist and are valid
3. Verify ExprTK is compiled in (check for `EXPRTK_AVAILABLE` in compilation output)
4. Check log file for error messages
