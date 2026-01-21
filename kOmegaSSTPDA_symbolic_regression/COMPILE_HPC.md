# HPC Compilation Guide for kOmegaSSTPDA with ExprTK

This guide provides step-by-step commands to compile the kOmegaSSTPDA turbulence model with ExprTK support on HPC systems.

## Prerequisites

1. OpenFOAM environment loaded (typically via `module load` or `source`)
2. Git available for downloading ExprTK
3. Write access to the compilation directory

## Step-by-Step Compilation

### Step 1: Navigate to the Symbolic Regression Directory

```bash
cd kOmegaSSTPDA_symbolic_regression
```

### Step 2: Download ExprTK

```bash
# Create external directory if it doesn't exist
mkdir -p external
cd external

# Download ExprTK (header-only library)
git clone https://github.com/ArashPartow/exprtk.git

# Verify download
ls exprtk/exprtk.hpp

# Return to symbolic regression directory
cd ..
```

### Step 3: Configure Makefile

**Option A: Use the provided ExprTK options file**

```bash
# Copy the ExprTK options file
cp Make/options.exprtk Make/options

# Edit Make/options to set EXPRTK_DIR (if needed)
# The default assumes: EXPRTK_DIR = ./external/exprtk
# If your path is different, edit Make/options and update EXPRTK_DIR
```

**Option B: Manually configure Make/options**

Edit `Make/options` and ensure it contains:

```makefile
# ExprTK include path
EXPRTK_DIR = ./external/exprtk

EXE_INC = \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/transportModels \
    -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/incompressible/lnInclude \
    -I./ \
    -I$(EXPRTK_DIR) \
    -DEXPRTK_AVAILABLE

LIB_LIBS = \
    -lturbulenceModels \
    -lfiniteVolume \
    -lmeshTools \
    -lincompressibleTransportModels \
    -lincompressibleTurbulenceModels \
    -lstdc++
```

### Step 4: Verify Configuration

```bash
# Check that ExprTK header exists
ls external/exprtk/exprtk.hpp

# Check Make/options contains ExprTK configuration
grep EXPRTK Make/options
```

### Step 5: Clean Previous Build (Optional but Recommended)

```bash
# Remove previous compilation artifacts
wclean libso
# Or if that doesn't work:
rm -rf Make/linux* Make/darwin* Make/lnInclude
```

### Step 6: Compile the Library

```bash
# Compile the library
wmake libso
```

### Step 7: Verify Compilation

```bash
# Check that library was created
ls $FOAM_USER_LIBBIN/libkOmegaSSTPDA*.so

# Or check in the Make directory
ls Make/*/libkOmegaSSTPDA*.so
```

## Complete Command Sequence (Copy-Paste Ready)

```bash
# Navigate to symbolic regression directory
cd kOmegaSSTPDA_symbolic_regression

# Download ExprTK
mkdir -p external
cd external
git clone https://github.com/ArashPartow/exprtk.git
cd ..

# Configure Makefile
cp Make/options.exprtk Make/options

# Clean previous build (optional)
wclean libso 2>/dev/null || rm -rf Make/linux* Make/darwin* Make/lnInclude

# Compile
wmake libso

# Verify
ls $FOAM_USER_LIBBIN/libkOmegaSSTPDA*.so 2>/dev/null || ls Make/*/libkOmegaSSTPDA*.so
```

## Troubleshooting

### Error: exprtk.hpp not found

**Solution**:
```bash
# Verify ExprTK is downloaded
ls external/exprtk/exprtk.hpp

# If missing, re-download:
cd external
rm -rf exprtk
git clone https://github.com/ArashPartow/exprtk.git
cd ..
```

### Error: EXPRTK_DIR not set

**Solution**:
```bash
# Edit Make/options and ensure EXPRTK_DIR is set:
# EXPRTK_DIR = ./external/exprtk
```

### Error: Compilation fails with ExprTK errors

**Solution**:
- Check ExprTK version compatibility
- Try updating ExprTK: `cd external/exprtk && git pull`
- Verify C++ standard compatibility (ExprTK requires C++11 or later)

### Error: Library not found after compilation

**Solution**:
```bash
# Check compilation output for errors
# Verify library location:
find . -name "libkOmegaSSTPDA*.so"

# Check OpenFOAM library path:
echo $FOAM_USER_LIBBIN
```

## Alternative: Compile Without ExprTK (Fallback Mode)

If you don't want to use ExprTK, the code will still compile but expressions won't be evaluated:

```bash
# Use standard options (without ExprTK)
# The code will read expressions but fall back to hardcoded coefficients
# No special setup needed - just compile normally
wmake libso
```

## Verification After Compilation

Run a test case and check the log for:

```
Symbolic regression: Separation expression compiled successfully.
Symbolic regression: Tij2 expression compiled successfully.
```

If you see warnings about ExprTK not being available, check your Make/options configuration.

## Notes for HPC Systems

1. **Module Systems**: Ensure OpenFOAM module is loaded:
   ```bash
   module load openfoam  # or your system's module name
   ```

2. **Path Issues**: If compilation fails due to paths, use absolute paths:
   ```makefile
   EXPRTK_DIR = /full/path/to/kOmegaSSTPDA/kOmegaSSTPDA_symbolic_regression/external/exprtk
   ```

3. **Permissions**: Ensure you have write access to:
   - The compilation directory
   - `$FOAM_USER_LIBBIN` (or equivalent)

4. **Parallel Compilation**: wmake automatically uses available cores. To limit:
   ```bash
   export WM_NCOMPPROCS=4  # Use 4 cores
   wmake libso
   ```

## Quick Reference

| Task | Command |
|------|---------|
| Download ExprTK | `git clone https://github.com/ArashPartow/exprtk.git external/exprtk` |
| Configure | `cp Make/options.exprtk Make/options` |
| Clean | `wclean libso` |
| Compile | `wmake libso` |
| Verify | `ls $FOAM_USER_LIBBIN/libkOmegaSSTPDA*.so` |
