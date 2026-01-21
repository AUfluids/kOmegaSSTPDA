# ExprTK Setup Guide

This guide explains how to set up ExprTK for the kOmegaSSTPDA symbolic regression integration.

## Quick Setup

### Step 1: Download ExprTK

```bash
cd kOmegaSSTPDA_symbolic_regression
mkdir -p external
cd external
git clone https://github.com/ArashPartow/exprtk.git
cd ../..
```

### Step 2: Configure Makefile

Edit `Make/options` (or copy from `Make/options.exprtk`):

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
```

**Note**: If you're using `Make/options.exprtk`, uncomment the `EXPRTK_DIR` line and the `-DEXPRTK_AVAILABLE` flag.

### Step 3: Verify Setup

Check that `external/exprtk/exprtk.hpp` exists:

```bash
ls external/exprtk/exprtk.hpp
```

If the file exists, you're ready to compile!

### Step 4: Compile

```bash
wmake libso
```

## Verification

After compilation, when you run a simulation with `useSymbolicRegression true`, you should see:

```
Symbolic regression: Separation expression compiled successfully.
Symbolic regression: Tij2 expression compiled successfully.
```

If ExprTK is not available, you'll see warnings:

```
Warning: ExprTK not available. Expression will not be evaluated.
```

## Troubleshooting

### ExprTK Not Found

**Error**: `fatal error: exprtk.hpp: No such file or directory`

**Solution**: 
1. Verify ExprTK is downloaded: `ls external/exprtk/exprtk.hpp`
2. Check `EXPRTK_DIR` in `Make/options` points to correct path
3. Ensure `-I$(EXPRTK_DIR)` is in `EXE_INC`

### Expression Not Compiled

**Error**: `Failed to compile expression`

**Solution**:
1. Check expression syntax (must be C++ compatible)
2. Verify all variables are registered
3. Check expression string doesn't contain invalid characters
4. Look for parser error messages in log

### Performance Issues

If expression evaluation is slow:
1. Simplify expressions (shorter expressions evaluate faster)
2. Consider caching evaluated fields if expressions don't change
3. Profile to identify bottlenecks

## Alternative: Without ExprTK

If you don't want to use ExprTK, the code will:
- Still read expressions from dictionaries
- Show debug output
- Fall back to hardcoded coefficients
- Display warnings that expressions are not evaluated

This is useful for testing/debugging but expressions won't affect the simulation.

## References

- [ExprTK GitHub](https://github.com/ArashPartow/exprtk)
- [ExprTK Documentation](https://github.com/ArashPartow/exprtk#documentation)
