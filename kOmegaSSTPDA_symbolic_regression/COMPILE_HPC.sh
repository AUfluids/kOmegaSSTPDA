#!/bin/bash
# HPC Compilation Script for kOmegaSSTPDA with ExprTK
# Usage: ./COMPILE_HPC.sh

set -e  # Exit on error

echo "=========================================="
echo "kOmegaSSTPDA Symbolic Regression Compilation"
echo "=========================================="
echo ""

# Step 1: Navigate to symbolic regression directory
echo "Step 1: Navigating to symbolic regression directory..."
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# If we're not in the right directory, try to find it
if [ ! -f "expressionParser.H" ]; then
    echo "WARNING: Not in kOmegaSSTPDA_symbolic_regression directory"
    echo "Current directory: $(pwd)"
    echo "Looking for kOmegaSSTPDA_symbolic_regression..."
    
    # Try parent directory
    if [ -d "../kOmegaSSTPDA_symbolic_regression" ]; then
        cd ../kOmegaSSTPDA_symbolic_regression
        echo "Found in parent directory"
    # Try current directory name
    elif [ -d "kOmegaSSTPDA_symbolic_regression" ]; then
        cd kOmegaSSTPDA_symbolic_regression
        echo "Found in current directory"
    else
        echo "ERROR: Could not find kOmegaSSTPDA_symbolic_regression directory"
        echo "Please run this script from the kOmegaSSTPDA_symbolic_regression directory"
        exit 1
    fi
fi

echo "Working directory: $(pwd)"

# Step 2: Download ExprTK
echo ""
echo "Step 2: Downloading ExprTK..."
if [ ! -d "external/exprtk" ]; then
    mkdir -p external
    cd external
    echo "Cloning ExprTK repository..."
    git clone https://github.com/ArashPartow/exprtk.git
    cd ..
    echo "✓ ExprTK downloaded successfully"
else
    echo "✓ ExprTK already exists, skipping download"
fi

# Verify ExprTK header exists
if [ ! -f "external/exprtk/exprtk.hpp" ]; then
    echo "ERROR: ExprTK header not found at external/exprtk/exprtk.hpp"
    exit 1
fi

# Step 3: Configure Makefile
echo ""
echo "Step 3: Configuring Makefile..."
if [ -f "Make/options.exprtk" ]; then
    cp Make/options.exprtk Make/options
    echo "✓ Copied Make/options.exprtk to Make/options"
    
    # Uncomment EXPRTK_DIR if commented
    sed -i 's/^# EXPRTK_DIR =/EXPRTK_DIR =/' Make/options || \
    sed -i.bak 's/^# EXPRTK_DIR =/EXPRTK_DIR =/' Make/options
    
    echo "✓ Configured EXPRTK_DIR in Make/options"
else
    echo "WARNING: Make/options.exprtk not found. Using existing Make/options"
fi

# Verify EXPRTK_DIR is set
if ! grep -q "EXPRTK_DIR" Make/options || grep -q "^#.*EXPRTK_DIR" Make/options; then
    echo "WARNING: EXPRTK_DIR may not be set correctly in Make/options"
    echo "Please verify Make/options contains:"
    echo "  EXPRTK_DIR = ./external/exprtk"
    echo "  And -I\$(EXPRTK_DIR) -DEXPRTK_AVAILABLE in EXE_INC"
fi

# Step 4: Clean previous build
echo ""
echo "Step 4: Cleaning previous build..."
wclean libso 2>/dev/null || {
    echo "wclean libso failed, trying manual cleanup..."
    rm -rf Make/linux* Make/darwin* Make/lnInclude 2>/dev/null || true
}
echo "✓ Clean completed"

# Step 5: Compile
echo ""
echo "Step 5: Compiling library..."
echo "This may take several minutes..."
wmake libso

# Step 6: Verify compilation
echo ""
echo "Step 6: Verifying compilation..."
if [ -f "$FOAM_USER_LIBBIN/libPDAIncompressibleTurbulenceModels.so" ]; then
    echo "✓ Library found at: $FOAM_USER_LIBBIN/libPDAIncompressibleTurbulenceModels.so"
    ls -lh "$FOAM_USER_LIBBIN/libPDAIncompressibleTurbulenceModels.so"
elif [ -f "Make/"*"/libPDAIncompressibleTurbulenceModels.so" ]; then
    LIB_FILE=$(find Make -name "libPDAIncompressibleTurbulenceModels.so" | head -1)
    echo "✓ Library found at: $LIB_FILE"
    ls -lh "$LIB_FILE"
else
    echo "WARNING: Library file not found. Check compilation output for errors."
    exit 1
fi

echo ""
echo "=========================================="
echo "Compilation completed successfully!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Verify expressions are configured in your case's turbulenceProperties"
echo "2. Run your simulation"
echo "3. Check log for: 'Symbolic regression: ... expression compiled successfully'"
echo ""
