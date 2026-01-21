# Integration Guide: PySR with kOmegaSSTPDA

This guide provides step-by-step instructions for integrating PySR symbolic regression expressions into the kOmegaSSTPDA turbulence model. The integration allows you to modify expressions without recompiling the turbulence model.

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Workflow Overview](#workflow-overview)
4. [Step-by-Step Integration](#step-by-step-integration)
5. [Dictionary Format Reference](#dictionary-format-reference)
6. [Running Multiple PySR Cases (Batch Processing)](#running-multiple-pysr-cases-batch-processing)
7. [Current Implementation Status](#current-implementation-status)
8. [Troubleshooting](#troubleshooting)
9. [Quick Reference for AI Automation](#quick-reference-for-ai-automation)

## Overview

The symbolic regression integration enables:

- ✅ **Runtime Expression Evaluation**: Expressions are read from dictionaries and evaluated at runtime using ExprTK
- ✅ **No Recompilation**: Modify expressions by editing dictionary files
- ✅ **PySR Compatibility**: Direct conversion from PySR output to OpenFOAM format
- ✅ **Debug Output**: Monitor expression loading and field statistics
- ✅ **ExprTK Integration**: Full expression evaluation implemented and functional
- ✅ **Batch Processing**: Scripts and examples for processing multiple cases

## Prerequisites

1. **OpenFOAM** (tested with v2506)
2. **Python 3** with PySR installed:
   ```bash
   pip install pysr
   ```
3. **kOmegaSSTPDA** turbulence model compiled with symbolic regression support
4. **ExprTK** (required for expression evaluation):
   ```bash
   cd kOmegaSSTPDA_symbolic_regression
   mkdir -p external
   cd external
   git clone https://github.com/ArashPartow/exprtk.git
   ```
   See `SETUP_EXPRTK.md` for detailed setup instructions.
   **Note**: The compilation script (`COMPILE_HPC.sh`) automatically downloads ExprTK if not present.
5. **Test Case**: Reference implementation in `testCases/SD_test/`
6. **Python Script Location**: The `pysr_to_openfoam.py` converter script is located at:
   ```
   kOmegaSSTPDA/kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py
   ```
   For batch processing, set `KOMEGASSTPDA_ROOT` environment variable or use relative paths.

## Workflow Overview

```
┌─────────────┐
│ Train PySR  │ → Generate JSON output
└──────┬──────┘
       │
       ▼
┌──────────────────┐
│ Convert to OF    │ → Generate dictionary files
│ format (Python)  │
└──────┬───────────┘
       │
       ▼
┌──────────────────┐
│ Place in case    │ → constant/separationExpression.dict
│ directory        │   constant/anisotropyExpressions.dict
└──────┬───────────┘
       │
       ▼
┌──────────────────┐
│ Configure        │ → turbulenceProperties
│ turbulenceModel │
└──────┬───────────┘
       │
       ▼
┌──────────────────┐
│ Run Simulation   │ → Expressions loaded at runtime
└──────────────────┘
```

## Step-by-Step Integration

### Step 1: Train Your PySR Model

Create a Python script to train symbolic regression models for separation and/or anisotropy corrections:

```python
import numpy as np
import pysr
from pysr import PySRRegressor
import json

# ============================================================
# SEPARATION CORRECTION MODEL
# ============================================================
# Load training data
# X: invariants (I1, I2, I3, I4, I5) - shape (n_samples, 5)
# y: target separation factor values - shape (n_samples,)
X_separation = np.load('training_data/invariants_separation.npy')
y_separation = np.load('training_data/separation_factor.npy')

# Train separation model
model_separation = PySRRegressor(
    niterations=100,
    binary_operators=["+", "-", "*", "/"],
    unary_operators=["exp", "log", "sqrt", "abs"],
    complexity_of_operators={"+": 1, "-": 1, "*": 2, "/": 2},
    maxsize=20,
    populations=10,
    ncyclesperiteration=550,
    timeout_in_seconds=3600,
)

model_separation.fit(X_separation, y_separation)

# Get best equation
best_eq_separation = model_separation.get_best()
print(f"Best separation equation: {best_eq_separation}")

# Save to JSON
with open('pysr_separation.json', 'w') as f:
    json.dump({
        'equation': str(best_eq_separation['equation']),
        'complexity': best_eq_separation['complexity'],
        'loss': best_eq_separation['loss'],
    }, f, indent=2)

# ============================================================
# ANISOTROPY CORRECTION MODEL (for each tensor)
# ============================================================
# Train models for each tensor (Tij2, Tij3, etc.)
tensor_models = {}

for tensor_idx in [2, 3, 4, 5]:  # Add more as needed
    # Load training data for this tensor
    X_anisotropy = np.load(f'training_data/invariants_anisotropy_Tij{tensor_idx}.npy')
    y_anisotropy = np.load(f'training_data/anisotropy_Tij{tensor_idx}.npy')
    
    # Train model
    model = PySRRegressor(
        niterations=100,
        binary_operators=["+", "-", "*", "/"],
        unary_operators=["exp", "log", "sqrt", "abs"],
        complexity_of_operators={"+": 1, "-": 1, "*": 2, "/": 2},
        maxsize=20,
        populations=10,
    )
    
    model.fit(X_anisotropy, y_anisotropy)
    
    # Get best equation
    best_eq = model.get_best()
    print(f"Best Tij{tensor_idx} equation: {best_eq}")
    
    # Save to JSON
    tensor_models[tensor_idx] = {
        'equation': str(best_eq['equation']),
        'complexity': best_eq['complexity'],
        'loss': best_eq['loss'],
    }

# Save all anisotropy models
with open('pysr_anisotropy.json', 'w') as f:
    json.dump(tensor_models, f, indent=2)
```

**Note**: PySR typically outputs expressions with numeric coefficients embedded (e.g., `"-2.070 + 1.119*x1 - 0.215*x2"`). The conversion script handles this automatically.

### Step 2: Convert PySR Output to OpenFOAM Format

Use the `pysr_to_openfoam.py` script to convert PySR JSON output to OpenFOAM dictionary format.

#### 2.1: Locate the Conversion Script

The script is located in `kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py`. For HPC environments, the script can be dynamically located using the `convert_pysr.sh` helper script.

#### 2.2: Convert Separation Expression

```bash
# From your case directory or from the kOmegaSSTPDA repository root
python kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py \
    pysr_separation.json \
    -o constant/separationExpression.dict \
    -t separation
```

This generates `constant/separationExpression.dict` with:
- Expression string with normalized invariants
- Normalization constants (mean/std values)

#### 2.3: Convert Anisotropy Expressions

For consolidated format (recommended - all tensors in one file):

```bash
# Convert all anisotropy models to a single consolidated dictionary
python kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py \
    pysr_anisotropy.json \
    -o constant/anisotropyExpressions.dict \
    -t anisotropy \
    --consolidate \
    --tensor-indices 2 3 4 5
```

This generates `constant/anisotropyExpressions.dict` containing all tensor expressions in one file.

**Alternative**: Convert individual tensors (if not using consolidated format):

```bash
# For individual tensor files (not recommended)
python kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py \
    pysr_anisotropy_Tij2.json \
    -o constant/anisotropyExpression_Tij2.dict \
    -t anisotropy \
    --tensor-index 2
```

### Step 3: Place Dictionary Files in Your Case

Copy the generated dictionary files to your OpenFOAM case directory:

```bash
# Example case structure
your_case/
├── constant/
│   ├── separationExpression.dict      # ← Place here
│   ├── anisotropyExpressions.dict     # ← Place here
│   └── turbulenceProperties           # ← Configure next
├── system/
└── 0/
```

### Step 4: Configure turbulenceProperties

Edit `constant/turbulenceProperties` in your case directory:

```cpp
simulationType  RAS;

RAS
{
    RASModel        kOmegaSSTPDA;
    turbulence      on;
    printCoeffs     on;

    kOmegaSSTPDACoeffs
    {
        // ============================================================
        // SYMBOLIC REGRESSION CONFIGURATION
        // ============================================================
        // Enable symbolic regression (set to true to read expressions from dictionaries)
        useSymbolicRegression    true;  // Set to true to read expressions
        
        // Enable debug output for symbolic regression
        // Shows expressions being loaded and alpha_S/alpha_A statistics every 100 iterations
        debugSymbolicRegression  true;  // Set to true for debug output
        
        // Path to separation correction expression dictionary
        separationExpressionDict "constant/separationExpression.dict";
        
        // Path to consolidated anisotropy correction expressions dictionary
        // Contains all tensor expressions (Tij2, Tij3, Tij4, etc.) in one file
        anisotropyExpressionsDict "constant/anisotropyExpressions.dict";
        
        // Optional: Write invariants and tensors for post-processing
        writeInvariantsAndTensors true;
        
        // ============================================================
        // STANDARD k-omega-SST COEFFICIENTS
        // ============================================================
        alphaK1         0.85;
        alphaK2         1.0;
        alphaOmega1     0.5;
        alphaOmega2     0.856;
        beta1           0.075;
        beta2           0.0828;
        betaStar        0.09;
        gamma1          0.5556;
        gamma2          0.44;
        a1              0.31;
        b1              1.0;
        c1              10.0;
        F3              no;
        
        // ============================================================
        // SEPARATION CORRECTION (fallback if symbolic regression disabled)
        // ============================================================
        separationCorrection      true;
        lambda1                  18.622;
        lambda2                  4.698;
        C0                       -2.070;  // Used if symbolic regression disabled
        C1                       1.119;   // Used if symbolic regression disabled
        C2                       -0.215;  // Used if symbolic regression disabled
        C3                       0.0;
        C4                       0.0;
        C5                       0.0;
        separationRelaxation      1.0;
        
        // ============================================================
        // ANISOTROPY CORRECTION (fallback if symbolic regression disabled)
        // ============================================================
        anisotropyCorrection       true;
        anisotropyRelaxation      0.6;
        
        // Anisotropy coefficients for Tij2 (used if symbolic regression disabled)
        A0_2                      -1.584;
        A1_2                      -0.685;
        A2_2                      -0.178;
        A3_2                      0.0;
        A4_2                      0.0;
        A5_2                      0.0;
        
        // Add coefficients for other tensors (Tij3-Tij10) as needed
    }
}
```

**Key Configuration Options**:

- `useSymbolicRegression`: Set to `true` to read expressions from dictionaries
- `debugSymbolicRegression`: Set to `true` to enable debug output (shows expressions and statistics)
- `separationExpressionDict`: Path to separation expression dictionary
- `anisotropyExpressionsDict`: Path to consolidated anisotropy expressions dictionary

### Step 5: Run Your Simulation

```bash
# From your case directory
simpleFoam  # or your solver
```

**What You'll See**:

1. **At Startup** (if `useSymbolicRegression true`):
   ```
   Symbolic regression: Loaded separation expression:
       "-2.070 + 1.119*((I1 - I1_mean_separation)/I1_std_separation) - 0.215*..."
   Symbolic regression: Loaded Tij2 expression:
       "-1.584 - 0.685*((I1 - I1_mean_anisotropy)/I1_std_anisotropy) - 0.178*..."
   ```

2. **During Simulation** (if `debugSymbolicRegression true`, every 100 iterations):
   ```
   [Debug] Symbolic regression evaluation (iteration 100):
       alpha_S: min=..., max=..., mean=...
       alpha_A_2: min=..., max=..., mean=...
       Expression evaluated successfully
   ```

3. **Expression Evaluation**: Expressions are evaluated at runtime using ExprTK. If evaluation fails, the code falls back to hardcoded coefficients with a warning.

## Dictionary Format Reference

### Separation Expression Dictionary

**File**: `constant/separationExpression.dict`

```cpp
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      separationExpression;
}

// PySR-generated expression with embedded coefficients
// Variables are automatically normalized: (I - mean) / std
expression     "-2.070 + 1.119*((I1 - I1_mean_separation)/I1_std_separation) - 0.215*((I2 - I2_mean_separation)/I2_std_separation)";

// Normalization constants
variables
{
    I1_mean_separation    0.029745472322525918;
    I1_std_separation     0.01781867158784395;
    I2_mean_separation    -0.024867181279038093;
    I2_std_separation     0.01800275771769106;
    I3_mean_separation    7.935733464624863e-06;
    I3_std_separation     0.00010183240848778372;
    I4_mean_separation    -2.3092425370425096e-07;
    I4_std_separation     2.7771796282956458e-05;
    I5_mean_separation    -0.0004981212641134251;
    I5_std_separation     0.0004679466893069297;
}
```

**Key Points**:
- Expression uses normalized invariants: `(I1 - I1_mean_separation)/I1_std_separation`
- PySR coefficients are embedded directly in the expression string
- Normalization constants are provided in the `variables` section

### Anisotropy Expressions Dictionary

**File**: `constant/anisotropyExpressions.dict`

```cpp
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      anisotropyExpressions;
}

// Normalization constants (shared by all tensors)
variables
{
    I1_mean_anisotropy    0.03679851253346419;
    I1_std_anisotropy     0.016109597689079515;
    I2_mean_anisotropy    -0.03681156960681101;
    I2_std_anisotropy     0.01612305604796805;
    I3_mean_anisotropy    0.0;
    I3_std_anisotropy     1.0;
    I4_mean_anisotropy    7.606158379651854e-05;
    I4_std_anisotropy     0.00014946062409414713;
    I5_mean_anisotropy    -0.0008195433737371881;
    I5_std_anisotropy     0.00043533017048806246;
}

// Tensor expressions
tensors
{
    // Tij2 expression
    Tij2
    {
        expression     "-1.584 - 0.685*((I1 - I1_mean_anisotropy)/I1_std_anisotropy) - 0.178*((I2 - I2_mean_anisotropy)/I2_std_anisotropy)";
    }
    
    // Tij3 expression (uncomment to use)
    // Tij3
    // {
    //     expression     "your_pysr_expression_here";
    // }
    
    // Add more tensors (Tij4-Tij10) as needed
}
```

**Key Points**:
- All tensor expressions are in one file (consolidated format)
- Each tensor has its own `expression` string
- Normalization constants are shared across all tensors
- Only tensors with expressions are used (commented out tensors are ignored)

## Running Multiple PySR Cases (Batch Processing)

This section provides instructions for automating the PySR workflow across multiple cases, suitable for AI automation or batch processing.

### Prerequisites for Batch Processing

1. **Repository Structure**: Ensure you have access to the `kOmegaSSTPDA` repository with the symbolic regression module
2. **Python Environment**: PySR and required dependencies installed
3. **OpenFOAM Cases**: Multiple case directories ready for PySR integration
4. **PySR Models**: Trained PySR models (JSON files) for each case

### Batch Processing Workflow

```
┌─────────────────────┐
│ For each case:      │
│ 1. Train PySR       │ → Generate JSON files
│ 2. Convert to OF    │ → Generate dictionaries
│ 3. Configure case   │ → Update turbulenceProperties
│ 4. Run simulation   │ → Execute OpenFOAM solver
└─────────────────────┘
```

### Step-by-Step Batch Processing

#### Step 1: Organize Your Cases

Create a directory structure for batch processing:

```bash
batch_processing/
├── cases/
│   ├── case_001/
│   │   ├── pysr_separation.json
│   │   ├── pysr_anisotropy.json
│   │   └── [OpenFOAM case files]
│   ├── case_002/
│   │   ├── pysr_separation.json
│   │   ├── pysr_anisotropy.json
│   │   └── [OpenFOAM case files]
│   └── ...
├── batch_convert.sh          # Conversion script
└── batch_run.sh              # Execution script
```

#### Step 2: Batch Conversion Script

Create `batch_convert.sh` to convert all PySR outputs:

```bash
#!/bin/bash
# Batch conversion script for multiple PySR cases
# Usage: ./batch_convert.sh [KOMEGASSTPDA_ROOT]

set -e

# Get kOmegaSSTPDA root directory
if [ -n "$1" ]; then
    KOMEGASSTPDA_ROOT="$1"
elif [ -n "$KOMEGASSTPDA_ROOT" ]; then
    # Use environment variable
    :
else
    # Try to find it relative to script
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
    KOMEGASSTPDA_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
fi

# Path to converter script
CONVERTER="$KOMEGASSTPDA_ROOT/kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py"

if [ ! -f "$CONVERTER" ]; then
    echo "Error: Converter script not found at: $CONVERTER"
    echo "Please set KOMEGASSTPDA_ROOT or provide as argument"
    exit 1
fi

# Process each case directory
for CASE_DIR in cases/*/; do
    if [ ! -d "$CASE_DIR" ]; then
        continue
    fi
    
    CASE_NAME=$(basename "$CASE_DIR")
    echo "=========================================="
    echo "Processing case: $CASE_NAME"
    echo "=========================================="
    
    cd "$CASE_DIR"
    
    # Convert separation expression
    if [ -f "pysr_separation.json" ]; then
        echo "Converting separation expression..."
        python3 "$CONVERTER" \
            pysr_separation.json \
            -o constant/separationExpression.dict \
            -t separation
        echo "✓ Separation expression converted"
    else
        echo "⚠ Warning: pysr_separation.json not found, skipping"
    fi
    
    # Convert anisotropy expressions (consolidated format)
    if [ -f "pysr_anisotropy.json" ]; then
        echo "Converting anisotropy expressions..."
        # Extract tensor indices from JSON (assumes structure: {"2": {...}, "3": {...}, ...})
        TENSOR_INDICES=$(python3 -c "
import json
with open('pysr_anisotropy.json', 'r') as f:
    data = json.load(f)
    indices = [str(k) for k in data.keys() if k.isdigit()]
    print(','.join(sorted(indices)))
" 2>/dev/null || echo "2,3,4,5")
        
        python3 "$CONVERTER" \
            pysr_anisotropy.json \
            -o constant/anisotropyExpressions.dict \
            -t anisotropy \
            --consolidate \
            --tensor-indices "$TENSOR_INDICES"
        echo "✓ Anisotropy expressions converted (tensors: $TENSOR_INDICES)"
    else
        echo "⚠ Warning: pysr_anisotropy.json not found, skipping"
    fi
    
    # Verify outputs
    if [ -f "constant/separationExpression.dict" ] || [ -f "constant/anisotropyExpressions.dict" ]; then
        echo "✓ Case $CASE_NAME processed successfully"
    else
        echo "✗ Error: No expression dictionaries created for case $CASE_NAME"
    fi
    
    cd - > /dev/null
    echo ""
done

echo "=========================================="
echo "Batch conversion complete!"
echo "=========================================="
```

#### Step 3: Batch Configuration Script

Create `batch_configure.sh` to update `turbulenceProperties` for all cases:

```bash
#!/bin/bash
# Batch configuration script to enable symbolic regression in all cases
# Usage: ./batch_configure.sh

set -e

for CASE_DIR in cases/*/; do
    if [ ! -d "$CASE_DIR" ]; then
        continue
    fi
    
    CASE_NAME=$(basename "$CASE_DIR")
    TURB_PROPS="$CASE_DIR/constant/turbulenceProperties"
    
    if [ ! -f "$TURB_PROPS" ]; then
        echo "⚠ Warning: turbulenceProperties not found in $CASE_NAME, skipping"
        continue
    fi
    
    echo "Configuring $CASE_NAME..."
    
    # Check if symbolic regression is already enabled
    if grep -q "useSymbolicRegression.*true" "$TURB_PROPS"; then
        echo "  ✓ Symbolic regression already enabled"
    else
        # Add symbolic regression configuration
        # This is a simplified example - adjust based on your turbulenceProperties format
        python3 << EOF
import re

with open("$TURB_PROPS", 'r') as f:
    content = f.read()

# Check if kOmegaSSTPDACoeffs section exists
if 'kOmegaSSTPDACoeffs' not in content:
    print("  ⚠ Warning: kOmegaSSTPDACoeffs section not found")
else:
    # Add useSymbolicRegression if not present
    if 'useSymbolicRegression' not in content:
        content = re.sub(
            r'(kOmegaSSTPDACoeffs\s*\{)',
            r'\1\n        useSymbolicRegression    true;\n        debugSymbolicRegression  false;',
            content
        )
    
    # Ensure expression dictionaries are referenced
    if 'separationExpressionDict' not in content:
        content = re.sub(
            r'(useSymbolicRegression\s+true;)',
            r'\1\n        separationExpressionDict "constant/separationExpression.dict";',
            content
        )
    
    if 'anisotropyExpressionsDict' not in content:
        content = re.sub(
            r'(separationExpressionDict[^;]+;)',
            r'\1\n        anisotropyExpressionsDict "constant/anisotropyExpressions.dict";',
            content
        )
    
    with open("$TURB_PROPS", 'w') as f:
        f.write(content)
    
    print("  ✓ Configuration updated")

EOF
    fi
done

echo "Batch configuration complete!"
```

#### Step 4: Batch Execution Script

Create `batch_run.sh` to run simulations for all cases:

```bash
#!/bin/bash
# Batch execution script for multiple OpenFOAM cases
# Usage: ./batch_run.sh [solver_name]

set -e

SOLVER="${1:-simpleFoam}"

echo "=========================================="
echo "Batch execution: $SOLVER"
echo "=========================================="

for CASE_DIR in cases/*/; do
    if [ ! -d "$CASE_DIR" ]; then
        continue
    fi
    
    CASE_NAME=$(basename "$CASE_DIR")
    echo ""
    echo "=========================================="
    echo "Running case: $CASE_NAME"
    echo "=========================================="
    
    cd "$CASE_DIR"
    
    # Verify case is configured
    if [ ! -f "constant/turbulenceProperties" ]; then
        echo "⚠ Warning: turbulenceProperties not found, skipping"
        cd - > /dev/null
        continue
    fi
    
    # Check if symbolic regression is enabled
    if ! grep -q "useSymbolicRegression.*true" constant/turbulenceProperties; then
        echo "⚠ Warning: Symbolic regression not enabled, skipping"
        cd - > /dev/null
        continue
    fi
    
    # Verify expression dictionaries exist
    if [ ! -f "constant/separationExpression.dict" ] && [ ! -f "constant/anisotropyExpressions.dict" ]; then
        echo "⚠ Warning: No expression dictionaries found, skipping"
        cd - > /dev/null
        continue
    fi
    
    # Run simulation
    echo "Starting simulation..."
    $SOLVER > log.$SOLVER 2>&1
    
    if [ $? -eq 0 ]; then
        echo "✓ Case $CASE_NAME completed successfully"
    else
        echo "✗ Error: Case $CASE_NAME failed (check log.$SOLVER)"
    fi
    
    cd - > /dev/null
done

echo ""
echo "=========================================="
echo "Batch execution complete!"
echo "=========================================="
```

### Python-Based Batch Processing

For more complex automation, use Python:

```python
#!/usr/bin/env python3
"""
Batch processing script for PySR integration with kOmegaSSTPDA
Usage: python batch_process.py --cases-dir cases/ --komegasstpda-root /path/to/kOmegaSSTPDA
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


def find_converter(komegasstpda_root):
    """Find the pysr_to_openfoam.py converter script"""
    converter = Path(komegasstpda_root) / "kOmegaSSTPDA_symbolic_regression" / "pysr_to_openfoam.py"
    if not converter.exists():
        raise FileNotFoundError(f"Converter not found at: {converter}")
    return converter


def convert_pysr_output(converter, input_file, output_file, expr_type, tensor_indices=None):
    """Convert PySR JSON to OpenFOAM dictionary"""
    cmd = [
        sys.executable,
        str(converter),
        str(input_file),
        "-o", str(output_file),
        "-t", expr_type
    ]
    
    if expr_type == "anisotropy":
        if tensor_indices:
            cmd.extend(["--consolidate", "--tensor-indices", tensor_indices])
        else:
            cmd.extend(["--tensor-index", "2"])
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0, result.stdout, result.stderr


def extract_tensor_indices(json_file):
    """Extract tensor indices from PySR anisotropy JSON"""
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    indices = [k for k in data.keys() if str(k).isdigit()]
    return sorted([int(i) for i in indices])


def process_case(case_dir, converter, verbose=False):
    """Process a single case directory"""
    case_path = Path(case_dir)
    case_name = case_path.name
    
    print(f"\n{'='*50}")
    print(f"Processing case: {case_name}")
    print(f"{'='*50}")
    
    # Ensure constant directory exists
    constant_dir = case_path / "constant"
    constant_dir.mkdir(exist_ok=True)
    
    success = True
    
    # Convert separation expression
    sep_json = case_path / "pysr_separation.json"
    sep_dict = constant_dir / "separationExpression.dict"
    
    if sep_json.exists():
        print(f"  Converting separation expression...")
        ok, stdout, stderr = convert_pysr_output(
            converter, sep_json, sep_dict, "separation"
        )
        if ok:
            print(f"  ✓ Separation expression converted")
        else:
            print(f"  ✗ Error converting separation: {stderr}")
            success = False
    else:
        print(f"  ⚠ pysr_separation.json not found, skipping")
    
    # Convert anisotropy expressions
    aniso_json = case_path / "pysr_anisotropy.json"
    aniso_dict = constant_dir / "anisotropyExpressions.dict"
    
    if aniso_json.exists():
        print(f"  Converting anisotropy expressions...")
        try:
            tensor_indices = extract_tensor_indices(aniso_json)
            indices_str = ",".join(str(i) for i in tensor_indices)
            print(f"    Found tensors: {indices_str}")
            
            ok, stdout, stderr = convert_pysr_output(
                converter, aniso_json, aniso_dict, "anisotropy",
                tensor_indices=indices_str
            )
            if ok:
                print(f"  ✓ Anisotropy expressions converted")
            else:
                print(f"  ✗ Error converting anisotropy: {stderr}")
                success = False
        except Exception as e:
            print(f"  ✗ Error processing anisotropy JSON: {e}")
            success = False
    else:
        print(f"  ⚠ pysr_anisotropy.json not found, skipping")
    
    # Verify outputs
    if not sep_dict.exists() and not aniso_dict.exists():
        print(f"  ⚠ Warning: No expression dictionaries created")
        success = False
    
    return success


def main():
    parser = argparse.ArgumentParser(
        description="Batch process PySR outputs for kOmegaSSTPDA"
    )
    parser.add_argument(
        "--cases-dir",
        type=str,
        required=True,
        help="Directory containing case subdirectories"
    )
    parser.add_argument(
        "--komegasstpda-root",
        type=str,
        help="Root directory of kOmegaSSTPDA repository (or set KOMEGASSTPDA_ROOT env var)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    
    args = parser.parse_args()
    
    # Get kOmegaSSTPDA root
    komegasstpda_root = args.komegasstpda_root or Path.cwd().parent
    if "KOMEGASSTPDA_ROOT" in os.environ:
        komegasstpda_root = Path(os.environ["KOMEGASSTPDA_ROOT"])
    
    converter = find_converter(komegasstpda_root)
    print(f"Using converter: {converter}")
    
    # Process all cases
    cases_dir = Path(args.cases_dir)
    if not cases_dir.exists():
        print(f"Error: Cases directory not found: {cases_dir}")
        sys.exit(1)
    
    case_dirs = [d for d in cases_dir.iterdir() if d.is_dir()]
    
    if not case_dirs:
        print(f"Warning: No case directories found in {cases_dir}")
        sys.exit(0)
    
    print(f"\nFound {len(case_dirs)} case(s) to process")
    
    successful = 0
    failed = 0
    
    for case_dir in sorted(case_dirs):
        if process_case(case_dir, converter, args.verbose):
            successful += 1
        else:
            failed += 1
    
    print(f"\n{'='*50}")
    print(f"Batch processing complete!")
    print(f"  Successful: {successful}")
    print(f"  Failed: {failed}")
    print(f"{'='*50}")
    
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    import os
    main()
```

### Verification Checklist

After batch processing, verify each case:

```bash
#!/bin/bash
# Verification script for batch-processed cases

for CASE_DIR in cases/*/; do
    CASE_NAME=$(basename "$CASE_DIR")
    echo "Verifying $CASE_NAME..."
    
    # Check expression dictionaries
    if [ -f "$CASE_DIR/constant/separationExpression.dict" ]; then
        echo "  ✓ separationExpression.dict exists"
    else
        echo "  ✗ separationExpression.dict missing"
    fi
    
    if [ -f "$CASE_DIR/constant/anisotropyExpressions.dict" ]; then
        echo "  ✓ anisotropyExpressions.dict exists"
    else
        echo "  ✗ anisotropyExpressions.dict missing"
    fi
    
    # Check turbulenceProperties
    if grep -q "useSymbolicRegression.*true" "$CASE_DIR/constant/turbulenceProperties" 2>/dev/null; then
        echo "  ✓ Symbolic regression enabled"
    else
        echo "  ✗ Symbolic regression not enabled"
    fi
    
    # Check expression syntax (basic)
    if [ -f "$CASE_DIR/constant/separationExpression.dict" ]; then
        if grep -q "expression" "$CASE_DIR/constant/separationExpression.dict"; then
            echo "  ✓ Separation expression found"
        else
            echo "  ✗ Separation expression missing"
        fi
    fi
    
    echo ""
done
```

## Current Implementation Status

### ✅ Completed

1. **Expression Reading**: Expressions are read from dictionaries at runtime
2. **Dictionary Format**: PySR output converted to OpenFOAM dictionary format
3. **Debug Output**: Statistics and expression loading information
4. **Fallback Mechanism**: Falls back to hardcoded coefficients if expressions not found
5. **ExprTK Integration**: Expression evaluation fully implemented and functional
6. **Runtime Evaluation**: Expressions are evaluated at runtime using ExprTK parser
7. **Batch Processing Support**: Scripts and examples for processing multiple cases

### 📋 Usage Notes

- **ExprTK Required**: For expression evaluation, ExprTK must be installed and compiled with the library
- **Expression Format**: Expressions use C++/ExprTK syntax with normalized invariants
- **Debug Mode**: Enable `debugSymbolicRegression true` to monitor expression evaluation

## Troubleshooting

### Expression Not Found

**Symptoms**: Warning message about expression dictionary not found

**Solutions**:
- Verify dictionary file paths in `turbulenceProperties`
- Check files exist in `constant/` directory
- Verify file permissions (readable)
- Check for typos in dictionary file names

### Expression Not Loaded

**Symptoms**: No output about expressions being loaded at startup

**Solutions**:
- Set `useSymbolicRegression true` in `turbulenceProperties`
- Verify dictionary files contain `expression` key
- Check dictionary syntax (proper OpenFOAM format)
- Look for error messages in log file

### Debug Output Not Appearing

**Symptoms**: No debug statistics printed during simulation

**Solutions**:
- Set `debugSymbolicRegression true` in `turbulenceProperties`
- Wait for 100 iterations (debug output appears every 100 iterations)
- Check log file for output

### PySR Conversion Errors

**Symptoms**: `pysr_to_openfoam.py` script fails

**Solutions**:
- Verify PySR JSON format matches expected structure
- Check Python version (requires Python 3)
- Ensure all required Python packages installed
- Check script path is correct

### Expression Syntax Errors

**Symptoms**: Warnings about expression format

**Solutions**:
- Verify expression uses C++ compatible syntax
- Check for proper parentheses matching
- Ensure operators are supported (`+`, `-`, `*`, `/`, `^`)
- Verify variable names match expected format (`I1`, `I2`, etc.)

## Example Reference Files

Reference implementation files are available in `testCases/SD_test/constant/`:

- `turbulenceProperties` - Complete configuration example
- `separationExpression.dict` - Separation expression format
- `anisotropyExpressions.dict` - Anisotropy expressions format

## Additional Resources

- **PySR Documentation**: https://github.com/MilesCranmer/PySR
- **ExprTK Documentation**: https://github.com/ArashPartow/exprtk
- **OpenFOAM Documentation**: https://www.openfoam.com/documentation/
- **kOmegaSSTPDA Repository**: See `README_SYMBOLIC_REGRESSION.md` for detailed documentation

## Support

For issues or questions:
1. Check this guide for common solutions
2. Review `README_SYMBOLIC_REGRESSION.md` for detailed documentation
3. Examine example files in `testCases/SD_test/constant/`
4. Check simulation log for error messages

---

**Last Updated**: 2024-01-01  
**Version**: 2.0.0  
**Author**: Mario Javier Rincon, PhD (Aarhus University/Kamstrup A/S)  
**Status**: ExprTK integration complete, batch processing support added
