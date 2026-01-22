# Symbolic Regression Integration for kOmegaSSTPDA

This directory contains the symbolic regression version of the kOmegaSSTPDA turbulence model, designed to work with [PySR](https://github.com/MilesCranmer/PySR) for runtime expression evaluation without recompilation.

## Overview

The symbolic regression integration allows you to:
1. Train symbolic expressions using PySR
2. Export expressions to OpenFOAM dictionary format
3. Use expressions in the turbulence model without recompiling

## Architecture

### Approach 1: Expression Dictionary (Recommended for Simple Expressions)

Store PySR-generated expressions as strings in OpenFOAM dictionaries. The turbulence model reads and evaluates these expressions at runtime.

**Pros:**
- No external dependencies
- Simple to implement
- Fast iteration

**Cons:**
- Limited to simple mathematical expressions
- Requires custom parser for complex expressions

### Approach 2: ExprTK Integration (Recommended for Complex Expressions)

Use the [ExprTK](https://github.com/ArashPartow/exprtk) header-only library for fast expression parsing and evaluation.

**Pros:**
- Fast evaluation
- Supports complex mathematical expressions
- Header-only (no linking required)

**Cons:**
- Requires downloading ExprTK header
- Slightly more complex integration

## Setup

### Option A: Using Expression Dictionaries (Current Implementation)

1. Train your model with PySR:
```python
import pysr
from pysr import PySRRegressor

model = PySRRegressor(
    niterations=100,
    binary_operators=["+", "-", "*", "/"],
    unary_operators=["exp", "log", "sqrt"],
)

model.fit(X, y)
```

2. Export expression to OpenFOAM format:
```bash
python pysr_to_openfoam.py pysr_model.json -o separationExpression.dict -t separation
```

3. Use in turbulenceProperties:
```cpp
kOmegaSSTPDABaseCoeffs
{
    // ... other coefficients ...
    
    // Symbolic regression expression
    separationExpression
    {
        type            expression;
        expression      "C0 + C1*I1 + C2*I2";
        // Variables will be registered automatically
    }
}
```

### Option B: Using ExprTK (Requires Setup)

1. Download ExprTK:
```bash
cd kOmegaSSTPDA_symbolic_regression
mkdir -p external
cd external
git clone https://github.com/ArashPartow/exprtk.git
```

2. Update `Make/options`:
```cpp
EXE_INC = \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -Iexternal/exprtk \
    ...
```

3. Include ExprTK in `expressionParser.H`:
```cpp
#include "exprtk.hpp"
```

4. Implement full parser functionality in `expressionParser.C`

## Usage

### Training with PySR

Example Python script for training separation correction:

```python
import numpy as np
import pysr
from pysr import PySRRegressor

# Load your training data
# X should contain invariants: I1, I2, I3, I4, I5
# y should contain target separation factor
X = np.load('invariants.npy')  # Shape: (n_samples, 5)
y = np.load('separation_factor.npy')  # Shape: (n_samples,)

# Train model
model = PySRRegressor(
    niterations=100,
    binary_operators=["+", "-", "*", "/"],
    unary_operators=["exp", "log", "sqrt", "abs"],
    complexity_of_operators={"+": 1, "-": 1, "*": 2, "/": 2},
    maxsize=20,
    populations=10,
)

model.fit(X, y)

# Get best equation
best_equation = model.get_best()
print(f"Best equation: {best_equation}")

# Export to JSON
import json
with open('pysr_separation_model.json', 'w') as f:
    json.dump({
        'equation': str(best_equation),
        'complexity': model.get_best()['complexity'],
        'loss': model.get_best()['loss'],
    }, f, indent=2)
```

### Converting PySR Output

```bash
# Convert separation model
python pysr_to_openfoam.py pysr_separation_model.json \
    -o constant/separationExpression.dict \
    -t separation

# Convert anisotropy model for Tij2
python pysr_to_openfoam.py pysr_anisotropy_Tij2.json \
    -o constant/anisotropyExpression_Tij2.dict \
    -t anisotropy \
    --tensor-index 2
```

### Using in OpenFOAM Case

1. Place expression dictionaries in `constant/` directory
2. Update `turbulenceProperties`:

```cpp
RAS
{
    RASModel        kOmegaSSTPDA;

    turbulence      on;
    printCoeffs     on;

    kOmegaSSTPDACoeffs
    {
        // Standard coefficients
        alphaK1         0.85;
        alphaK2         1.0;
        // ... etc ...

        // Enable symbolic regression
        useSymbolicRegression    true;
        
        // Path to expression dictionaries
        separationExpressionDict "constant/separationExpression.dict";
        anisotropyExpressionDicts
        {
            Tij2    "constant/anisotropyExpression_Tij2.dict";
            Tij3    "constant/anisotropyExpression_Tij3.dict";
            // ... etc ...
        }
    }
}
```

## Implementation Details

### Expression Parser

The `expressionParser` class provides:
- Expression compilation from strings
- Variable registration (fields, scalars, constants)
- Field-wise evaluation
- Dictionary I/O

### Variable Mapping

PySR variables are mapped to OpenFOAM fields:
- `I1` → `(I1 - I1_mean)/I1_std` (normalised)
- `I2` → `(I2 - I2_mean)/I2_std`
- etc.

Constants (means, stds) are registered automatically from the dictionary.

### Performance Considerations

- Expressions are compiled once at model initialization
- Evaluation happens field-wise (cell-by-cell)
- For large meshes, consider caching evaluated fields

## Future Improvements

1. **Full ExprTK Integration**: Complete the ExprTK-based parser for maximum performance
2. **Expression Caching**: Cache evaluated expressions for repeated use
3. **Gradient Computation**: Automatic differentiation for adjoint methods
4. **Multi-threading**: Parallel evaluation across cells
5. **Expression Validation**: Check expression validity before compilation

## References

- [PySR GitHub](https://github.com/MilesCranmer/PySR)
- [ExprTK GitHub](https://github.com/ArashPartow/exprtk)
- [OpenFOAM Documentation](https://www.openfoam.com/documentation/)

## Troubleshooting

### Expression Not Evaluating

- Check expression syntax matches ExprTK/C++ format
- Verify all variables are registered
- Check dictionary file paths are correct

### Performance Issues

- Consider simplifying expressions
- Use ExprTK for better performance
- Profile expression evaluation

### Compilation Errors

- Ensure ExprTK header is in include path (if using)
- Check C++ standard compatibility (C++11 or later)
