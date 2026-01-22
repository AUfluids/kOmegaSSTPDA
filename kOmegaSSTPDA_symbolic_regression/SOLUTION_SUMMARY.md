# Solution Summary: PySR Integration Without Recompilation

## Problem

You want to use PySR-generated symbolic regression expressions in your turbulence model, but don't want to recompile every time you change an expression.

## Solution Architecture

The solution uses **runtime expression evaluation** with the following components:

1. **Expression Parser** (`expressionParser.H/C`): Wraps ExprTK for fast expression evaluation
2. **PySR Converter** (`pysr_to_openfoam.py`): Converts PySR output to OpenFOAM dictionary format
3. **Dictionary-Based Configuration**: Expressions stored in dictionaries, read at runtime
4. **Integration Points**: Modified turbulence model reads and evaluates expressions

## Key Benefits

✅ **No Recompilation**: Change expressions by editing dictionary files  
✅ **Fast Evaluation**: ExprTK provides optimised expression parsing  
✅ **Flexible**: Support complex mathematical expressions  
✅ **Easy Iteration**: Train → Convert → Deploy workflow  

## Implementation Status

### ✅ Completed

1. Expression parser class structure (`expressionParser.H/C`)
2. PySR to OpenFOAM converter script (`pysr_to_openfoam.py`)
3. Documentation and integration guides
4. Makefile configuration for ExprTK

### 🔄 To Complete

1. **Full ExprTK Integration**: Complete the `expressionParser.C` implementation with actual ExprTK code
2. **Turbulence Model Integration**: Modify `kOmegaSSTPDABase.C` to use expression parsers
3. **Testing**: Verify expressions work correctly with your invariants

## Next Steps

### Option 1: Complete ExprTK Integration (Recommended)

1. Download ExprTK:
```bash
cd kOmegaSSTPDA_symbolic_regression
mkdir -p external
cd external
git clone https://github.com/ArashPartow/exprtk.git
```

2. Complete `expressionParser.C` implementation:
   - Include ExprTK header: `#include "exprtk.hpp"`
   - Implement `compileExpression()` with ExprTK parser
   - Implement `evaluate()` and `evaluateField()` methods

3. Integrate into `kOmegaSSTPDABase`:
   - Add `expressionParser` members for separation and anisotropy
   - Read expressions from dictionaries in `read()` function
   - Replace hardcoded calculations with expression evaluation in `correct()`

### Option 2: Simple Dictionary-Based Approach (Quick Start)

For immediate use without ExprTK:

1. Store expressions as strings in dictionaries
2. Use OpenFOAM's `stringOps` or create a simple evaluator
3. Parse and evaluate expressions using basic string operations

This is simpler but less performant and limited in expression complexity.

## Code Modifications Needed

### In `kOmegaSSTPDABase.H`:

```cpp
// Add expression parser members
expressionParser separationExpression_;
std::map<label, expressionParser> anisotropyExpressions_;  // Map tensor index to parser
```

### In `kOmegaSSTPDABase.C` constructor:

```cpp
// Read expression dictionaries if symbolic regression is enabled
if (useSymbolicRegression_)
{
    // Read separation expression
    if (dict.found("separationExpressionDict"))
    {
        fileName exprFile(dict.lookup("separationExpressionDict"));
        dictionary exprDict(IOdictionary(IOobject(exprFile, ...)));
        separationExpression_.read(exprDict);
    }
    
    // Read anisotropy expressions
    if (dict.found("anisotropyExpressionDicts"))
    {
        dictionary anisotropyDicts(dict.subDict("anisotropyExpressionDicts"));
        // Read each tensor expression
    }
}
```

### In `kOmegaSSTPDABase.C::correct()`:

```cpp
// Replace hardcoded calculation:
// alpha_S_ = C0_ + C1_*(I1_ - I1_mean_separation.value())/I1_std_separation.value() + ...

// With expression evaluation:
if (useSymbolicRegression_ && separationExpression_.isValid())
{
    // Register variables
    separationExpression_.registerFieldVariable("I1", I1_);
    separationExpression_.registerFieldVariable("I2", I2_);
    // ... etc ...
    
    // Evaluate expression
    separationExpression_.evaluateField(alpha_S_);
}
else
{
    // Fall back to hardcoded calculation
    alpha_S_ = C0_ + C1_*(I1_ - I1_mean_separation.value())/I1_std_separation.value() + ...;
}
```

## Example Workflow

1. **Train Model**:
```python
model = PySRRegressor(...)
model.fit(X, y)
# Save to JSON
```

2. **Convert Expression**:
```bash
python pysr_to_openfoam.py model.json -o constant/expression.dict
```

3. **Configure Case**:
```cpp
kOmegaSSTPDACoeffs
{
    useSymbolicRegression true;
    separationExpressionDict "constant/expression.dict";
}
```

4. **Run Simulation**:
```bash
# No recompilation needed!
simpleFoam
```

5. **Update Expression**:
```bash
# Modify expression.dict or regenerate from PySR
# Re-run simulation - still no recompilation!
```

## Performance Considerations

- **ExprTK**: ~10-100x faster than simple parsers
- **Caching**: Cache evaluated fields if expressions don't change
- **Compilation**: Expressions compiled once at initialization
- **Parallelization**: ExprTK supports multi-threaded evaluation

## References

- [PySR Documentation](https://github.com/MilesCranmer/PySR)
- [ExprTK Documentation](https://github.com/ArashPartow/exprtk)
- [OpenFOAM Expression Evaluation](https://www.openfoam.com/documentation/)

## Support

For issues or questions:
1. Check `README_SYMBOLIC_REGRESSION.md` for detailed documentation
2. Review `INTEGRATION_GUIDE.md` for step-by-step instructions
3. Examine example expressions in `constant/` directory
