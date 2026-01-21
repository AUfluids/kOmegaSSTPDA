#!/bin/bash
#
# ===============================================================================
# FILE: debug_symbolic_regression.sh
# ===============================================================================
#
# DESCRIPTION:
#     Run and debug the kOmegaSSTPDA model with symbolic regression
#     This script runs the simulation and filters output for symbolic regression
#     debug information
#
# AUTHOR INFORMATION:
#     Name:           Mario Javier Rincon, PhD
#     Title:          Research Scientist
#     Affiliation:    Aarhus University/Kamstrup A/S
#     Email:          mjrp@mpe.au.dk/mjrp@kamstrup.com
#
# METADATA:
#     Created:        2024-01-01
#     Last Modified:  2024-01-01
#     Version:        1.0.0
#
# LICENSE & COPYRIGHT:
#     Copyright:      © 2024 M. J. Rincon. All rights reserved.
# ===============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "================================================================================"
echo "kOmegaSSTPDA Symbolic Regression Debug Script"
echo "================================================================================"
echo ""

# Check if turbulenceProperties exists
if [ ! -f "constant/turbulenceProperties" ]; then
    echo -e "${RED}ERROR: constant/turbulenceProperties not found!${NC}"
    exit 1
fi

# Check if symbolic regression is enabled
if grep -q "useSymbolicRegression.*true" constant/turbulenceProperties; then
    echo -e "${GREEN}✓ Symbolic regression is enabled${NC}"
else
    echo -e "${YELLOW}⚠ Warning: Symbolic regression is not enabled in turbulenceProperties${NC}"
    echo "  Set 'useSymbolicRegression true;' in kOmegaSSTPDACoeffs"
fi

# Check if debug is enabled
if grep -q "debugSymbolicRegression.*true" constant/turbulenceProperties; then
    echo -e "${GREEN}✓ Debug output is enabled${NC}"
else
    echo -e "${YELLOW}⚠ Warning: Debug output is not enabled${NC}"
    echo "  Set 'debugSymbolicRegression true;' in kOmegaSSTPDACoeffs"
fi

# Check expression dictionaries and display expressions
echo ""
echo "Checking expression dictionaries..."
echo ""

# Check separation expression
if [ -f "constant/separationExpression.dict" ]; then
    echo -e "${GREEN}✓ Found: constant/separationExpression.dict${NC}"
    if grep -q "expression" constant/separationExpression.dict; then
        echo "  Expression found in dictionary:"
        # Extract expression - look for line with "expression" followed by quoted string
        # Format: expression     "-2.070 + 1.119*I1 - 0.215*I2";
        SEP_EXPR=$(grep "expression" constant/separationExpression.dict | grep -v "^[[:space:]]*//" | sed 's/.*expression[[:space:]]*//' | sed 's/^["'\''"]//' | sed 's/["'\''"];$//' | sed 's/;$//' | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//')
        
        # Alternative: extract quoted string after expression keyword using awk
        if [ -z "$SEP_EXPR" ] || echo "$SEP_EXPR" | grep -q "^//"; then
            SEP_EXPR=$(awk '/expression/ && !/^[[:space:]]*\/\// {match($0, /["'\''"]([^"'\''"]+)["'\''"]/, arr); if (arr[1] != "") print arr[1]}' constant/separationExpression.dict | head -1)
        fi
        
        # Check if expression contains operators or numbers (validate it's a real expression)
        if [ -n "$SEP_EXPR" ] && ! echo "$SEP_EXPR" | grep -q "^//"; then
            # Check for operators or numbers (put - at end of character class to avoid range interpretation)
            if echo "$SEP_EXPR" | grep -qE '[+*/\-]|[0-9]'; then
                echo -e "  ${BLUE}$SEP_EXPR${NC}"
            else
                echo -e "  ${YELLOW}(Could not extract expression - check file format)${NC}"
            fi
        else
            echo -e "  ${YELLOW}(Could not extract expression - check file format)${NC}"
        fi
    else
        echo -e "${RED}  ERROR: No 'expression' key found!${NC}"
    fi
else
    echo -e "${RED}✗ Missing: constant/separationExpression.dict${NC}"
fi

echo ""

# Check anisotropy expressions
if [ -f "constant/anisotropyExpressions.dict" ]; then
    echo -e "${GREEN}✓ Found: constant/anisotropyExpressions.dict${NC}"
    if grep -q "tensors" constant/anisotropyExpressions.dict; then
        echo "  Tensors section found"
        echo ""
        echo "  Loaded anisotropy expressions:"
        
        # Extract expressions for each tensor
        for tensor in Tij2 Tij3 Tij4 Tij5 Tij6 Tij7 Tij8 Tij9 Tij10; do
            if grep -q "$tensor" constant/anisotropyExpressions.dict; then
                # Extract expression for this tensor - look for quoted string after "expression" within tensor block
                # Format: expression     "-1.584 - 0.685*I1 - 0.178*I2";
                TENSOR_EXPR=$(sed -n "/$tensor/,/^[[:space:]]*}/p" constant/anisotropyExpressions.dict | grep "expression" | grep -v "^[[:space:]]*//" | sed 's/.*expression[[:space:]]*//' | sed 's/^["'\''"]//' | sed 's/["'\''"];$//' | sed 's/;$//' | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//' | head -1)
                
                # Alternative: extract quoted string using awk
                if [ -z "$TENSOR_EXPR" ] || echo "$TENSOR_EXPR" | grep -q "^//"; then
                    TENSOR_EXPR=$(sed -n "/$tensor/,/^[[:space:]]*}/p" constant/anisotropyExpressions.dict | awk '/expression/ && !/^[[:space:]]*\/\// {match($0, /["'\''"]([^"'\''"]+)["'\''"]/, arr); if (arr[1] != "") print arr[1]}' | head -1)
                fi
                
                # Check if expression contains operators or numbers (validate it's a real expression)
                if [ -n "$TENSOR_EXPR" ] && ! echo "$TENSOR_EXPR" | grep -q "^//" && [ "$TENSOR_EXPR" != "}" ]; then
                    # Check for operators or numbers (put - at end of character class to avoid range interpretation)
                    if echo "$TENSOR_EXPR" | grep -qE '[+*/\-]|[0-9]'; then
                        echo -e "    ${BLUE}$tensor: $TENSOR_EXPR${NC}"
                    fi
                fi
            fi
        done
    else
        echo -e "${RED}  ERROR: No 'tensors' section found!${NC}"
    fi
else
    echo -e "${RED}✗ Missing: constant/anisotropyExpressions.dict${NC}"
fi

# Determine solver
SOLVER="simpleFoam"
if [ -f "system/controlDict" ]; then
    SOLVER=$(grep "^application" system/controlDict | sed 's/.*application[[:space:]]*\([^;]*\);/\1/' | tr -d ' ')
    if [ -z "$SOLVER" ]; then
        SOLVER="simpleFoam"
    fi
fi

echo ""
echo "================================================================================"
echo "Running simulation: $SOLVER"
echo "================================================================================"
echo ""
echo "This will run the simulation and filter output for symbolic regression messages."
echo "Look for:"
echo "  - 'Symbolic regression enabled' messages at startup"
echo "  - 'Loaded ... expression' messages"
echo "  - '[Debug] Symbolic Regression' messages every 100 iterations"
echo ""
echo "Press Ctrl+C to stop the simulation"
echo ""

# Create log file
LOG_FILE="log.${SOLVER}.symbolic_regression"
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
FULL_LOG="log.${SOLVER}.${TIMESTAMP}"

echo "Full log will be written to: $FULL_LOG"
echo "Filtered symbolic regression log: $LOG_FILE"
echo ""

# Run solver and redirect output to log file
echo "Starting simulation..."
echo "Running $SOLVER (output redirected to log file)..."
echo ""

# Run solver and redirect all output to log file
$SOLVER > "$FULL_LOG" 2>&1
SOLVER_EXIT_CODE=$?

echo ""
echo "Simulation finished with exit code: $SOLVER_EXIT_CODE"

# Extract and analyze symbolic regression messages
echo ""
echo "================================================================================"
echo "Analyzing symbolic regression output..."
echo "================================================================================"

# Extract all symbolic regression related messages
grep -E "(Symbolic regression|\[CONSTRUCTOR\].*SYMBOLIC|\[CONSTRUCTOR\].*Loaded|\[CONSTRUCTOR\].*expression|\[CONSTRUCTOR\].*ERROR|\[CONSTRUCTOR\].*WARNING|\[DEBUG\].*symbolic|\[DEBUG\].*alpha_S|\[DEBUG\].*alpha_A|\[DEBUG\].*Using symbolic|expression compiled|Failed to compile|useSymbolicRegression|debugSymbolicRegression|Attempting to read|Separation correction check)" "$FULL_LOG" > "symbolic_regression_messages.txt" || true

# Display expressions that were loaded from the log
echo ""
echo "Expressions loaded during simulation:"
echo "======================================"

# Extract separation expression from log
if grep -q "Loaded separation expression" "$FULL_LOG"; then
    echo ""
    echo -e "${GREEN}Separation Expression:${NC}"
    # The expression is on the line after "Loaded separation expression:"
    # Look for lines that contain the expression pattern (contains operators and numbers)
    SEP_EXPR_LOG=$(grep -A 5 "Loaded separation expression" "$FULL_LOG" | grep -E '[+\-*/].*[0-9]|[0-9].*[+\-*/]' | head -1 | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//' | sed 's/^["'\''"]//' | sed 's/["'\''"]$//')
    
    # Alternative: get any non-empty line after "Loaded separation expression" that looks like an expression
    if [ -z "$SEP_EXPR_LOG" ] || [ "$SEP_EXPR_LOG" = ":" ]; then
        SEP_EXPR_LOG=$(grep -A 5 "Loaded separation expression" "$FULL_LOG" | grep -v "Loaded separation expression" | grep -v "^$" | grep -v "\[CONSTRUCTOR\]" | grep -E '[+\-*/]|[0-9]' | head -1 | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//' | sed 's/^["'\''"]//' | sed 's/["'\''"]$//')
    fi
    
    if [ -n "$SEP_EXPR_LOG" ] && [ "$SEP_EXPR_LOG" != ":" ]; then
        echo -e "  ${BLUE}$SEP_EXPR_LOG${NC}"
    else
        echo -e "  ${YELLOW}(Expression string not found in log - check log file directly)${NC}"
    fi
fi

# Extract anisotropy expressions from log
echo ""
echo -e "${GREEN}Anisotropy Expressions:${NC}"
for tensor in Tij2 Tij3 Tij4 Tij5 Tij6 Tij7 Tij8 Tij9 Tij10; do
    if grep -q "Loaded $tensor expression" "$FULL_LOG"; then
        # Similar extraction logic for anisotropy expressions
        TENSOR_EXPR_LOG=$(grep -A 5 "Loaded $tensor expression" "$FULL_LOG" | grep -E '[+\-*/].*[0-9]|[0-9].*[+\-*/]' | head -1 | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//' | sed 's/^["'\''"]//' | sed 's/["'\''"]$//')
        
        if [ -z "$TENSOR_EXPR_LOG" ] || [ "$TENSOR_EXPR_LOG" = ":" ]; then
            TENSOR_EXPR_LOG=$(grep -A 5 "Loaded $tensor expression" "$FULL_LOG" | grep -v "Loaded $tensor expression" | grep -v "^$" | grep -v "\[CONSTRUCTOR\]" | grep -E '[+\-*/]|[0-9]' | head -1 | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//' | sed 's/^["'\''"]//' | sed 's/["'\''"]$//')
        fi
        
        if [ -n "$TENSOR_EXPR_LOG" ] && [ "$TENSOR_EXPR_LOG" != ":" ]; then
            echo -e "  ${BLUE}$tensor: $TENSOR_EXPR_LOG${NC}"
        fi
    fi
done

echo ""

if [ -s "symbolic_regression_messages.txt" ]; then
    echo -e "${GREEN}✓ Found symbolic regression messages${NC}"
    echo ""
    
    # Check for successful compilation
    if grep -q "expression compiled successfully" "symbolic_regression_messages.txt"; then
        echo -e "${GREEN}✓ Expressions compiled successfully${NC}"
        grep "expression compiled successfully" "symbolic_regression_messages.txt" | sed 's/^/  /'
    else
        echo -e "${RED}✗ No successful compilation messages found${NC}"
    fi
    
    echo ""
    
    # Check for expression evaluation
    if grep -q "\[DEBUG\].*Using symbolic regression\|Separation correction check" "symbolic_regression_messages.txt"; then
        echo -e "${GREEN}✓ Found symbolic regression evaluation messages${NC}"
        echo ""
        echo "Evaluation status:"
        if grep -q "\[DEBUG\].*Using symbolic regression" "symbolic_regression_messages.txt"; then
            echo -e "  ${GREEN}✓ Expressions are being evaluated${NC}"
            echo ""
            echo "Expression evaluation details:"
            grep "\[DEBUG\].*Using symbolic regression" "symbolic_regression_messages.txt" | sed 's/^/  /'
        fi
        
        # Check for separation correction check messages
        if grep -q "Separation correction check" "symbolic_regression_messages.txt"; then
            echo ""
            echo "Separation correction check results:"
            grep -A 5 "Separation correction check" "symbolic_regression_messages.txt" | head -10 | sed 's/^/  /'
        fi
    else
        echo -e "${YELLOW}⚠ No expression evaluation messages found${NC}"
        echo "  This might mean:"
        echo "    - Expressions are compiled but not being evaluated"
        echo "    - The condition check is failing"
        echo "    - Simulation didn't reach the evaluation point"
        echo "    - Debug output is not being printed"
        echo ""
        echo "  Checking for separationCorrection flag..."
        if grep -qi "separationCorrection" "$FULL_LOG"; then
            echo "  Found separationCorrection references in log"
        else
            echo -e "  ${YELLOW}No separationCorrection found - check if it's enabled in turbulenceProperties${NC}"
        fi
    fi
    
    echo ""
    
    # Check for alpha_S and alpha_A values
    if grep -q "alpha_S" "symbolic_regression_messages.txt"; then
        echo -e "${GREEN}✓ Found alpha_S statistics${NC}"
        echo ""
        echo "alpha_S values:"
        grep "alpha_S" "symbolic_regression_messages.txt" | head -10 | sed 's/^/  /'
    fi
    
    if grep -q "alpha_A" "symbolic_regression_messages.txt"; then
        echo -e "${GREEN}✓ Found alpha_A statistics${NC}"
        echo ""
        echo "alpha_A values:"
        grep "alpha_A" "symbolic_regression_messages.txt" | head -10 | sed 's/^/  /'
    fi
    
    echo ""
    echo "All symbolic regression messages saved to: symbolic_regression_messages.txt"
    echo ""
    echo "First 100 lines of messages:"
    head -100 "symbolic_regression_messages.txt" | sed 's/^/  /'
else
    echo -e "${YELLOW}⚠ No symbolic regression messages found in log${NC}"
    echo "  This might mean:"
    echo "    - Symbolic regression is not enabled"
    echo "    - Expressions are not being loaded"
    echo "    - Debug output is disabled"
    echo ""
    echo "Checking log file for any relevant messages..."
    grep -i "symbolic\|expression\|separation\|anisotropy" "$FULL_LOG" | head -20 | sed 's/^/  /' || echo "  No relevant messages found"
fi

echo ""
echo "================================================================================"
echo "Checking for errors and warnings..."
echo "================================================================================"

# Check for errors
ERROR_COUNT=$(grep -ci "error\|failed\|fatal" "$FULL_LOG" || echo "0")
if [ "$ERROR_COUNT" -gt 0 ]; then
    echo -e "${RED}✗ Found $ERROR_COUNT error(s) in log${NC}"
    echo ""
    echo "First 20 errors:"
    grep -i "error\|failed\|fatal" "$FULL_LOG" | head -20 | sed 's/^/  /'
else
    echo -e "${GREEN}✓ No errors found in log${NC}"
fi

echo ""

# Check for warnings related to symbolic regression
SR_WARNINGS=$(grep -ci "warning.*symbolic\|warning.*expression\|warning.*separation\|warning.*anisotropy" "$FULL_LOG" 2>/dev/null | head -1 || echo "0")
if [ -n "$SR_WARNINGS" ] && [ "$SR_WARNINGS" -gt 0 ] 2>/dev/null; then
    echo -e "${YELLOW}⚠ Found $SR_WARNINGS symbolic regression related warning(s)${NC}"
    echo ""
    echo "Warnings:"
    grep -i "warning.*symbolic\|warning.*expression\|warning.*separation\|warning.*anisotropy" "$FULL_LOG" | head -20 | sed 's/^/  /'
else
    echo -e "${GREEN}✓ No symbolic regression warnings found${NC}"
fi

echo ""
echo "================================================================================"
echo "Summary"
echo "================================================================================"
echo "Full log: $FULL_LOG"
echo "Messages: symbolic_regression_messages.txt"
echo ""

# Final status check
echo "Final Status:"
echo "============="

# Check if expressions were loaded
if grep -q "Loaded.*expression" "$FULL_LOG"; then
    echo -e "${GREEN}✓ Expressions were loaded${NC}"
    LOADED_COUNT=$(grep -c "Loaded.*expression" "$FULL_LOG" || echo "0")
    echo "  Found $LOADED_COUNT expression(s)"
else
    echo -e "${RED}✗ No expressions were loaded${NC}"
fi

# Check if expressions were compiled
if grep -q "expression compiled successfully" "$FULL_LOG"; then
    echo -e "${GREEN}✓ Expressions were compiled${NC}"
    COMPILED_COUNT=$(grep -c "expression compiled successfully" "$FULL_LOG" || echo "0")
    echo "  Found $COMPILED_COUNT compiled expression(s)"
else
    echo -e "${RED}✗ No expressions were compiled successfully${NC}"
fi

# Check if expressions were evaluated
if grep -q "\[DEBUG\].*Using symbolic regression" "$FULL_LOG"; then
    echo -e "${GREEN}✓ Expressions are being evaluated${NC}"
    EVAL_COUNT=$(grep -c "\[DEBUG\].*Using symbolic regression" "$FULL_LOG" || echo "0")
    echo "  Found $EVAL_COUNT evaluation message(s)"
else
    echo -e "${YELLOW}⚠ Expressions may not be evaluated${NC}"
    echo "  Check the log for 'Separation correction check' messages"
fi

# Check simulation status
if [ "$SOLVER_EXIT_CODE" -eq 0 ]; then
    echo -e "${GREEN}✓ Simulation completed successfully${NC}"
else
    echo -e "${RED}✗ Simulation exited with code: $SOLVER_EXIT_CODE${NC}"
fi

echo ""
echo "Useful commands:"
echo "================"
echo "View all symbolic regression messages:"
echo "  cat symbolic_regression_messages.txt"
echo ""
echo "Check expression loading at startup:"
echo "  grep '\[CONSTRUCTOR\].*Loaded' $FULL_LOG"
echo ""
echo "Check expression compilation:"
echo "  grep 'expression compiled' $FULL_LOG"
echo ""
echo "Check expression evaluation:"
echo "  grep '\[DEBUG\].*Using symbolic' $FULL_LOG"
echo ""
echo "Check separation correction status:"
echo "  grep 'Separation correction check' $FULL_LOG"
echo ""
echo "View full log:"
echo "  less $FULL_LOG"
echo ""
