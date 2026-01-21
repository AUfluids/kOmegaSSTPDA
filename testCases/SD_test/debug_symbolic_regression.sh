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

# Check expression dictionaries
echo ""
echo "Checking expression dictionaries..."
if [ -f "constant/separationExpression.dict" ]; then
    echo -e "${GREEN}✓ Found: constant/separationExpression.dict${NC}"
    if grep -q "expression" constant/separationExpression.dict; then
        echo "  Expression found in dictionary"
    else
        echo -e "${RED}  ERROR: No 'expression' key found!${NC}"
    fi
else
    echo -e "${RED}✗ Missing: constant/separationExpression.dict${NC}"
fi

if [ -f "constant/anisotropyExpressions.dict" ]; then
    echo -e "${GREEN}✓ Found: constant/anisotropyExpressions.dict${NC}"
    if grep -q "tensors" constant/anisotropyExpressions.dict; then
        echo "  Tensors section found"
        if grep -A 5 "Tij2" constant/anisotropyExpressions.dict | grep -q "expression"; then
            echo "  Tij2 expression found"
        else
            echo -e "${YELLOW}  Warning: Tij2 expression not found${NC}"
        fi
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

# Function to filter symbolic regression output
filter_sr_output() {
    grep -E "(Symbolic regression|\[Debug\].*Symbolic|Loaded.*expression|expression compiled|Failed to compile|useSymbolicRegression|debugSymbolicRegression|Attempting to read)" || true
}

# Run solver and capture output
echo "Starting simulation..."
echo ""

# Run solver and tee output to both full log and filtered log
$SOLVER 2>&1 | tee "$FULL_LOG" | \
    while IFS= read -r line; do
        # Print all lines, but highlight symbolic regression messages
        if echo "$line" | grep -qE "(Symbolic regression|\[Debug\].*Symbolic|Loaded.*expression|expression compiled|Failed to compile)"; then
            echo -e "${GREEN}$line${NC}"
        elif echo "$line" | grep -qE "(Warning|Error|Failed)"; then
            echo -e "${RED}$line${NC}"
        else
            echo "$line"
        fi
    done | tee "$LOG_FILE"

# Extract symbolic regression messages to a separate file
echo ""
echo "================================================================================"
echo "Extracting symbolic regression messages..."
echo "================================================================================"
grep -E "(Symbolic regression|\[Debug\].*Symbolic|Loaded.*expression|expression compiled|Failed to compile|useSymbolicRegression|debugSymbolicRegression|Attempting to read)" "$FULL_LOG" > "symbolic_regression_messages.txt" || true

if [ -s "symbolic_regression_messages.txt" ]; then
    echo -e "${GREEN}✓ Found symbolic regression messages${NC}"
    echo ""
    echo "Messages saved to: symbolic_regression_messages.txt"
    echo ""
    echo "First 50 lines:"
    head -50 "symbolic_regression_messages.txt"
else
    echo -e "${YELLOW}⚠ No symbolic regression messages found in log${NC}"
    echo "  This might mean:"
    echo "    - Symbolic regression is not enabled"
    echo "    - Expressions are not being loaded"
    echo "    - Simulation didn't reach 100 iterations (debug output appears every 100 iterations)"
fi

echo ""
echo "================================================================================"
echo "Summary"
echo "================================================================================"
echo "Full log: $FULL_LOG"
echo "Filtered log: $LOG_FILE"
echo "Messages: symbolic_regression_messages.txt"
echo ""
echo "To view all symbolic regression messages:"
echo "  cat symbolic_regression_messages.txt"
echo ""
echo "To check if expressions were loaded at startup:"
echo "  grep 'Loaded.*expression' $FULL_LOG"
echo ""
echo "To check debug output (every 100 iterations):"
echo "  grep '\[Debug\].*Symbolic' $FULL_LOG"
echo ""
