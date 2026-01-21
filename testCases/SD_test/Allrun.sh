#!/bin/bash
#
# ===============================================================================
# FILE: Allrun.sh
# ===============================================================================
#
# DESCRIPTION:
#     Main run script for OpenFOAM case with symbolic regression support
#     Incorporates symbolic regression checks and enhanced logging
#
# AUTHOR INFORMATION:
#     Name:           Mario Javier Rincon, PhD
#     Title:          Research Scientist
#     Affiliation:    Aarhus University/Kamstrup A/S
#     Email:          mjrp@mpe.au.dk/mjrp@kamstrup.com
#
# METADATA:
#     Created:        2020-12-10
#     Last Modified:  2024-01-01
#     Version:        2.0.0
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
echo "OpenFOAM Case Runner"
echo "================================================================================"
echo ""

# Create logs directory
if [ -d "logs" ]; then
    echo -e "${YELLOW}⚠ Removing existing logs directory...${NC}"
    rm -rf logs
fi
mkdir -p logs

# Clean previous results
echo -e "${BLUE}Cleaning previous results...${NC}"
foamListTimes -rm 2>/dev/null || true

# Copy initial conditions
if [ -d "0.orig" ]; then
    echo -e "${BLUE}Copying initial conditions from 0.orig...${NC}"
    cp -r 0.orig 0
else
    echo -e "${YELLOW}⚠ Warning: 0.orig directory not found${NC}"
fi

# Check if turbulenceProperties exists
if [ ! -f "constant/turbulenceProperties" ]; then
    echo -e "${RED}ERROR: constant/turbulenceProperties not found!${NC}"
    exit 1
fi

# Check for symbolic regression configuration
SYMBOLIC_REGRESSION_ENABLED=false
if grep -q "useSymbolicRegression.*true" constant/turbulenceProperties; then
    SYMBOLIC_REGRESSION_ENABLED=true
    echo -e "${GREEN}✓ Symbolic regression is enabled${NC}"
    
    # Check expression dictionaries
    echo ""
    echo "Checking symbolic regression expression dictionaries..."
    
    if [ -f "constant/separationExpression.dict" ]; then
        echo -e "${GREEN}✓ Found: constant/separationExpression.dict${NC}"
        if grep -q "expression" constant/separationExpression.dict; then
            echo "  Separation expression dictionary is valid"
        else
            echo -e "${YELLOW}  ⚠ Warning: No 'expression' key found in separationExpression.dict${NC}"
        fi
    else
        echo -e "${YELLOW}⚠ Warning: constant/separationExpression.dict not found${NC}"
    fi
    
    if [ -f "constant/anisotropyExpressions.dict" ]; then
        echo -e "${GREEN}✓ Found: constant/anisotropyExpressions.dict${NC}"
        if grep -q "tensors" constant/anisotropyExpressions.dict; then
            echo "  Anisotropy expressions dictionary is valid"
            # Count active tensors
            TENSOR_COUNT=$(grep -c "Tij[0-9]" constant/anisotropyExpressions.dict | grep -v "^[[:space:]]*//" || echo "0")
            if [ "$TENSOR_COUNT" -gt 0 ]; then
                echo "  Found expressions for tensor(s)"
            fi
        else
            echo -e "${YELLOW}  ⚠ Warning: No 'tensors' section found${NC}"
        fi
    else
        echo -e "${YELLOW}⚠ Warning: constant/anisotropyExpressions.dict not found${NC}"
    fi
    
    # Check debug flag
    if grep -q "debugSymbolicRegression.*true" constant/turbulenceProperties; then
        echo -e "${GREEN}✓ Debug output is enabled${NC}"
    else
        echo -e "${BLUE}ℹ Debug output is disabled (set debugSymbolicRegression true to enable)${NC}"
    fi
else
    echo -e "${BLUE}ℹ Symbolic regression is not enabled${NC}"
fi

# Determine solver
SOLVER="simpleFoam"
if [ -f "system/controlDict" ]; then
    SOLVER=$(grep "^application" system/controlDict 2>/dev/null | sed 's/.*application[[:space:]]*\([^;]*\);/\1/' | tr -d ' ' || echo "simpleFoam")
    if [ -z "$SOLVER" ]; then
        SOLVER="simpleFoam"
    fi
fi

echo ""
echo "================================================================================"
echo "Running simulation: $SOLVER"
echo "================================================================================"
echo ""

# Create timestamped log file
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOG_FILE="logs/simulation.${TIMESTAMP}"
FULL_LOG="logs/simulation"

echo "Log file: $LOG_FILE"
echo "Starting simulation..."
echo ""

# Create start marker
touch STARTFILE

# Run solver and redirect output to log file
echo -e "${BLUE}Running $SOLVER...${NC}"
$SOLVER > "$LOG_FILE" 2>&1
SOLVER_EXIT_CODE=$?

# Also create a symlink or copy to the standard name for compatibility
cp "$LOG_FILE" "$FULL_LOG" 2>/dev/null || true

echo ""
if [ "$SOLVER_EXIT_CODE" -eq 0 ]; then
    echo -e "${GREEN}✓ Simulation completed successfully${NC}"
else
    echo -e "${RED}✗ Simulation exited with code: $SOLVER_EXIT_CODE${NC}"
fi

# Analyze symbolic regression output if enabled
if [ "$SYMBOLIC_REGRESSION_ENABLED" = true ]; then
    echo ""
    echo "================================================================================"
    echo "Symbolic Regression Analysis"
    echo "================================================================================"
    
    # Check for expression loading
    if grep -q "Loaded.*expression" "$LOG_FILE"; then
        echo -e "${GREEN}✓ Expressions were loaded${NC}"
        LOADED_COUNT=$(grep -c "Loaded.*expression" "$LOG_FILE" || echo "0")
        echo "  Found $LOADED_COUNT expression(s)"
    else
        echo -e "${YELLOW}⚠ No expressions were loaded${NC}"
    fi
    
    # Check for successful compilation
    if grep -q "expression compiled successfully" "$LOG_FILE"; then
        echo -e "${GREEN}✓ Expressions compiled successfully${NC}"
        COMPILED_COUNT=$(grep -c "expression compiled successfully" "$LOG_FILE" || echo "0")
        echo "  Found $COMPILED_COUNT compiled expression(s)"
    else
        echo -e "${YELLOW}⚠ No successful compilation messages found${NC}"
    fi
    
    # Check for expression evaluation
    if grep -q "\[DEBUG\].*Using symbolic regression" "$LOG_FILE"; then
        echo -e "${GREEN}✓ Expressions are being evaluated${NC}"
        EVAL_COUNT=$(grep -c "\[DEBUG\].*Using symbolic regression" "$LOG_FILE" || echo "0")
        echo "  Found $EVAL_COUNT evaluation message(s)"
    else
        echo -e "${BLUE}ℹ Expression evaluation messages not found (may be normal if debug is disabled)${NC}"
    fi
    
    # Extract symbolic regression messages to separate file
    grep -E "(Symbolic regression|\[CONSTRUCTOR\].*SYMBOLIC|\[CONSTRUCTOR\].*Loaded|\[CONSTRUCTOR\].*expression|\[DEBUG\].*symbolic|\[DEBUG\].*alpha_S|\[DEBUG\].*alpha_A|expression compiled|Failed to compile)" "$LOG_FILE" > "logs/symbolic_regression_messages.txt" 2>/dev/null || true
    
    if [ -s "logs/symbolic_regression_messages.txt" ]; then
        echo ""
        echo "Symbolic regression messages saved to: logs/symbolic_regression_messages.txt"
    fi
fi

# Check for errors
echo ""
echo "================================================================================"
echo "Error Check"
echo "================================================================================"
ERROR_COUNT=$(grep -ci "error\|failed\|fatal" "$LOG_FILE" 2>/dev/null || echo "0")
if [ "$ERROR_COUNT" -gt 0 ]; then
    echo -e "${RED}✗ Found $ERROR_COUNT error(s) in log${NC}"
    echo ""
    echo "First 10 errors:"
    grep -i "error\|failed\|fatal" "$LOG_FILE" | head -10 | sed 's/^/  /'
else
    echo -e "${GREEN}✓ No errors found in log${NC}"
fi

# Post-processing
echo ""
echo "================================================================================"
echo "Post-processing"
echo "================================================================================"

# Run foamLog if available
if command -v foamLog &> /dev/null; then
    echo -e "${BLUE}Running foamLog...${NC}"
    foamLog "$FULL_LOG" 2>/dev/null || echo -e "${YELLOW}  ⚠ foamLog completed with warnings${NC}"
else
    echo -e "${YELLOW}⚠ foamLog not available, skipping${NC}"
fi

# Run Python post-processing scripts if they exist
if [ -f "post_process.py" ]; then
    echo -e "${BLUE}Running post_process.py...${NC}"
    if command -v python3 &> /dev/null; then
        python3 post_process.py 2>&1 | tee logs/post_process.log || echo -e "${YELLOW}  ⚠ post_process.py completed with warnings${NC}"
    elif command -v python &> /dev/null; then
        python post_process.py 2>&1 | tee logs/post_process.log || echo -e "${YELLOW}  ⚠ post_process.py completed with warnings${NC}"
    else
        echo -e "${YELLOW}  ⚠ Python not found, skipping post_process.py${NC}"
    fi
else
    echo -e "${BLUE}ℹ post_process.py not found, skipping${NC}"
fi

if [ -f "plot_residuals.py" ]; then
    echo -e "${BLUE}Running plot_residuals.py...${NC}"
    if command -v python3 &> /dev/null; then
        python3 plot_residuals.py 2>&1 | tee logs/plot_residuals.log || echo -e "${YELLOW}  ⚠ plot_residuals.py completed with warnings${NC}"
    elif command -v python &> /dev/null; then
        python plot_residuals.py 2>&1 | tee logs/plot_residuals.log || echo -e "${YELLOW}  ⚠ plot_residuals.py completed with warnings${NC}"
    else
        echo -e "${YELLOW}  ⚠ Python not found, skipping plot_residuals.py${NC}"
    fi
else
    echo -e "${BLUE}ℹ plot_residuals.py not found, skipping${NC}"
fi

# Summary
echo ""
echo "================================================================================"
echo "Summary"
echo "================================================================================"
echo "Simulation log: $LOG_FILE"
echo "Standard log:   $FULL_LOG"

if [ "$SYMBOLIC_REGRESSION_ENABLED" = true ]; then
    echo "SR messages:     logs/symbolic_regression_messages.txt"
fi

if [ "$SOLVER_EXIT_CODE" -eq 0 ]; then
    echo ""
    echo -e "${GREEN}✓ Simulation completed successfully${NC}"
else
    echo ""
    echo -e "${RED}✗ Simulation exited with code: $SOLVER_EXIT_CODE${NC}"
    echo "  Check log file for details: $LOG_FILE"
fi

# Create done marker
touch DONEFILE

echo ""
echo "================================================================================"
echo "Useful Commands"
echo "================================================================================"
echo "View simulation log:"
echo "  less $LOG_FILE"
echo ""
if [ "$SYMBOLIC_REGRESSION_ENABLED" = true ]; then
    echo "View symbolic regression messages:"
    echo "  less logs/symbolic_regression_messages.txt"
    echo ""
    echo "Check expression loading:"
    echo "  grep '\[CONSTRUCTOR\].*Loaded' $LOG_FILE"
    echo ""
    echo "Check expression evaluation:"
    echo "  grep '\[DEBUG\].*Using symbolic' $LOG_FILE"
    echo ""
fi
echo "View residuals:"
echo "  foamLog $FULL_LOG"
echo ""
echo "================================================================================"
echo ""
echo -e "${GREEN}All done! :)${NC}"
echo ""
