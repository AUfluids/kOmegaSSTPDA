#!/bin/bash
#
# ===============================================================================
# FILE: convert_pysr.sh
# ===============================================================================
#
# DESCRIPTION:
#     Convert PySR model outputs to OpenFOAM dictionary format
#     Helper script for SD_test case
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

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Function to find pysr_to_openfoam.py script
find_pysr_converter() {
    local found_script=""
    
    # Strategy 1: Check environment variable (for HPC/custom setups)
    if [ -n "$KOMEGASSTPDA_ROOT" ]; then
        local candidate="$KOMEGASSTPDA_ROOT/kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py"
        if [ -f "$candidate" ]; then
            found_script="$candidate"
            return 0
        fi
    fi
    
    # Strategy 2: Try relative path from script location (standard repo structure)
    local candidate="$SCRIPT_DIR/../../kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py"
    if [ -f "$candidate" ]; then
        found_script="$candidate"
        return 0
    fi
    
    # Strategy 3: Search upward from script directory to find repo root
    local current_dir="$SCRIPT_DIR"
    while [ "$current_dir" != "/" ]; do
        # Check if we're in the kOmegaSSTPDA repo root
        if [ -f "$current_dir/kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py" ]; then
            found_script="$current_dir/kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py"
            return 0
        fi
        # Check if parent is the repo root
        local parent_dir="$(dirname "$current_dir")"
        if [ -f "$parent_dir/kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py" ]; then
            found_script="$parent_dir/kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py"
            return 0
        fi
        current_dir="$parent_dir"
    done
    
    # Strategy 4: Check if script exists in current directory (for copied cases)
    if [ -f "$SCRIPT_DIR/pysr_to_openfoam.py" ]; then
        found_script="$SCRIPT_DIR/pysr_to_openfoam.py"
        return 0
    fi
    
    return 1
}

# Find the converter script
if ! PYSR_CONVERTER=$(find_pysr_converter); then
    echo "Error: PySR converter script (pysr_to_openfoam.py) not found!"
    echo ""
    echo "Searched locations:"
    echo "  1. KOMEGASSTPDA_ROOT/kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py"
    echo "  2. ../../kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py (relative to script)"
    echo "  3. Upward search from script directory"
    echo "  4. ./pysr_to_openfoam.py (in case directory)"
    echo ""
    echo "Solutions:"
    echo "  - Set KOMEGASSTPDA_ROOT environment variable to repository root"
    echo "    export KOMEGASSTPDA_ROOT=/path/to/kOmegaSSTPDA"
    echo ""
    echo "  - Copy pysr_to_openfoam.py to this case directory"
    echo "    cp \$KOMEGASSTPDA_ROOT/kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py ."
    echo ""
    echo "  - Use Python directly with full path:"
    echo "    python /path/to/kOmegaSSTPDA/kOmegaSSTPDA_symbolic_regression/pysr_to_openfoam.py ..."
    exit 1
fi

# Verify the found script is executable/readable
if [ ! -r "$PYSR_CONVERTER" ]; then
    echo "Error: Found script at $PYSR_CONVERTER but cannot read it"
    exit 1
fi

# Default values
EXPRESSION_TYPE="separation"
TENSOR_INDEX=2

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--type)
            EXPRESSION_TYPE="$2"
            shift 2
            ;;
        -i|--tensor-index)
            TENSOR_INDEX="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS] INPUT_FILE [OUTPUT_FILE]"
            echo ""
            echo "Convert PySR JSON output to OpenFOAM dictionary format"
            echo ""
            echo "Options:"
            echo "  -t, --type TYPE          Expression type: 'separation' or 'anisotropy' (default: separation)"
            echo "  -i, --tensor-index N    Tensor index for anisotropy expressions (default: 2)"
            echo "  -h, --help              Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0 pysr_separation.json"
            echo "  $0 pysr_anisotropy.json -t anisotropy -i 2"
            exit 0
            ;;
        *)
            if [ -z "$INPUT_FILE" ]; then
                INPUT_FILE="$1"
            elif [ -z "$OUTPUT_FILE" ]; then
                OUTPUT_FILE="$1"
            else
                echo "Error: Too many arguments"
                exit 1
            fi
            shift
            ;;
    esac
done

# Check input file
if [ -z "$INPUT_FILE" ]; then
    echo "Error: Input file not specified"
    echo "Use -h or --help for usage information"
    exit 1
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
    exit 1
fi

# Determine output file
if [ -z "$OUTPUT_FILE" ]; then
    if [ "$EXPRESSION_TYPE" == "separation" ]; then
        OUTPUT_FILE="constant/separationExpression.dict"
    else
        OUTPUT_FILE="constant/anisotropyExpression_Tij${TENSOR_INDEX}.dict"
    fi
fi

# Create constant directory if it doesn't exist
mkdir -p "$SCRIPT_DIR/constant"

# Convert
echo "Converting PySR model..."
echo "  Input:  $INPUT_FILE"
echo "  Output: $OUTPUT_FILE"
echo "  Type:   $EXPRESSION_TYPE"
if [ "$EXPRESSION_TYPE" == "anisotropy" ]; then
    echo "  Tensor: Tij${TENSOR_INDEX}"
fi

# Show which script is being used (helpful for debugging)
if [ -n "$KOMEGASSTPDA_ROOT" ]; then
    echo "  Using converter from: $PYSR_CONVERTER (via KOMEGASSTPDA_ROOT)"
else
    echo "  Using converter from: $PYSR_CONVERTER"
fi
echo ""

python3 "$PYSR_CONVERTER" \
    "$INPUT_FILE" \
    -o "$SCRIPT_DIR/$OUTPUT_FILE" \
    -t "$EXPRESSION_TYPE" \
    $(if [ "$EXPRESSION_TYPE" == "anisotropy" ]; then echo "--tensor-index $TENSOR_INDEX"; fi)

if [ $? -eq 0 ]; then
    echo ""
    echo "Success! Expression dictionary created: $OUTPUT_FILE"
    echo ""
    echo "Next steps:"
    echo "  1. Review the expression in: $OUTPUT_FILE"
    echo "  2. Update constant/turbulenceProperties to enable symbolic regression"
    echo "  3. Run simulation: ./Allrun.sh"
else
    echo ""
    echo "Error: Conversion failed"
    exit 1
fi
