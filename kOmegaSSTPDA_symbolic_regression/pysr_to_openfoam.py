#!/usr/bin/env python3
"""
pysr_to_openfoam.py

Description:    Convert PySR symbolic regression output to OpenFOAM dictionary format
Author:         Mario Javier Rincon, PhD
Affiliation:    Aarhus University/Kamstrup A/S
Contact:        mjrp@mpe.au.dk/mjrp@kamstrup.com
Version:        1.0.0
Last Updated:   2024-01-01

This script converts PySR model outputs to OpenFOAM dictionary format
for use with the kOmegaSSTPDA turbulence model with symbolic regression.
"""

import json
import sys
import argparse
import re
from pathlib import Path


def pysr_equation_to_openfoam(equation_str, variable_map=None):
    """
    Convert PySR equation string to OpenFOAM-compatible format.
    
    PySR typically outputs expressions with numeric coefficients, e.g.:
    - "-2.070 + 1.119*x1 - 0.215*x2"
    - "2.5*x1*x2 + 0.8*exp(x3)"
    
    Parameters:
    -----------
    equation_str : str
        PySR equation string (e.g., "-2.070 + 1.119*x1 - 0.215*x2" or "C0 + C1*I1")
    variable_map : dict, optional
        Mapping from PySR variable names to OpenFOAM field names
        Common mappings: {'x1': 'I1', 'x2': 'I2', ...} or {'I1': 'normalised_I1', ...}
        
    Returns:
    --------
    str : OpenFOAM-compatible expression string
    """
    if variable_map is None:
        variable_map = {}
    
    # Replace variable names according to mapping
    # Use word boundaries to avoid partial replacements
    result = equation_str
    
    # Sort by length (longest first) to avoid partial replacements
    sorted_map = sorted(variable_map.items(), key=lambda x: len(x[0]), reverse=True)
    
    for pysr_name, of_name in sorted_map:
        # Replace with word boundaries to avoid partial matches
        # Handle cases like x1, x2, x10 (need to be careful with x1 vs x10)
        # Use word boundaries for all variables to avoid partial matches
        pattern = r'\b' + re.escape(pysr_name) + r'\b'
        result = re.sub(pattern, of_name, result)
    
    # PySR uses Python syntax, OpenFOAM uses C++ syntax
    # Most operators are the same, but handle special cases:
    # - Python: ** for power, OpenFOAM: ^ or pow()
    # - Python: abs(), OpenFOAM: abs() or fabs()
    result = result.replace('**', '^')  # Power operator
    
    return result


def create_separation_expression_dict(pysr_model, output_file):
    """
    Create OpenFOAM dictionary for separation correction expression.
    
    Parameters:
    -----------
    pysr_model : dict or str
        PySR model (either JSON dict or path to JSON file)
    output_file : str
        Path to output OpenFOAM dictionary file
    """
    # Load PySR model if it's a file path
    if isinstance(pysr_model, (str, Path)):
        with open(pysr_model, 'r') as f:
            model_data = json.load(f)
    else:
        model_data = pysr_model
    
    # Extract equation (PySR format may vary, adjust as needed)
    # PySR typically stores equations in different formats:
    # - model.get_best() returns dict with 'equation' key
    # - JSON export might have 'equation' or 'best' -> 'equation'
    # - Equation string might use x1, x2, ... or I1, I2, ... as variable names
    if 'equation' in model_data:
        equation = model_data['equation']
    elif 'best' in model_data:
        if isinstance(model_data['best'], dict) and 'equation' in model_data['best']:
            equation = model_data['best']['equation']
        else:
            equation = str(model_data['best'])
    elif isinstance(model_data, str):
        equation = model_data
    else:
        raise ValueError("Could not find equation in PySR model data. Expected 'equation' key or 'best' -> 'equation'")
    
    # Variable mapping: PySR variables to OpenFOAM invariant names
    # Note: Normalisation is handled automatically by the turbulence model
    # The model normalises invariants before registering them with the expression parser
    var_map = {
        # Common PySR variable names (x1, x2, ...)
        'x1': 'I1',
        'x2': 'I2',
        'x3': 'I3',
        'x4': 'I4',
        'x5': 'I5',
        # Alternative: if PySR uses I1, I2, ... directly
        'I1': 'I1',
        'I2': 'I2',
        'I3': 'I3',
        'I4': 'I4',
        'I5': 'I5',
    }
    
    # Convert equation
    of_equation = pysr_equation_to_openfoam(equation, var_map)
    
    # Create OpenFOAM dictionary content
    dict_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      separationExpression;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Symbolic regression expression for separation correction
// Generated from PySR model

expression     "{of_equation}";

// Variable definitions (normalisation constants)
// These are used by the turbulence model to normalise invariants: (I - mean) / std
// The normalisation is applied automatically before expression evaluation
// PySR expressions typically have coefficients embedded, so no separate coefficient definitions needed
variables
{{
    // Separation correction normalisation constants
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
}}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write(dict_content)
    
    print(f"Created separation expression dictionary: {output_file}")
    print(f"Expression: {of_equation}")


def create_consolidated_anisotropy_dict(pysr_models_dict, output_file):
    """
    Create consolidated OpenFOAM dictionary for all anisotropy correction expressions.
    
    Parameters:
    -----------
    pysr_models_dict : dict
        Dictionary mapping tensor indices (as int or str like "2", "3") to PySR models
        Example: {2: model_data, 3: model_data, ...}
    output_file : str
        Path to output OpenFOAM dictionary file
    """
    # Variable mapping: PySR variables to OpenFOAM invariant names
    # Note: Normalisation is handled automatically by the turbulence model
    # The model normalises invariants before registering them with the expression parser
    var_map = {
        # Common PySR variable names (x1, x2, ...)
        'x1': 'I1',
        'x2': 'I2',
        'x3': 'I3',
        'x4': 'I4',
        'x5': 'I5',
        # Alternative: if PySR uses I1, I2, ... directly
        'I1': 'I1',
        'I2': 'I2',
        'I3': 'I3',
        'I4': 'I4',
        'I5': 'I5',
    }
    
    # Build tensor expressions section
    tensor_sections = []
    for tensor_idx, model_data in sorted(pysr_models_dict.items()):
        tensor_idx = int(tensor_idx)  # Ensure integer
        
        # Extract equation
        if isinstance(model_data, dict):
            if 'equation' in model_data:
                equation = model_data['equation']
            elif 'best' in model_data:
                if isinstance(model_data['best'], dict) and 'equation' in model_data['best']:
                    equation = model_data['best']['equation']
                else:
                    equation = str(model_data['best'])
            else:
                print(f"Warning: Could not find equation for Tij{tensor_idx}, skipping")
                continue
        else:
            equation = str(model_data)
        
        # Convert equation
        of_equation = pysr_equation_to_openfoam(equation, var_map)
        
        # Create tensor section
        tensor_section = f"""    // Tij{tensor_idx} expression
    // PySR-generated expression with embedded coefficients
    Tij{tensor_idx}
    {{
        expression     "{of_equation}";
        
        // Note: Coefficients are embedded in the expression above
        // If you need to override coefficients, you can add them here:
        // coefficients
        // {{
        //     A0_{tensor_idx}                   0.0;
        //     A1_{tensor_idx}                   0.0;
        //     ...
        // }}
    }}"""
        tensor_sections.append(tensor_section)
    
    if not tensor_sections:
        raise ValueError("No valid tensor expressions found in input models")
    
    # Create OpenFOAM dictionary content
    dict_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      anisotropyExpressions;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Symbolic regression expressions for anisotropy correction
// Generated from PySR models
//
// Each tensor (Tij2, Tij3, etc.) can have its own expression.
// Expressions replace the hardcoded calculations:
//   alpha_A_i = A0_i + A1_i*I1 + A2_i*I2 + A3_i*I3 + A4_i*I4 + A5_i*I5
//
// Variables I1-I5 in expressions are automatically normalised by the turbulence model
// using mean and std values below before expression evaluation
// Expression syntax follows C++/ExprTK format - use I1, I2, etc. directly (not normalised form)

// Variable definitions (normalisation constants) - shared by all tensors
variables
{{
    // Anisotropy correction normalisation constants
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
}}

// Tensor expressions - add/remove tensors as needed
// Only tensors with non-empty expressions will be used
tensors
{{
{chr(10).join(tensor_sections)}
    
    // Add more tensors (Tij4-Tij10) as needed following the same pattern
}}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write(dict_content)
    
    print(f"Created consolidated anisotropy expressions dictionary: {output_file}")
    print(f"Added {len(tensor_sections)} tensor expression(s): {', '.join([f'Tij{int(k)}' for k in sorted(pysr_models_dict.keys())])}")


def create_anisotropy_expression_dict(pysr_models, output_file, tensor_index):
    """
    Create OpenFOAM dictionary for anisotropy correction expression.
    
    Parameters:
    -----------
    pysr_models : dict
        Dictionary mapping tensor indices to PySR models
    output_file : str
        Path to output OpenFOAM dictionary file
    tensor_index : int
        Tensor index (2, 3, 4, 5, etc.)
    """
    if tensor_index not in pysr_models:
        raise ValueError(f"No PySR model found for tensor index {tensor_index}")
    
    model_data = pysr_models[tensor_index]
    
    # Extract equation
    if isinstance(model_data, dict):
        if 'equation' in model_data:
            equation = model_data['equation']
        elif 'best' in model_data:
            if isinstance(model_data['best'], dict) and 'equation' in model_data['best']:
                equation = model_data['best']['equation']
            else:
                equation = str(model_data['best'])
        else:
            raise ValueError("Could not find equation in PySR model data")
    else:
        equation = str(model_data)
    
    # Variable mapping: PySR variables to OpenFOAM invariant names
    # Note: Normalisation is handled automatically by the turbulence model
    # The model normalises invariants before registering them with the expression parser
    var_map = {
        'x1': 'I1',
        'x2': 'I2',
        'x3': 'I3',
        'x4': 'I4',
        'x5': 'I5',
        'I1': 'I1',
        'I2': 'I2',
        'I3': 'I3',
        'I4': 'I4',
        'I5': 'I5',
    }
    
    # Convert equation
    of_equation = pysr_equation_to_openfoam(equation, var_map)
    
    # Create OpenFOAM dictionary content
    dict_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      anisotropyExpression_Tij{tensor_index};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Symbolic regression expression for anisotropy correction (Tij{tensor_index})
// Generated from PySR model

expression     "{of_equation}";

// Variable definitions (will be registered at runtime)
variables
{{
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
}}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write(dict_content)
    
    print(f"Created anisotropy expression dictionary: {output_file}")
    print(f"Expression: {of_equation}")


def main():
    parser = argparse.ArgumentParser(
        description='Convert PySR output to OpenFOAM dictionary format'
    )
    parser.add_argument(
        'input',
        type=str,
        help='Input PySR JSON file or equation string'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='expression.dict',
        help='Output OpenFOAM dictionary file'
    )
    parser.add_argument(
        '-t', '--type',
        type=str,
        choices=['separation', 'anisotropy'],
        default='separation',
        help='Expression type: separation or anisotropy'
    )
    parser.add_argument(
        '--tensor-index',
        type=int,
        default=2,
        help='Tensor index for anisotropy expressions (default: 2). Use --consolidate to create single file with multiple tensors.'
    )
    parser.add_argument(
        '--consolidate',
        action='store_true',
        help='Create consolidated anisotropy expressions file with multiple tensors. Requires --tensor-indices or multiple --tensor-index arguments.'
    )
    parser.add_argument(
        '--tensor-indices',
        type=str,
        help='Comma-separated list of tensor indices for consolidated file (e.g., "2,3,4"). Each index requires a corresponding input file or JSON key.'
    )
    
    args = parser.parse_args()
    
    # Create appropriate dictionary
    if args.type == 'separation':
        # Check if input is a file or a string
        input_path = Path(args.input)
        if input_path.exists():
            with open(input_path, 'r') as f:
                content = f.read()
                try:
                    model_data = json.loads(content)
                except json.JSONDecodeError:
                    # Assume it's a plain equation string
                    model_data = {'equation': content.strip()}
        else:
            # Assume it's an equation string
            model_data = {'equation': args.input}
        
        create_separation_expression_dict(model_data, args.output)
    
    elif args.consolidate:
        # Consolidated anisotropy expressions
        input_path = Path(args.input)
        
        if args.tensor_indices:
            # Multiple tensors from single JSON file
            tensor_indices = [int(x.strip()) for x in args.tensor_indices.split(',')]
            
            if not input_path.exists():
                raise ValueError(f"Input file not found: {args.input}")
            
            with open(input_path, 'r') as f:
                all_models = json.load(f)
            
            # Extract models for each tensor index
            models_dict = {}
            for idx in tensor_indices:
                if str(idx) in all_models:
                    models_dict[idx] = all_models[str(idx)]
                elif idx in all_models:
                    models_dict[idx] = all_models[idx]
                else:
                    print(f"Warning: No model found for tensor index {idx} in input file")
            
            if not models_dict:
                raise ValueError("No valid tensor models found in input file")
            
            create_consolidated_anisotropy_dict(models_dict, args.output)
        else:
            # Single tensor but use consolidated format
            if input_path.exists():
                with open(input_path, 'r') as f:
                    content = f.read()
                    try:
                        model_data = json.loads(content)
                    except json.JSONDecodeError:
                        model_data = {'equation': content.strip()}
            else:
                model_data = {'equation': args.input}
            
            create_consolidated_anisotropy_dict(
                {args.tensor_index: model_data},
                args.output
            )
    
    else:
        # Individual anisotropy expression dictionary
        input_path = Path(args.input)
        if input_path.exists():
            with open(input_path, 'r') as f:
                content = f.read()
                try:
                    model_data = json.loads(content)
                except json.JSONDecodeError:
                    model_data = {'equation': content.strip()}
        else:
            model_data = {'equation': args.input}
        
        create_anisotropy_expression_dict(
            {args.tensor_index: model_data},
            args.output,
            args.tensor_index
        )


if __name__ == '__main__':
    main()
