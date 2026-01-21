#!/usr/bin/env python3
"""
check_symbolic_regression.py

Description:    Check and verify symbolic regression configuration
Author:         Mario Javier Rincon, PhD
Affiliation:    Aarhus University/Kamstrup A/S
Contact:        mjrp@mpe.au.dk/mjrp@kamstrup.com
Version:        1.0.0
Last Updated:   2024-01-01

Check symbolic regression configuration and expression dictionaries.
"""

import os
import sys
from pathlib import Path


def check_turbulence_properties(case_dir):
    """Check turbulenceProperties configuration"""
    print("=" * 80)
    print("Checking turbulenceProperties")
    print("=" * 80)
    
    turb_props = Path(case_dir) / "constant" / "turbulenceProperties"
    
    if not turb_props.exists():
        print("[ERROR] constant/turbulenceProperties not found!")
        return False
    
    print(f"[OK] Found: {turb_props}")
    
    with open(turb_props, 'r') as f:
        content = f.read()
    
    # Check if symbolic regression is enabled
    if "useSymbolicRegression" in content:
        if "useSymbolicRegression" in content and "true" in content.split("useSymbolicRegression")[1].split(";")[0]:
            print("  [OK] Symbolic regression is enabled")
            sr_enabled = True
        else:
            print("  [WARNING] Symbolic regression is set but not enabled (should be 'true')")
            sr_enabled = False
    else:
        print("  [WARNING] Symbolic regression not found in turbulenceProperties")
        sr_enabled = False
    
    # Check if debug is enabled
    if "debugSymbolicRegression" in content:
        if "true" in content.split("debugSymbolicRegression")[1].split(";")[0]:
            print("  [OK] Debug output is enabled")
        else:
            print("  [WARNING] Debug output is not enabled (should be 'true')")
    else:
        print("  [WARNING] Debug output not found in turbulenceProperties")
    
    # Check dictionary paths
    if "separationExpressionDict" in content:
        path = content.split("separationExpressionDict")[1].split('"')[1]
        print(f"  [OK] Separation expression dict path: {path}")
    else:
        print("  [WARNING] separationExpressionDict not found")
    
    if "anisotropyExpressionsDict" in content:
        path = content.split("anisotropyExpressionsDict")[1].split('"')[1]
        print(f"  [OK] Anisotropy expressions dict path: {path}")
    else:
        print("  [WARNING] anisotropyExpressionsDict not found")
    
    return sr_enabled


def check_separation_expression(case_dir):
    """Check separation expression dictionary"""
    print("\n" + "=" * 80)
    print("Checking separationExpression.dict")
    print("=" * 80)
    
    sep_dict = Path(case_dir) / "constant" / "separationExpression.dict"
    
    if not sep_dict.exists():
        print("[ERROR] constant/separationExpression.dict not found!")
        return False
    
    print(f"[OK] Found: {sep_dict}")
    
    with open(sep_dict, 'r') as f:
        content = f.read()
    
    # Check for expression key
    if "expression" in content:
        # Extract expression
        expr_start = content.find("expression")
        expr_line = content[expr_start:expr_start+200].split("\n")[0]
        expr = expr_line.split('"')[1] if '"' in expr_line else "not found"
        print(f"  [OK] Expression found: {expr[:80]}...")
    else:
        print("  [ERROR] No 'expression' key found!")
        return False
    
    # Check for variables section
    if "variables" in content:
        print("  [OK] Variables section found")
        # Count normalization constants
        var_count = content.count("_mean_") + content.count("_std_")
        print(f"    Found {var_count} normalization constants")
    else:
        print("  [WARNING] No 'variables' section found")
    
    return True


def check_anisotropy_expressions(case_dir):
    """Check anisotropy expressions dictionary"""
    print("\n" + "=" * 80)
    print("Checking anisotropyExpressions.dict")
    print("=" * 80)
    
    aniso_dict = Path(case_dir) / "constant" / "anisotropyExpressions.dict"
    
    if not aniso_dict.exists():
        print("[ERROR] constant/anisotropyExpressions.dict not found!")
        return False
    
    print(f"[OK] Found: {aniso_dict}")
    
    with open(aniso_dict, 'r') as f:
        content = f.read()
    
    # Check for tensors section
    if "tensors" not in content:
        print("  [ERROR] No 'tensors' section found!")
        return False
    
    print("  [OK] Tensors section found")
    
    # Check for variables section
    if "variables" in content:
        print("  [OK] Variables section found")
    else:
        print("  [WARNING] No 'variables' section found")
    
    # Check for Tij2 expression
    if "Tij2" in content:
        print("  [OK] Tij2 found")
        # Extract Tij2 expression
        tij2_start = content.find("Tij2")
        tij2_section = content[tij2_start:tij2_start+300]
        if "expression" in tij2_section:
            expr_line = [l for l in tij2_section.split("\n") if "expression" in l][0]
            expr = expr_line.split('"')[1] if '"' in expr_line else "not found"
            print(f"    Expression: {expr[:80]}...")
        else:
            print("    [WARNING] Tij2 expression not found")
    else:
        print("  [WARNING] Tij2 not found")
    
    # Count other tensors
    for i in range(3, 11):
        if f"Tij{i}" in content:
            print(f"  [OK] Tij{i} found")
    
    return True


def check_log_files(case_dir):
    """Check log files for symbolic regression messages"""
    print("\n" + "=" * 80)
    print("Checking log files")
    print("=" * 80)
    
    case_path = Path(case_dir)
    log_files = list(case_path.glob("log.*"))
    
    if not log_files:
        print("  [WARNING] No log files found")
        return
    
    # Find the most recent log file
    latest_log = max(log_files, key=lambda p: p.stat().st_mtime)
    print(f"  Checking: {latest_log.name}")
    
    with open(latest_log, 'r', errors='ignore') as f:
        content = f.read()
    
    # Check for symbolic regression messages
    sr_messages = []
    for line in content.split("\n"):
        if any(keyword in line for keyword in ["Symbolic regression", "[Debug]", "Loaded", "expression compiled"]):
            sr_messages.append(line.strip())
    
    if sr_messages:
        print(f"  [OK] Found {len(sr_messages)} symbolic regression messages")
        print("\n  First 10 messages:")
        for msg in sr_messages[:10]:
            print(f"    {msg[:100]}...")
    else:
        print("  [WARNING] No symbolic regression messages found in log")
        print("    This might mean:")
        print("      - Simulation hasn't started yet")
        print("      - Symbolic regression is not enabled")
        print("      - Expressions are not being loaded")


def main():
    """Main function"""
    # Get case directory
    if len(sys.argv) > 1:
        case_dir = sys.argv[1]
    else:
        case_dir = Path(__file__).parent
    
    case_dir = Path(case_dir).resolve()
    
    print("\n" + "=" * 80)
    print("kOmegaSSTPDA Symbolic Regression Configuration Checker")
    print("=" * 80)
    print(f"\nCase directory: {case_dir}\n")
    
    # Run checks
    sr_enabled = check_turbulence_properties(case_dir)
    sep_ok = check_separation_expression(case_dir)
    aniso_ok = check_anisotropy_expressions(case_dir)
    check_log_files(case_dir)
    
    # Summary
    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    
    if sr_enabled and sep_ok and aniso_ok:
        print("[OK] All checks passed! Symbolic regression should be working.")
        print("\nTo verify:")
        print("  1. Run the simulation")
        print("  2. Check log file for 'Symbolic regression: Loaded ... expression' messages")
        print("  3. Wait for 100+ iterations to see debug output")
    else:
        print("[WARNING] Some issues found. Please fix the above errors before running.")
    
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
