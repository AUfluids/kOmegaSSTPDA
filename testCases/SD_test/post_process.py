#!/usr/bin/env python3
"""
post_process_SD.py

Description:    Post-processing script for SD test case. Computes error metrics
                and generates visualisation plots using PyVista with a slice at
                x=0 plane, comparing RANS results with LES reference data.

Author:         Mario Javier Rincon, PhD
Affiliation:    Aarhus University/Kamstrup A/S
Contact:        mjrp@mpe.au.dk/mjrp@kamstrup.com
Version:        2.0.0
Last Updated:   2024-01-01
"""

import json
import logging
from pathlib import Path
from typing import Tuple, Optional

import click
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from scipy.interpolate import griddata, LinearNDInterpolator, interp1d
from scipy.spatial.distance import cdist
from scipy.stats import beta
from matplotlib.lines import Line2D
from rich.logging import RichHandler

# Configure logging using RichHandler for coloured output
logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)],
)
logger = logging.getLogger(__name__)

# Constants
UB = 0.4820072072378321  # Bulk velocity
H = 0.5  # Channel half-height

# Matplotlib configuration
plt.rcParams.update({
    "text.usetex": False,
    'font.size': 7,
    'axes.linewidth': 2,
    "grid.color": "#cccccc",
    "axes.grid": True,
    "grid.alpha": 0.8,
    "grid.linewidth": 0.5,
    "grid.linestyle": '-',
    "axes.grid.which": 'both',
    "axes.spines.right": False,
    "axes.spines.top": False,
    'axes.axisbelow': True,
})

cm = 1 / 2.54  # inches to cm
figSizeMedium = (14 * cm, 14 * cm * 3 / 4)

# Base figure size constants (matching SD_test_multi_.png)
# Multi plot uses: 1 row × 3 cols = figSizeMedium
BASE_WIDTH_PER_COL = 14 * cm / 3  # Width per column
BASE_HEIGHT_PER_ROW = 14 * cm * 3 / 4  # Height per row
# Maximum width (same as multi plot width: 3 columns)
MAX_FIGURE_WIDTH = 14 * cm  # Maximum width for all figures (matches multi plot)
# Reduced height per row for better vertical spacing (reduces gap between rows)
HEIGHT_REDUCTION_FACTOR = 0.85  # Reduce height by 15% to bring plots closer


def get_figure_size(ncols: int, nrows: int) -> Tuple[float, float]:
    """
    Calculate figure size based on number of columns and rows.
    
    All figures use the same maximum width as the multi plot (14 * cm).
    Width is capped at MAX_FIGURE_WIDTH regardless of number of columns.
    Height is reduced to improve vertical spacing between plots.
    
    Args:
        ncols: Number of columns
        nrows: Number of rows
        
    Returns:
        Tuple of (width, height) in inches
    """
    # Use maximum width for all figures (same as multi plot)
    fig_width = MAX_FIGURE_WIDTH
    # Apply height reduction factor to bring plots closer together vertically
    fig_height = BASE_HEIGHT_PER_ROW * nrows * HEIGHT_REDUCTION_FACTOR
    return (fig_width, fig_height)

# Colour palette (matching reference code)
maxRGB = 255
colours = np.array([
    [0, 61 / maxRGB, 115 / maxRGB],  # navy blue
    [226 / maxRGB, 0, 26 / maxRGB],  # red
    [238 / maxRGB, 127 / maxRGB, 0],  # orange
    [139 / maxRGB, 173 / maxRGB, 63 / maxRGB],  # green
    [101 / maxRGB, 90 / maxRGB, 159 / maxRGB],  # lilla
    [250 / maxRGB, 187 / maxRGB, 0],  # yellow
    [55 / maxRGB, 160 / maxRGB, 203 / maxRGB],  # cyan
    [0, 171 / maxRGB, 164 / maxRGB],  # turkis
    [226 / maxRGB, 0, 112 / maxRGB],  # magenta
    [135 / maxRGB, 135 / maxRGB, 135 / maxRGB],  # grey
])


def get_cell_volumes(mesh: pv.PolyData) -> np.ndarray:
    """
    Get cell volumes from PyVista mesh.

    Args:
        mesh: PyVista mesh

    Returns:
        Cell volumes array
    """
    # Try to get volumes from cell data first
    if "V" in mesh.cell_data:
        return np.array(mesh.cell_data["V"])
    elif "Volume" in mesh.cell_data:
        return np.array(mesh.cell_data["Volume"])
    else:
        # Compute volumes using PyVista's built-in method
        logger.info("Computing cell volumes from mesh geometry...")
        volumes = np.zeros(mesh.n_cells)
        for i in range(mesh.n_cells):
            cell = mesh.extract_cells(i)
            if cell.n_cells > 0:
                volumes[i] = cell.volume
        return volumes


def point_to_cell_data(
    mesh: pv.PolyData, point_field: np.ndarray
) -> np.ndarray:
    """
    Convert point data to cell data by averaging point values at cell vertices.

    Args:
        mesh: PyVista mesh
        point_field: Field values at points (1D or 2D array)

    Returns:
        Field values at cells (1D or 2D array, matching input shape)
    """
    # Handle both scalar and vector fields
    if point_field.ndim == 1:
        # Scalar field
        cell_data = np.zeros(mesh.n_cells)
        for i in range(mesh.n_cells):
            cell = mesh.get_cell(i)
            cell_point_indices = cell.point_ids
            cell_data[i] = np.mean(point_field[cell_point_indices])
        return cell_data
    else:
        # Vector field (n_points, n_components)
        n_components = point_field.shape[1]
        cell_data = np.zeros((mesh.n_cells, n_components))
        for i in range(mesh.n_cells):
            cell = mesh.get_cell(i)
            cell_point_indices = cell.point_ids
            cell_data[i] = np.mean(point_field[cell_point_indices], axis=0)
        return cell_data


def find_latest_time_folder(directory: Path) -> Optional[str]:
    """
    Find the latest time folder in the directory.

    Args:
        directory: Directory path to search

    Returns:
        Name of the latest time folder, or None if no numeric folders found
    """
    numeric_folders = []
    for item in directory.iterdir():
        if item.is_dir() and item.name != "0" and item.name.isnumeric():
            numeric_folders.append(int(item.name))

    if not numeric_folders:
        return None

    return str(max(numeric_folders))


def find_foam_file(case_dir: Path) -> Path:
    """
    Find the .foam file in the case directory.

    Args:
        case_dir: Case directory path

    Returns:
        Path to the .foam file

    Raises:
        FileNotFoundError: If no .foam file is found
    """
    foam_files = list(case_dir.glob("*.foam"))
    if not foam_files:
        raise FileNotFoundError(f"No .foam file found in {case_dir}")

    logger.info(f"Using foam file: {foam_files[0].name}")
    return foam_files[0]


def reorder_pyvista_tensor_to_openfoam(rij: np.ndarray) -> np.ndarray:
    """
    Reorder symmetric tensor from PyVista format to OpenFOAM format.
    
    PyVista format: [Rxx, Ryy, Rzz, Rxy, Ryz, Rxz] -> indices [0, 1, 2, 3, 4, 5]
    OpenFOAM format: [Rxx, Rxy, Rxz, Ryy, Ryz, Rzz] -> indices [0, 1, 2, 3, 4, 5]
    
    Mapping: [0, 3, 5, 1, 4, 2]
    - Rxx: PyVista[0] -> OpenFOAM[0]
    - Rxy: PyVista[3] -> OpenFOAM[1]
    - Rxz: PyVista[5] -> OpenFOAM[2]
    - Ryy: PyVista[1] -> OpenFOAM[3]
    - Ryz: PyVista[4] -> OpenFOAM[4]
    - Rzz: PyVista[2] -> OpenFOAM[5]
    
    Args:
        rij: Reynolds stress tensor in PyVista format
             Shape: (n_points, 6) or (6, n_points) or (n_points, 9)
    
    Returns:
        Rij tensor in OpenFOAM format with same shape
    """
    if rij.ndim == 1:
        # Single tensor: (6,)
        if rij.shape[0] == 6:
            return rij[[0, 3, 5, 1, 4, 2]]
        elif rij.shape[0] == 9:
            # Full tensor, extract symmetric components first
            rij_6 = np.array([rij[0], rij[4], rij[8], rij[1], rij[5], rij[2]])
            return rij_6[[0, 3, 5, 1, 4, 2]]
        else:
            return rij
    elif rij.ndim == 2:
        if rij.shape[1] == 6:
            # Row-based: (n_points, 6)
            return rij[:, [0, 3, 5, 1, 4, 2]]
        elif rij.shape[0] == 6:
            # Column-based: (6, n_points)
            return rij[[0, 3, 5, 1, 4, 2], :]
        elif rij.shape[1] == 9:
            # Full tensor: (n_points, 9), extract symmetric components first
            rij_6 = np.column_stack([
                rij[:, 0], rij[:, 4], rij[:, 8],  # Rxx, Ryy, Rzz
                rij[:, 1], rij[:, 5], rij[:, 2],  # Rxy, Ryz, Rxz
            ])
            return rij_6[:, [0, 3, 5, 1, 4, 2]]
        else:
            return rij
    else:
        return rij


def load_openfoam_mesh(
    case_dir: Path, time_value: str = "latest_time"
) -> Tuple[pv.POpenFOAMReader, pv.MultiBlock, pv.PolyData]:
    """
    Load OpenFOAM mesh using PyVista.

    Args:
        case_dir: Path to OpenFOAM case directory
        time_value: Time value to load ('latest_time' or numeric string)

    Returns:
        Tuple of (reader, mesh MultiBlock, internal mesh PolyData)
    """
    foam_file = find_foam_file(case_dir)

    logger.info(f"Loading OpenFOAM case from: {foam_file}")
    reader = pv.POpenFOAMReader(str(foam_file))

    # Get available time values
    time_values = reader.time_values
    logger.info(f"Available time values: {time_values}")

    # Select time value
    if time_value == "latest_time":
        selected_time = time_values[-1]
    else:
        selected_time = float(time_value)

    logger.info(f"Using time value: {selected_time}")

    # Configure reader and load mesh
    reader.set_active_time_value(selected_time)
    reader.cell_to_point_creation = True

    mesh = reader.read()
    internal_mesh = mesh["internalMesh"]

    logger.info(
        f"Mesh loaded with {internal_mesh.n_points} points "
        f"and {internal_mesh.n_cells} cells"
    )

    return reader, mesh, internal_mesh


def slice_mesh_at_x_plane(
    mesh: pv.PolyData, x_value: Optional[float] = None
) -> pv.PolyData:
    """
    Slice mesh at a plane perpendicular to x-axis.

    Args:
        mesh: PyVista mesh to slice
        x_value: X-coordinate value for the slice plane.
                 If None, uses the midpoint between min and max x.

    Returns:
        Sliced mesh
    """
    # Calculate x midpoint if not provided
    if x_value is None:
        points = mesh.points
        x_min = np.min(points[:, 0])
        x_max = np.max(points[:, 0])
        x_value = (x_min + x_max) / 2.0
        logger.info(
            f"Calculated x midpoint: {x_value:.6f} (range: {x_min:.6f} to {x_max:.6f})"
        )

    # Create slice origin at x=x_value, using mesh centre for y and z
    centre = np.array(mesh.center)
    centre[0] = x_value

    logger.info(f"Slicing mesh at x={x_value:.6f} with normal: x")
    sliced = mesh.slice(normal="x", origin=centre)

    logger.info(f"Sliced mesh has {sliced.n_points} points")
    return sliced


def get_field_from_mesh(
    mesh: pv.PolyData, field_name: str, component: Optional[int] = None
) -> np.ndarray:
    """
    Get field from mesh point_data or cell_data.

    Args:
        mesh: PyVista mesh
        field_name: Name of the field
        component: Optional component index for vector fields (0, 1, or 2)

    Returns:
        Field data array

    Raises:
        KeyError: If field is not found
    """
    # Try point data first
    if field_name in mesh.point_data:
        field = np.array(mesh.point_data[field_name])
        if field.ndim > 1 and component is not None:
            return field[:, component]
        return field.flatten() if field.ndim > 1 else field

    # Try cell data
    if field_name in mesh.cell_data:
        field = np.array(mesh.cell_data[field_name])
        if field.ndim > 1 and component is not None:
            return field[:, component]
        return field.flatten() if field.ndim > 1 else field

    raise KeyError(f"Field '{field_name}' not found in mesh point_data or cell_data")


def interpolate_field_to_mesh(
    source_mesh: pv.PolyData,
    target_mesh: pv.PolyData,
    field_name: str,
    component: Optional[int] = None,
) -> np.ndarray:
    """
    Interpolate a field from source mesh to target mesh.

    Args:
        source_mesh: Source mesh with the field
        target_mesh: Target mesh to interpolate to
        field_name: Name of the field to interpolate
        component: Optional component index for vector fields

    Returns:
        Interpolated field on target mesh
    """
    # Get field from source mesh
    source_field = get_field_from_mesh(source_mesh, field_name, component)

    # Get source points
    source_points = source_mesh.points

    # Get target points (use cell centers if field is cell data)
    if field_name in source_mesh.cell_data:
        target_points = target_mesh.cell_centers().points
    else:
        target_points = target_mesh.points

    # Use nearest neighbor interpolation
    distances = cdist(target_points, source_points)
    nearest_indices = np.argmin(distances, axis=1)

    return source_field[nearest_indices]


def compute_error_metrics(
    current_field: np.ndarray,
    reference_field: np.ndarray,
    volume: np.ndarray,
    baseline_error: float,
) -> float:
    """
    Compute normalised error metric.

    Args:
        current_field: Current simulation field
        reference_field: Reference (LES) field
        volume: Cell volumes for volumetric averaging
        baseline_error: Baseline error for normalisation

    Returns:
        Normalised error metric
    """
    diff = np.abs(current_field - reference_field)
    vol_avg_diff = volumetric_average(diff, volume)
    return vol_avg_diff / baseline_error if baseline_error > 0 else 1000.0


def compute_rmse(
    current_field: np.ndarray,
    reference_field: np.ndarray,
    volume: np.ndarray,
) -> float:
    """
    Compute Root Mean Square Error (RMSE) with volumetric weighting.

    Args:
        current_field: Current simulation field
        reference_field: Reference (HF) field
        volume: Cell volumes for volumetric averaging

    Returns:
        Volumetrically-weighted RMSE
    """
    diff_squared = (current_field - reference_field) ** 2
    vol_avg_diff_squared = volumetric_average(diff_squared, volume)
    return np.sqrt(vol_avg_diff_squared)


def compute_mae(
    current_field: np.ndarray,
    reference_field: np.ndarray,
    volume: np.ndarray,
) -> float:
    """
    Compute Mean Absolute Error (MAE) with volumetric weighting.

    Args:
        current_field: Current simulation field
        reference_field: Reference (HF) field
        volume: Cell volumes for volumetric averaging

    Returns:
        Volumetrically-weighted MAE
    """
    diff = np.abs(current_field - reference_field)
    return volumetric_average(diff, volume)


def compute_relative_rmse(
    current_field: np.ndarray,
    reference_field: np.ndarray,
    volume: np.ndarray,
) -> float:
    """
    Compute relative RMSE (normalised by reference field magnitude).

    Args:
        current_field: Current simulation field
        reference_field: Reference (HF) field
        volume: Cell volumes for volumetric averaging

    Returns:
        Relative RMSE (as percentage)
    """
    rmse = compute_rmse(current_field, reference_field, volume)
    ref_rms = np.sqrt(volumetric_average(reference_field ** 2, volume))
    return (rmse / ref_rms * 100.0) if ref_rms > 0 else 1000.0


def compute_relative_mae(
    current_field: np.ndarray,
    reference_field: np.ndarray,
    volume: np.ndarray,
) -> float:
    """
    Compute relative MAE (normalised by reference field magnitude).

    Args:
        current_field: Current simulation field
        reference_field: Reference (HF) field
        volume: Cell volumes for volumetric averaging

    Returns:
        Relative MAE (as percentage)
    """
    mae = compute_mae(current_field, reference_field, volume)
    ref_mean_abs = volumetric_average(np.abs(reference_field), volume)
    return (mae / ref_mean_abs * 100.0) if ref_mean_abs > 0 else 1000.0


def volumetric_average(
    field: np.ndarray, volume: np.ndarray
) -> float:
    """
    Compute volumetric average of a field.

    Args:
        field: Field values (1D array)
        volume: Cell volumes (1D array)

    Returns:
        Volumetric average value
    """
    return np.sum(np.multiply(np.squeeze(field), volume)) / np.sum(volume)


def save_objective_function_results(
    output_path: Path,
    error_velocity: float,
    error_vorticity: float,
    error_combined: float,
    rmse_velocity: Optional[float] = None,
    rmse_vorticity: Optional[float] = None,
    mae_velocity: Optional[float] = None,
    mae_vorticity: Optional[float] = None,
    relative_rmse_velocity: Optional[float] = None,
    relative_rmse_vorticity: Optional[float] = None,
    relative_mae_velocity: Optional[float] = None,
    relative_mae_vorticity: Optional[float] = None,
    rmse_velocity_kosst: Optional[float] = None,
    rmse_vorticity_kosst: Optional[float] = None,
    mae_velocity_kosst: Optional[float] = None,
    mae_vorticity_kosst: Optional[float] = None,
    converged: bool = True,
) -> None:
    """
    Save objective function results to JSON file in human-readable format.

    Args:
        output_path: Path to output JSON file
        error_velocity: Velocity error metric (j1) - normalised vs kOmegaSST baseline
        error_vorticity: Vorticity error metric (j2) - normalised vs kOmegaSST baseline
        error_combined: Combined error metric (J)
        rmse_velocity: Root Mean Square Error for velocity vs HF
        rmse_vorticity: Root Mean Square Error for vorticity vs HF
        mae_velocity: Mean Absolute Error for velocity vs HF
        mae_vorticity: Mean Absolute Error for vorticity vs HF
        relative_rmse_velocity: Relative RMSE for velocity vs HF (%)
        relative_rmse_vorticity: Relative RMSE for vorticity vs HF (%)
        relative_mae_velocity: Relative MAE for velocity vs HF (%)
        relative_mae_vorticity: Relative MAE for vorticity vs HF (%)
        rmse_velocity_kosst: RMSE for kOmegaSST velocity vs HF (for comparison)
        rmse_vorticity_kosst: RMSE for kOmegaSST vorticity vs HF (for comparison)
        mae_velocity_kosst: MAE for kOmegaSST velocity vs HF (for comparison)
        mae_vorticity_kosst: MAE for kOmegaSST vorticity vs HF (for comparison)
    """
    # Build results dictionary with values first, then descriptions
    results = {}
    descriptions = {}
    
    # Add all numeric values first
    results["j1"] = error_velocity
    results["j2"] = error_vorticity
    results["J"] = error_combined
    descriptions["j1"] = "Normalised velocity error metric (vs kOmegaSST baseline)"
    descriptions["j2"] = "Normalised vorticity error metric (vs kOmegaSST baseline)"
    descriptions["J"] = "Combined error metric (average of j1 and j2)"
    
    # Add RMSE metrics if provided
    if rmse_velocity is not None:
        results["rmse_velocity"] = rmse_velocity
        descriptions["rmse_velocity"] = "Root Mean Square Error for velocity vs HF (m/s)"
    
    if rmse_vorticity is not None:
        results["rmse_vorticity"] = rmse_vorticity
        descriptions["rmse_vorticity"] = "Root Mean Square Error for vorticity vs HF (1/s)"
    
    # Add MAE metrics if provided
    if mae_velocity is not None:
        results["mae_velocity"] = mae_velocity
        descriptions["mae_velocity"] = "Mean Absolute Error for velocity vs HF (m/s)"
    
    if mae_vorticity is not None:
        results["mae_vorticity"] = mae_vorticity
        descriptions["mae_vorticity"] = "Mean Absolute Error for vorticity vs HF (1/s)"
    
    # Add relative RMSE metrics if provided
    if relative_rmse_velocity is not None:
        results["relative_rmse_velocity"] = relative_rmse_velocity
        descriptions["relative_rmse_velocity"] = "Relative RMSE for velocity vs HF (absolute percentage, 0-100 scale)"
    
    if relative_rmse_vorticity is not None:
        results["relative_rmse_vorticity"] = relative_rmse_vorticity
        descriptions["relative_rmse_vorticity"] = "Relative RMSE for vorticity vs HF (absolute percentage, 0-100 scale)"
    
    # Add relative MAE metrics if provided
    if relative_mae_velocity is not None:
        results["relative_mae_velocity"] = relative_mae_velocity
        descriptions["relative_mae_velocity"] = "Relative MAE for velocity vs HF (absolute percentage, 0-100 scale)"
    
    if relative_mae_vorticity is not None:
        results["relative_mae_vorticity"] = relative_mae_vorticity
        descriptions["relative_mae_vorticity"] = "Relative MAE for vorticity vs HF (absolute percentage, 0-100 scale)"
    
    # Add kOmegaSST baseline metrics for comparison if provided
    if rmse_velocity_kosst is not None:
        results["rmse_velocity_kosst"] = rmse_velocity_kosst
        descriptions["rmse_velocity_kosst"] = "RMSE for kOmegaSST velocity vs HF (m/s) - baseline"
    
    if rmse_vorticity_kosst is not None:
        results["rmse_vorticity_kosst"] = rmse_vorticity_kosst
        descriptions["rmse_vorticity_kosst"] = "RMSE for kOmegaSST vorticity vs HF (1/s) - baseline"
    
    if mae_velocity_kosst is not None:
        results["mae_velocity_kosst"] = mae_velocity_kosst
        descriptions["mae_velocity_kosst"] = "MAE for kOmegaSST velocity vs HF (m/s) - baseline"
    
    if mae_vorticity_kosst is not None:
        results["mae_vorticity_kosst"] = mae_vorticity_kosst
        descriptions["mae_vorticity_kosst"] = "MAE for kOmegaSST vorticity vs HF (1/s) - baseline"
    
    # Add converged flag
    results["converged"] = converged
    descriptions["converged"] = "Whether the simulation converged (False if max iterations reached)"
    
    # Add descriptions section at the end
    results["description"] = descriptions
    
    with open(output_path, "w") as f:
        json.dump(results, f, indent=4)
    logger.info(f"Objective function results saved to: {output_path}")


def barycentric_calculation_hf(
    Rij: np.ndarray,
    Coff: float = 0.65,
    Cexp: float = 5
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate barycentric coordinates (eta, zeta) and colours from Reynolds stress tensor for HF data.
    
    Matches baryCentricCalculationHF reference implementation.
    Uses column-based indexing (Rij shape is (6, n_points)).

    Args:
        Rij: Reynolds stress tensor (6, n_points) where rows are:
             [Rxx, Rxy, Rxz, Ryy, Ryz, Rzz]
        Coff: Colour offset for enhanced colourmap (default: 0.65)
        Cexp: Colour exponent for enhanced colourmap (default: 5)

    Returns:
        Tuple of (eta, zeta, colours) arrays
    """
    num_elements = Rij.shape[1]
    
    # TKE: k = 0.5 * trace(Rij)
    k = 0.5 * (Rij[0] + Rij[3] + Rij[5])
    
    # Build anisotropy tensor bij: bij = Rij / (2k) - I/3
    bij = np.zeros((num_elements, 3, 3))
    bij[:, 0, 0] = (Rij[0] / k / 2 - 1 / 3).T
    bij[:, 0, 1] = (Rij[1] / k / 2).T
    bij[:, 0, 2] = (Rij[2] / k / 2).T
    bij[:, 1, 0] = (Rij[1] / k / 2).T
    bij[:, 1, 1] = (Rij[3] / k / 2 - 1 / 3).T
    bij[:, 1, 2] = (Rij[4] / k / 2).T
    bij[:, 2, 0] = (Rij[2] / k / 2).T
    bij[:, 2, 1] = (Rij[4] / k / 2).T
    bij[:, 2, 2] = (Rij[5] / k / 2 - 1 / 3).T
    
    # Calculate eigenvalues
    eigvals = np.linalg.eigvals(bij)
    
    # Sort eigenvalues in descending order for each element
    indices = np.argsort(-eigvals)
    w = np.zeros((num_elements, 3))
    for i in range(num_elements):
        w[i] = eigvals[i][indices[i]]
    
    alpha = np.array([w[:, 0], w[:, 1], w[:, 2]])
    
    # Calculate barycentric coefficients
    C_1cc = alpha[0, :] - alpha[1, :]
    C_2cc = 2 * (alpha[1, :] - alpha[2, :])
    C_3cc = 3 * alpha[2, :] + 1
    
    # Consistency check
    checkC = C_1cc + C_2cc + C_3cc
    if (checkC - 1 > 1e-3).any():
        logger.warning("Reynolds Stress Tensor not consistent (C1+C2+C3 != 1)")
    
    # Calculate colours (Emory & Iaccarino enhanced colourmap)
    colour_prime = np.hstack((
        (C_1cc.reshape(-1, 1) + Coff) ** Cexp,
        (C_2cc.reshape(-1, 1) + Coff) ** Cexp,
        (C_3cc.reshape(-1, 1) + Coff) ** Cexp
    ))
    colours = np.minimum(colour_prime, 1)
    
    # Calculate barycentric coordinates (eta, zeta)
    eta = C_1cc + C_3cc * 0.5
    zeta = C_3cc * np.sqrt(3) / 2
    
    return eta, zeta, colours


def barycentric_calculation_kOmegaSST(
    Rij: np.ndarray,
    Coff: float = 0.65,
    Cexp: float = 5
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate barycentric coordinates (eta, zeta) and colours from Reynolds stress tensor.
    
    Matches baryCentricCalculationOfpp reference implementation.
    Uses row-based indexing: Rij shape is (n_points, 6).

    Args:
        Rij: Reynolds stress tensor (n_points, 6) where columns are:
             [Rxx, Rxy, Rxz, Ryy, Ryz, Rzz]
        Coff: Colour offset for enhanced colourmap (default: 0.65)
        Cexp: Colour exponent for enhanced colourmap (default: 5)

    Returns:
        Tuple of (eta, zeta, colours) arrays
    """
    num_elements = Rij.shape[0]
    
    # TKE: k = 0.5 * trace(Rij)
    k = 0.5 * (Rij[:, 0] + Rij[:, 3] + Rij[:, 5])
    
    # Build anisotropy tensor bij: bij = Rij / (2k) - I/3
    bij = np.zeros((num_elements, 3, 3))
    bij[:, 0, 0] = (Rij[:, 0] / k / 2 - 1 / 3).T
    bij[:, 0, 1] = (Rij[:, 1] / k / 2).T
    bij[:, 0, 2] = (Rij[:, 2] / k / 2).T
    bij[:, 1, 0] = (Rij[:, 1] / k / 2).T
    bij[:, 1, 1] = (Rij[:, 3] / k / 2 - 1 / 3).T
    bij[:, 1, 2] = (Rij[:, 4] / k / 2).T
    bij[:, 2, 0] = (Rij[:, 2] / k / 2).T
    bij[:, 2, 1] = (Rij[:, 4] / k / 2).T
    bij[:, 2, 2] = (Rij[:, 5] / k / 2 - 1 / 3).T
    
    # Calculate eigenvalues
    eigvals = np.linalg.eigvals(bij)
    
    # Sort eigenvalues in descending order for each element
    indices = np.argsort(-eigvals)
    w = np.zeros((num_elements, 3))
    for i in range(num_elements):
        w[i] = eigvals[i][indices[i]]
    
    alpha = np.array([w[:, 0], w[:, 1], w[:, 2]])
    
    # Calculate barycentric coefficients
    C_1cc = alpha[0, :] - alpha[1, :]
    C_2cc = 2 * (alpha[1, :] - alpha[2, :])
    C_3cc = 3 * alpha[2, :] + 1
    
    # Consistency check
    checkC = C_1cc + C_2cc + C_3cc
    if (checkC - 1 > 1e-3).any():
        logger.warning("Reynolds Stress Tensor not consistent (C1+C2+C3 != 1)")
    
    # Calculate colours (Emory & Iaccarino enhanced colourmap)
    colour_prime = np.hstack((
        (C_1cc.reshape(-1, 1) + Coff) ** Cexp,
        (C_2cc.reshape(-1, 1) + Coff) ** Cexp,
        (C_3cc.reshape(-1, 1) + Coff) ** Cexp
    ))
    colours = np.minimum(colour_prime, 1)
    
    # Calculate barycentric coordinates (eta, zeta)
    eta = C_1cc + C_3cc * 0.5
    zeta = C_3cc * np.sqrt(3) / 2
    
    return eta, zeta, colours


def plot_barycentric_triangle(ax: plt.Axes) -> None:
    """
    Draw the barycentric triangle frame on an axes.

    Args:
        ax: Matplotlib axes to draw on
    """
    # Triangle vertices
    a_1cc = np.array([1, 0])
    a_2cc = np.array([0, 0])
    a_3cc = np.array([0.5, np.sqrt(3) / 2])

    # Triangle outline
    triangle_x = np.array([a_1cc[0], a_2cc[0], a_3cc[0], a_1cc[0]])
    triangle_y = np.array([a_1cc[1], a_2cc[1], a_3cc[1], a_1cc[1]])

    # Plane strain line
    strain_x = np.array([a_3cc[0], 1 / 3])
    strain_y = np.array([a_3cc[1], 0])

    ax.axis('off')
    ax.plot(triangle_x, triangle_y, lw=2, c='k')
    ax.plot(strain_x, strain_y, lw=1, ls=':', c='k')
    
    # Vertex labels
    ax.text(1.05, -0.05, r'$\hat{x}_{1c}$', size=7)
    ax.text(-0.1, -0.05, r'$\hat{x}_{2c}$', size=7)
    ax.text(0.47, a_3cc[1] + 0.025, r'$\hat{x}_{3c}$', size=7)
    ax.set_aspect('equal')


def create_barycentric_plot(
    eta_kosst: np.ndarray,
    zeta_kosst: np.ndarray,
    colours_kosst: np.ndarray,
    eta_current: np.ndarray,
    zeta_current: np.ndarray,
    colours_current: np.ndarray,
    eta_hf: np.ndarray,
    zeta_hf: np.ndarray,
    colours_hf: np.ndarray,
    output_path: Path,
) -> None:
    """
    Create barycentric triangle plots comparing kOmegaSST, current model, and HF.

    Args:
        eta_kosst: Eta coordinates for kOmegaSST
        zeta_kosst: Zeta coordinates for kOmegaSST
        colours_kosst: Colours for kOmegaSST points
        eta_current: Eta coordinates for current model
        zeta_current: Zeta coordinates for current model
        colours_current: Colours for current model points
        eta_hf: Eta coordinates for HF
        zeta_hf: Zeta coordinates for HF
        colours_hf: Colours for HF points
        output_path: Path to save the plot
    """
    # Create figure with 3 subplots (1 row × 3 cols)
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=get_figure_size(3, 1), constrained_layout=False
    )
    
    # Draw triangle frames
    plot_barycentric_triangle(ax1)
    plot_barycentric_triangle(ax2)
    plot_barycentric_triangle(ax3)
    
    # kOmegaSST plot
    if len(eta_kosst) > 0:
        # Ensure colours are in [0, 1] range
        colours_kosst_safe = np.clip(colours_kosst, 0, 1)
        ax1.scatter(
            eta_kosst, zeta_kosst, marker='.', s=5,
            c=colours_kosst_safe, alpha=0.1
        )
    ax1.set_title(r"$k-\omega \ \text{SST}$", fontsize=10)
    
    # Current model plot
    if len(eta_current) > 0:
        # Ensure colours are in [0, 1] range
        colours_current_safe = np.clip(colours_current, 0, 1)
        ax2.scatter(
            eta_current, zeta_current, marker='.', s=5,
            c=colours_current_safe, alpha=0.1
        )
    ax2.set_title(r"$k-\omega \ \text{SST-PDA}$", fontsize=10)
    
    # HF plot
    if len(eta_hf) > 0:
        # Ensure colours are in [0, 1] range
        colours_hf_safe = np.clip(colours_hf, 0, 1)
        ax3.scatter(
            eta_hf, zeta_hf, marker='.', s=5,
            c=colours_hf_safe, alpha=0.1
        )
    ax3.set_title("High Fidelity", fontsize=10)
    
    # Save figure with improved spacing
    fig.tight_layout(pad=1.0, h_pad=0.5, w_pad=0.5)
    fig.savefig(output_path, dpi=600, transparent=False, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Barycentric plot saved to: {output_path}")


def create_anisotropy_contour_plot(
    sliced_mesh: pv.PolyData,
    rij_kosst: np.ndarray,
    rij_current: np.ndarray,
    rij_hf: np.ndarray,
    output_path: Path,
    resolution: int = 100,
) -> None:
    """
    Create anisotropy contour plots showing barycentric colours in physical space.

    This plots the y-z plane with each point coloured according to its
    turbulence anisotropy state (using the Emory & Iaccarino colourmap).

    Args:
        sliced_mesh: Sliced PyVista mesh (y-z plane)
        rij_kosst: Reynolds stress tensor for k-omega SST (n_points, 6)
        rij_current: Reynolds stress tensor for k-omega SST-PDA (n_points, 6)
        rij_hf: Reynolds stress tensor for HF (n_points, 6)
        output_path: Path to save the plot
        resolution: Grid resolution for interpolation
    """
    # Get coordinates from sliced mesh
    points = sliced_mesh.points
    y_coords = points[:, 1]
    z_coords = points[:, 2]

    # Create regular grid for interpolation
    y_grid = np.linspace(y_coords.min(), y_coords.max(), resolution)
    z_grid = np.linspace(z_coords.min(), z_coords.max(), resolution)
    yg, zg = np.meshgrid(y_grid, z_grid, indexing='ij')

    # Points for interpolation
    points_2d = np.column_stack((y_coords, z_coords))

    # Figure setup (1 row × 3 cols)
    fig, axes = plt.subplots(
        figsize=get_figure_size(3, 1),
        ncols=3,
        nrows=1,
        sharex=True,
        sharey=True,
    )

    titles = [r'$k-\omega \ \text{SST}$', r'$k-\omega \ \text{SST-PDA}$', 'High Fidelity']
    rij_data = [rij_kosst, rij_current, rij_hf]

    for idx, (ax, title, rij) in enumerate(zip(axes, titles, rij_data)):
        if rij is None or len(rij) == 0:
            ax.set_title(title)
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            continue

        # Interpolate each component of Rij onto the regular grid
        rij_interp = np.zeros((resolution, resolution, 6))
        for comp in range(6):
            rij_interp[:, :, comp] = griddata(
                points_2d, rij[:, comp], (yg, zg), method='linear'
            )

        # Compute barycentric colours for each grid point
        for i in range(resolution):
            # Get Rij values along this row (shape: resolution, 6)
            rij_row = rij_interp[i, :, :]

            # Skip rows with NaN values
            valid_mask = ~np.isnan(rij_row).any(axis=1)
            if not valid_mask.any():
                continue

            # Compute barycentric coordinates and colours for valid points
            rij_valid = rij_row[valid_mask]
            _, _, colours = barycentric_calculation_kOmegaSST(rij_valid)

            # Clip colours to valid range
            colours_safe = np.clip(colours, 0, 1)

            # Plot scatter points with barycentric colours
            ax.scatter(
                yg[i, valid_mask] / H,
                zg[i, valid_mask] / H,
                marker='.',
                s=1,
                c=colours_safe,
                alpha=1,
                rasterized=True,
            )

        ax.set_title(title, fontsize=10)
        ax.set_aspect('equal', 'box')

        # Set axis limits and ticks
        y_min, y_max = y_coords.min() / H, y_coords.max() / H
        z_min, z_max = z_coords.min() / H, z_coords.max() / H

        ax.set_xlim([y_min, y_max])
        ax.set_ylim([z_min, z_max])

        # Set ticks
        delta = 0.5
        xticks = np.arange(np.ceil(y_min / delta) * delta, y_max + delta, delta)
        yticks = np.arange(np.ceil(z_min / delta) * delta, z_max + delta, delta)
        ax.set_xticks(xticks)
        ax.set_xticklabels([f'{t:.1f}' for t in xticks])

    # Set y-axis ticks only for first subplot
    axes[0].set_yticks(yticks)
    axes[0].set_yticklabels([f'{t:.1f}' for t in yticks])

    # Labels
    fig.supxlabel(r'$y/H$', fontsize=10)
    axes[0].set_ylabel(r'$z/H$', fontsize=10)

    # Add colour legend (triangle reference)
    # Create a small inset axes for the triangle legend
    legend_ax = fig.add_axes([0.92, 0.3, 0.08, 0.4])
    plot_barycentric_triangle(legend_ax)
    legend_ax.set_title('Anisotropy', fontsize=7)

    # Save figure with improved spacing
    fig.tight_layout(rect=[0, 0, 0.9, 1], pad=1.0, h_pad=0.3, w_pad=0.5)
    fig.savefig(output_path, dpi=600, transparent=False, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Anisotropy contour plot saved to: {output_path}")


def create_velocity_profiles_plot(
    internal_mesh: pv.PolyData,
    ref_internal_mesh: pv.PolyData,
    u_current: np.ndarray,
    u_kosst: np.ndarray,
    u_hf: np.ndarray,
    output_path: Path,
    profiles: np.ndarray = np.array([0.25, 0.5, 0.75]),
    scale: np.ndarray = np.array([1, 10, 10]),
    Ny: int = 100,
    NyDNS: int = 50,
    AR: float = 1.0,
) -> None:
    """
    Create velocity profile plots at specific y/H locations.

    Args:
        internal_mesh: Current simulation mesh
        ref_internal_mesh: Reference mesh (for kOmegaSST and HF data)
        u_current: Current velocity field (n_points, 3) [u, v, w]
        u_kosst: kOmegaSST velocity field (n_points, 3) [u, v, w]
        u_hf: High Fidelity velocity field (n_points, 3) [u, v, w]
        output_path: Path to save the plot
        profiles: y/H locations for profiles (default: [0.25, 0.5, 0.75])
        scale: Scaling factors for each velocity component (default: [1, 10, 10])
        Ny: Number of points for RANS profiles
        NyDNS: Number of points for DNS/HF profiles
        AR: Aspect ratio (default: 1.0)
    """
    # Get cell centers and coordinates
    current_points = internal_mesh.cell_centers().points
    ref_points = ref_internal_mesh.cell_centers().points

    # Normalize coordinates
    C_y_current = current_points[:, 1] / H
    C_z_current = current_points[:, 2] / H
    C_y_ref = ref_points[:, 1] / H
    C_z_ref = ref_points[:, 2] / H

    # Compute bulk velocity from HF data
    if "V" in ref_internal_mesh.cell_data:
        volume = np.array(ref_internal_mesh.cell_data["V"])
    elif "Volume" in ref_internal_mesh.cell_data:
        volume = np.array(ref_internal_mesh.cell_data["Volume"])
    else:
        # Compute volumes if not available
        volume = np.ones(ref_internal_mesh.n_cells)
        logger.warning("Volume data not found, using uniform weights for bulk velocity")

    Ub = volumetric_average(u_hf[:, 0], volume)

    # Create figure (1 row × 3 cols)
    fig, axes = plt.subplots(
        figsize=get_figure_size(3, 1),
        ncols=3,
        nrows=1,
        sharey=True,
        constrained_layout=False,
    )

    # Line styles and colours
    lw = 1
    COLORS = np.array([
        [[0, 0, 0], colours[1], [0, 0, 0]],  # kOmegaSST: black, Current: red, HF: black
        [[0, 0, 0], colours[3], [0, 0, 0]],  # Alternative if needed
    ])
    ls = np.array(['--', ':'])

    # Velocity component names
    comp_names = ['x', 'y', 'z']
    comp_labels = [
        r"$\langle u_{1} \rangle/u_{b} + y/H$",
        r"$10\langle u_{2} \rangle/u_{b} + y/H$",
        r"$10\langle u_{3} \rangle/u_{b} + y/H$",
    ]

    # Generate y profiles using beta distribution
    yProfile = (beta.cdf(np.linspace(0, 1, Ny), 3, 1) - 1) * AR
    yProfileDNS = (beta.cdf(np.linspace(0, 1, NyDNS), 3, 1) - 1) * AR

    # Plot for each velocity component
    for k, comp_name in enumerate(comp_names):
        ax = axes[k]

        # Get velocity components for each model
        u_comp_kosst = u_kosst[:, k]
        u_comp_current = u_current[:, k]
        u_comp_hf = u_hf[:, k]

        # Plot profiles at each y/H location
        for i, profile_y in enumerate(profiles):
            # kOmegaSST profile
            interp_kosst = LinearNDInterpolator(
                list(zip(C_y_ref, C_z_ref)), u_comp_kosst
            )
            xProfile_kosst = np.ones(Ny) * profile_y
            uProfile_kosst = interp_kosst(xProfile_kosst, yProfile)
            # Remove NaN values
            valid_kosst = ~np.isnan(uProfile_kosst)
            if valid_kosst.any():
                ax.plot(
                    uProfile_kosst[valid_kosst] * scale[k] / Ub - profile_y,
                    yProfile[valid_kosst],
                    c=COLORS[0][0],
                    lw=lw * 0.75,
                    ls='-',
                    zorder=2,
                )

            # Current model profile
            interp_current = LinearNDInterpolator(
                list(zip(C_y_current, C_z_current)), u_comp_current
            )
            xProfile_current = np.ones(Ny) * profile_y
            uProfile_current = interp_current(xProfile_current, yProfile)
            valid_current = ~np.isnan(uProfile_current)
            if valid_current.any():
                ax.plot(
                    uProfile_current[valid_current] * scale[k] / Ub - profile_y,
                    yProfile[valid_current],
                    c=COLORS[0][1],  # red
                    lw=lw,
                    ls=ls[0],  # dashed
                    zorder=3,
                )

            # High Fidelity profile
            interp_hf = LinearNDInterpolator(
                list(zip(C_y_ref, C_z_ref)), u_comp_hf
            )
            xProfile_hf = np.ones(NyDNS) * profile_y
            uProfile_hf = interp_hf(xProfile_hf, yProfileDNS)
            valid_hf = ~np.isnan(uProfile_hf)
            if valid_hf.any():
                ax.plot(
                    uProfile_hf[valid_hf] * scale[k] / Ub - profile_y,
                    yProfileDNS[valid_hf],
                    c=COLORS[0][2],
                    lw=lw,
                    marker='o',
                    markevery=1,
                    markerfacecolor='None',
                    markeredgecolor=COLORS[0][2],
                    markersize=2,
                    markeredgewidth=lw / 2,
                    ls='None',
                    zorder=1,
                )

        # Set labels
        ax.set_xlabel(comp_labels[k], fontsize=10)
        ax.set_box_aspect(1)

        # Set ticks
        if k == 0:
            xticks = np.arange(-1, 1 + 0.5, 0.5)
            xticklabels = np.round(np.arange(0, 2 + 0.5, 0.5), 2)
        else:
            xticks = np.arange(-1, 0 + 0.25, 0.25)
            xticklabels = np.round(np.arange(0, 1 + 0.25, 0.25), 2)
        ax.set_xticks(xticks)
        ax.set_xticklabels([f'{t:.2f}' for t in xticklabels])

    # Y-axis settings
    deltay = 0.5
    yticks = np.arange(-AR, 0 + deltay, deltay)
    yticklabels = np.round(np.arange(0, AR + deltay, deltay), 2)
    axes[0].set_yticks(yticks)
    axes[0].set_yticklabels([f'{t:.2f}' for t in yticklabels])
    axes[0].set_ylabel(r"$z/H$", fontsize=10)

    # Add legend
    lines = [
        Line2D(
            [0], [0], marker='o', color='k', markersize=2, ls='None',
            markerfacecolor='None', markeredgewidth=lw / 2, markeredgecolor='k'
        ),
        Line2D([0], [0], color='k', linewidth=lw * 0.75, linestyle='-'),
        Line2D([0], [0], color=colours[1], linewidth=lw, linestyle='--'),
    ]
    labels = ['High Fidelity', r'$k-\omega \ \text{SST}$', r'$k-\omega \ \text{SST-PDA}$']
    fig.legend(lines, labels, ncol=3, frameon=True, fancybox=True, bbox_to_anchor=(0.75, 0.76))

    # Save figure with improved spacing
    fig.tight_layout(pad=1.0, h_pad=0.3, w_pad=0.5)
    fig.savefig(output_path, dpi=1000, transparent=False, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Velocity profiles plot saved to: {output_path}")


def create_diagonal_velocity_profile_plot(
    internal_mesh: pv.PolyData,
    ref_internal_mesh: pv.PolyData,
    u_current: np.ndarray,
    u_kosst: np.ndarray,
    u_hf: np.ndarray,
    output_path: Path,
    u_tau: Optional[float] = None,
    Ny: int = 200,
) -> None:
    """
    Create velocity profile plot along the diagonal (y=z) of the duct.
    
    Plots u₂/u_τ vs y/H along the diagonal profile where y=z.

    Args:
        internal_mesh: Current simulation mesh
        ref_internal_mesh: Reference mesh (for kOmegaSST and HF data)
        u_current: Current velocity field (n_points, 3) [u, v, w]
        u_kosst: kOmegaSST velocity field (n_points, 3) [u, v, w]
        u_hf: High Fidelity velocity field (n_points, 3) [u, v, w]
        output_path: Path to save the plot
        u_tau: Friction velocity (u_τ). If None, will try to load from mesh or use default
        Ny: Number of points along the diagonal profile
    """
    # Get cell centers and coordinates
    current_points = internal_mesh.cell_centers().points
    ref_points = ref_internal_mesh.cell_centers().points

    # Normalize coordinates
    C_y_current = current_points[:, 1] / H
    C_z_current = current_points[:, 2] / H
    C_y_ref = ref_points[:, 1] / H
    C_z_ref = ref_points[:, 2] / H

    # Try to get u_tau from mesh if not provided
    if u_tau is None:
        # Try to load u_tau from reference mesh (HF data)
        if "u_tau" in ref_internal_mesh.cell_data:
            u_tau_array = np.array(ref_internal_mesh.cell_data["u_tau"])
            u_tau = np.mean(u_tau_array)
            logger.info(f"Using u_tau from mesh: {u_tau:.6f}")
        elif "u_tau" in ref_internal_mesh.point_data:
            u_tau_array = np.array(ref_internal_mesh.point_data["u_tau"])
            u_tau = np.mean(u_tau_array)
            logger.info(f"Using u_tau from mesh: {u_tau:.6f}")
        elif "utau" in ref_internal_mesh.cell_data:
            u_tau_array = np.array(ref_internal_mesh.cell_data["utau"])
            u_tau = np.mean(u_tau_array)
            logger.info(f"Using utau from mesh: {u_tau:.6f}")
        elif "utau" in ref_internal_mesh.point_data:
            u_tau_array = np.array(ref_internal_mesh.point_data["utau"])
            u_tau = np.mean(u_tau_array)
            logger.info(f"Using utau from mesh: {u_tau:.6f}")
        else:
            # Calculate u_tau from bulk velocity (approximate for channel flow)
            # u_tau ≈ 0.05 * U_b for typical channel flows
            if "V" in ref_internal_mesh.cell_data:
                volume = np.array(ref_internal_mesh.cell_data["V"])
            elif "Volume" in ref_internal_mesh.cell_data:
                volume = np.array(ref_internal_mesh.cell_data["Volume"])
            else:
                volume = np.ones(ref_internal_mesh.n_cells)
            
            Ub = volumetric_average(u_hf[:, 0], volume)
            u_tau = 0.05 * Ub  # Approximate relationship
            logger.warning(f"u_tau not found in mesh, using approximate value: {u_tau:.6f} (0.05 * U_b)")

    # Create diagonal profile points: from p0(1, -1) to p1(0, 0) in (y, z) coordinates
    # The diagonal follows z = -y
    p0_y, p0_z = 1.0, -1.0  # Start point
    p1_y, p1_z = 0.0, 0.0   # End point
    
    # Generate points along the diagonal from p0 to p1
    diagonal_y = np.linspace(p0_y, p1_y, Ny)
    diagonal_z = np.linspace(p0_z, p1_z, Ny)  # z = -y
    
    # Get u_2 component (y-component, index 1) for each model
    u_2_kosst = u_kosst[:, 1]
    u_2_current = u_current[:, 1]
    u_2_hf = u_hf[:, 1]
    
    # Interpolate u_2 along the diagonal for each model
    interp_kosst = LinearNDInterpolator(
        list(zip(C_y_ref, C_z_ref)), u_2_kosst, fill_value=np.nan
    )
    interp_current = LinearNDInterpolator(
        list(zip(C_y_current, C_z_current)), u_2_current, fill_value=np.nan
    )
    interp_hf = LinearNDInterpolator(
        list(zip(C_y_ref, C_z_ref)), u_2_hf, fill_value=np.nan
    )
    
    # Extract profiles along diagonal
    u_2_profile_kosst = interp_kosst(diagonal_y, diagonal_z)
    u_2_profile_current = interp_current(diagonal_y, diagonal_z)
    u_2_profile_hf = interp_hf(diagonal_y, diagonal_z)

    u_2_profile_kosst = -u_2_profile_kosst
    u_2_profile_current = -u_2_profile_current
    u_2_profile_hf = -u_2_profile_hf
    
    # Reverse arrays so y goes from 0 to 1 for plotting (diagonal_y goes from 1 to 0)
    y_profile_kosst = np.flip(diagonal_y)
    y_profile_current = np.flip(diagonal_y)
    y_profile_hf = np.flip(diagonal_y)
    u_2_profile_kosst = np.flip(u_2_profile_kosst)
    u_2_profile_current = np.flip(u_2_profile_current)
    u_2_profile_hf = np.flip(u_2_profile_hf)

    # Remove NaN values
    valid_kosst = ~np.isnan(u_2_profile_kosst)
    valid_current = ~np.isnan(u_2_profile_current)
    valid_hf = ~np.isnan(u_2_profile_hf)

    # Normalize by u_tau
    u_2_norm_kosst = u_2_profile_kosst[valid_kosst] / u_tau
    u_2_norm_current = u_2_profile_current[valid_current] / u_tau
    u_2_norm_hf = u_2_profile_hf[valid_hf] / u_tau

    y_profile_kosst = y_profile_kosst[valid_kosst]
    y_profile_current = y_profile_current[valid_current]
    y_profile_hf = y_profile_hf[valid_hf]

    # Create figure (single plot)
    fig, ax = plt.subplots(
        figsize=get_figure_size(3, 1),
        ncols=1,
        nrows=1,
        constrained_layout=False,
    )

    # Line styles and colours
    lw = 1
    ax.plot(
        y_profile_hf,
        u_2_norm_hf,
        c='k',
        lw=lw,
        marker='o',
        markevery=max(1, len(y_profile_hf) // 20),
        markerfacecolor='None',
        markeredgecolor='k',
        markersize=3,
        markeredgewidth=lw / 2,
        ls='None',
        label='High Fidelity',
        zorder=1,
    )
    ax.plot(
        y_profile_kosst,
        u_2_norm_kosst,
        c=colours[0],  # navy blue
        lw=lw * 0.75,
        ls='-',
        label=r'$k-\omega \ \text{SST}$',
        zorder=2,
    )
    ax.plot(
        y_profile_current,
        u_2_norm_current,
        c=colours[1],  # red
        lw=lw,
        ls='--',
        label=r'$k-\omega \ \text{SST-PDA}$',
        zorder=3,
    )

    # Set labels
    ax.set_xlabel(r'$y/H$', fontsize=10)
    ax.set_ylabel(r'$\langle u_{2} \rangle/u_{\tau}$', fontsize=10)

    # Set axis limits
    ax.set_xlim([0, 1])
    
    # Add grid
    ax.grid(True, alpha=0.3)

    # Add legend
    ax.legend(frameon=True, fancybox=True, fontsize=9)

    # Save figure with improved spacing
    fig.tight_layout(pad=1.0, h_pad=0.5, w_pad=0.5)
    fig.savefig(output_path, dpi=1000, transparent=False, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Diagonal velocity profile plot saved to: {output_path}")


def create_rij_profiles_plot(
    internal_mesh: pv.PolyData,
    ref_internal_mesh: pv.PolyData,
    rij_current: np.ndarray,
    rij_kosst: np.ndarray,
    rij_hf: np.ndarray,
    output_path: Path,
    profiles: np.ndarray = np.array([0.25, 0.5, 0.75]),
    scale_rij: np.ndarray = np.array([5, 30, 25, 25, 100, 20]),
    Ny: int = 100,
    NyDNS: int = 50,
    AR: float = 1.0,
) -> None:
    """
    Create Reynolds stress tensor component profile plots at specific y/H locations.

    Args:
        internal_mesh: Current simulation mesh
        ref_internal_mesh: Reference mesh (for kOmegaSST and HF data)
        rij_current: Current Rij tensor (n_points, 6) [Rxx, Rxy, Rxz, Ryy, Ryz, Rzz]
        rij_kosst: kOmegaSST Rij tensor (n_points, 6) [Rxx, Rxy, Rxz, Ryy, Ryz, Rzz]
        rij_hf: High Fidelity Rij tensor (n_points, 6) [Rxx, Rxy, Rxz, Ryy, Ryz, Rzz]
        output_path: Path to save the plot
        profiles: y/H locations for profiles (default: [0.25, 0.5, 0.75])
        scale_rij: Scaling factors for each Rij component (default: [5, 35, 25, 25, 150, 20])
        Ny: Number of points for RANS profiles
        NyDNS: Number of points for DNS/HF profiles
        AR: Aspect ratio (default: 1.0)
    """
    # Get cell centers and coordinates
    current_points = internal_mesh.cell_centers().points
    ref_points = ref_internal_mesh.cell_centers().points

    # Normalize coordinates
    C_y_ref = ref_points[:, 1] / H
    C_z_ref = ref_points[:, 2] / H
    C_y_current = current_points[:, 1] / H
    C_z_current = current_points[:, 2] / H

    # Interpolate current model Rij to reference mesh coordinates for consistency
    # This ensures all three models are compared on the same grid
    logger.info("Interpolating current model Rij to reference mesh coordinates...")
    distances = cdist(ref_points, current_points)
    nearest_indices = np.argmin(distances, axis=1)
    rij_current_interp = rij_current[nearest_indices]

    # Compute bulk velocity from HF data
    if "U_HF" in ref_internal_mesh.point_data:
        u_hf_full = np.array(ref_internal_mesh.point_data["U_HF"])
        u_x_hf = u_hf_full[:, 0]
    elif "U_HF" in ref_internal_mesh.cell_data:
        u_hf_full = np.array(ref_internal_mesh.cell_data["U_HF"])
        u_x_hf = u_hf_full[:, 0]
    else:
        raise KeyError("U_HF not found for bulk velocity calculation")

    if "V" in ref_internal_mesh.cell_data:
        volume = np.array(ref_internal_mesh.cell_data["V"])
    elif "Volume" in ref_internal_mesh.cell_data:
        volume = np.array(ref_internal_mesh.cell_data["Volume"])
    else:
        volume = np.ones(ref_internal_mesh.n_cells)
        logger.warning("Volume data not found, using uniform weights for bulk velocity")

    if u_x_hf.ndim == 0 or len(u_x_hf) != len(volume):
        # Convert point data to cell data if needed
        if "U_HF" in ref_internal_mesh.point_data:
            u_x_hf = point_to_cell_data(ref_internal_mesh, u_x_hf)

    Ub = volumetric_average(u_x_hf, volume)

    # Create figure with 2x3 subplots (2 rows × 3 cols)
    fig, axes = plt.subplots(
        figsize=get_figure_size(3, 2),
        ncols=3,
        nrows=2,
        sharey=True,
        sharex=True,  # Match reference: sharex=True
        constrained_layout=False,
    )

    # Line styles and colours
    lw = 1
    ls = np.array(['--', ':'])

    # Rij component names and labels (dynamically generated from scale_rij)
    # Component indices: 0=R11, 1=R12, 2=R13, 3=R22, 4=R23, 5=R33
    rij_component_names = [
        [r"R_{11}", r"R_{12}", r"R_{13}"],
        [r"R_{22}", r"R_{23}", r"R_{33}"],
    ]
    rij_labels = []
    for row in range(2):
        row_labels = []
        for col in range(3):
            comp_idx = row * 3 + col
            scale_val = int(scale_rij[comp_idx])
            comp_name = rij_component_names[row][col]
            label = rf"${scale_val}\langle {comp_name} \rangle/u_{{b}}^{{2}} + y/H$"
            row_labels.append(label)
        rij_labels.append(row_labels)

    # Generate y profiles using beta distribution
    yProfile = beta.cdf(np.linspace(0, 1, Ny), 3, 1) - 1
    yProfileDNS = beta.cdf(np.linspace(0, 1, NyDNS), 3, 1) - 1

    # Verify data shapes and ranges
    logger.info(f"Rij data shapes - Current: {rij_current.shape}, kOmegaSST: {rij_kosst.shape}, HF: {rij_hf.shape}")
    logger.info(f"Rij data ranges - Current: [{rij_current.min():.4e}, {rij_current.max():.4e}], "
                f"kOmegaSST: [{rij_kosst.min():.4e}, {rij_kosst.max():.4e}], "
                f"HF: [{rij_hf.min():.4e}, {rij_hf.max():.4e}]")

    # Plot for each Rij component (6 components total)
    for j in range(6):
        row = j // 3
        col = j % 3
        ax = axes[row, col]

        # Get Rij component for each model
        rij_comp_kosst = rij_kosst[:, j]
        rij_comp_hf = rij_hf[:, j]
        # Use interpolated current model data on reference mesh for consistency
        rij_comp_current = rij_current_interp[:, j]
        
        # Log component statistics for debugging
        logger.debug(f"Rij[{j}] - Current: [{rij_comp_current.min():.4e}, {rij_comp_current.max():.4e}], "
                    f"kOmegaSST: [{rij_comp_kosst.min():.4e}, {rij_comp_kosst.max():.4e}], "
                    f"HF: [{rij_comp_hf.min():.4e}, {rij_comp_hf.max():.4e}]")

        # Plot profiles at each y/H location
        for i, profile_y in enumerate(profiles):
            # kOmegaSST profile
            interp_kosst = LinearNDInterpolator(
                list(zip(C_y_ref, C_z_ref)), rij_comp_kosst
            )
            xProfile_kosst = np.ones(Ny) * profile_y
            rijProfile_kosst = interp_kosst(xProfile_kosst, yProfile)
            valid_kosst = ~np.isnan(rijProfile_kosst)
            if valid_kosst.any():
                ax.plot(
                    rijProfile_kosst[valid_kosst] * scale_rij[j] / Ub ** 2 + profile_y,
                    yProfile[valid_kosst],
                    c=colours[0],  # navy blue (black equivalent)
                    lw=lw * 0.75,
                    ls='-',
                    zorder=2,
                )

            # Current model profile (use interpolated data on reference mesh)
            interp_current = LinearNDInterpolator(
                list(zip(C_y_ref, C_z_ref)), rij_comp_current
            )
            xProfile_current = np.ones(Ny) * profile_y
            rijProfile_current = interp_current(xProfile_current, yProfile)
            valid_current = ~np.isnan(rijProfile_current)
            if valid_current.any():
                ax.plot(
                    rijProfile_current[valid_current] * scale_rij[j] / Ub ** 2 + profile_y,
                    yProfile[valid_current],
                    c=colours[1],  # red
                    lw=lw,
                    ls=ls[0],  # dashed
                    zorder=3,
                )

            # High Fidelity profile
            interp_hf = LinearNDInterpolator(
                list(zip(C_y_ref, C_z_ref)), rij_comp_hf
            )
            xProfile_hf = np.ones(NyDNS) * profile_y
            rijProfile_hf = interp_hf(xProfile_hf, yProfileDNS)
            valid_hf = ~np.isnan(rijProfile_hf)
            if valid_hf.any():
                ax.plot(
                    rijProfile_hf[valid_hf] * scale_rij[j] / Ub ** 2 + profile_y,
                    yProfileDNS[valid_hf],
                    c='k',
                    lw=lw,
                    marker='o',
                    markevery=1,
                    markerfacecolor='None',
                    markeredgecolor='k',
                    markersize=2,
                    markeredgewidth=lw / 2,
                    ls='None',
                    zorder=1,
                )

        # Set labels
        ax.set_xlabel(rij_labels[row][col], fontsize=10)
        ax.set_box_aspect(1)

        # Set ticks and limits (matching reference code with user's 0-1 range requirement)
        # Reference code uses sharex=True and + profile_y for horizontal separation
        # User wants all plots to show 0-1 range (not 0-2 for first column)
        # Set xlim to [0, 1] for all plots
        ax.set_xlim([0, 1])
        xticks = np.arange(0, 1 + 0.25, 0.25)
        xticklabels = np.round(np.arange(0, 1 + 0.25, 0.25), 2)
        
        # Only set ticks on bottom row to avoid duplication (sharex=True)
        if row == 1:
            ax.set_xticks(xticks)
            ax.set_xticklabels([f'{t:.2f}' for t in xticklabels])

    # Y-axis settings
    deltay = 0.5
    yticks = np.arange(-AR, 0 + deltay, deltay)
    yticklabels = np.round(np.arange(0, AR + deltay, deltay), 2)
    for ax in axes.flat:
        ax.set_yticks(yticks)
        ax.set_yticklabels([f'{t:.2f}' for t in yticklabels])

    fig.supylabel(r"$z/H$", fontsize=10)

    # Add legend
    lines = [
        Line2D(
            [0], [0], marker='o', color='k', markersize=2, ls='None',
            markerfacecolor='None', markeredgewidth=lw / 2, markeredgecolor='k'
        ),
        Line2D([0], [0], color=colours[0], linewidth=lw * 0.75, linestyle='-'),
        Line2D([0], [0], color=colours[1], linewidth=lw, linestyle='--'),
    ]
    labels = ['High Fidelity', r'$k-\omega \ \text{SST}$', r'$k-\omega \ \text{SST-PDA}$']
    fig.legend(lines, labels, ncol=3, frameon=True, fancybox=True, bbox_to_anchor=(0.75, 0.96))

    # Save figure with improved spacing
    fig.tight_layout(pad=1.0, h_pad=0.3, w_pad=0.5)
    fig.savefig(output_path, dpi=1000, transparent=False, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Rij profiles plot saved to: {output_path}")


def create_comparison_contour_plot(
    sliced_mesh: pv.PolyData,
    normalised_velocity_kosst: np.ndarray,
    normalised_velocity_hf: np.ndarray,
    normalised_velocity_current: np.ndarray,
    normalised_error_velocity_kosst: np.ndarray,
    normalised_error_velocity_hf: np.ndarray,
    normalised_error_velocity_current: np.ndarray,
    normalised_error_vorticity_kosst: np.ndarray,
    normalised_error_vorticity_hf: np.ndarray,
    normalised_error_vorticity_current: np.ndarray,
    u_y_kosst: Optional[np.ndarray] = None,
    u_z_kosst: Optional[np.ndarray] = None,
    u_y_hf: Optional[np.ndarray] = None,
    u_z_hf: Optional[np.ndarray] = None,
    u_y_current: Optional[np.ndarray] = None,
    u_z_current: Optional[np.ndarray] = None,
) -> Tuple[plt.Figure, np.ndarray]:
    """
    Create comparison contour plots with streamlines for kOmegaSST, HF, and kOmegaSSTPDA.

    Args:
        sliced_mesh: Sliced PyVista mesh
        normalised_velocity_*: Normalised velocity fields for each model
        normalised_error_velocity_*: Normalised velocity error fields for each model
        normalised_error_vorticity_*: Normalised vorticity error fields for each model
        u_y_*, u_z_*: Optional velocity components for streamlines

    Returns:
        Matplotlib figure and axes array
    """
    # Get coordinates from sliced mesh
    points = sliced_mesh.points
    y_coords = points[:, 1]  # y coordinates
    z_coords = points[:, 2]  # z coordinates

    # Normalize coordinates to match postProcessMulti.py format
    y_norm = -y_coords / H
    z_norm = z_coords / H

    # Create regular grid for contour plotting
    resolution = 100
    y_grid = np.linspace(y_norm.min(), y_norm.max(), resolution)
    z_grid = np.linspace(z_norm.min(), z_norm.max(), resolution)
    yi, zi = np.meshgrid(y_grid, z_grid)

    # Interpolate fields onto regular grid
    points_2d = np.column_stack((y_norm, z_norm))
    
    # Interpolate all fields
    fields_to_interpolate = [
        (normalised_velocity_kosst, "normalised_velocity_kosst"),
        (normalised_velocity_hf, "normalised_velocity_hf"),
        (normalised_velocity_current, "normalised_velocity_current"),
        (normalised_error_velocity_kosst, "normalised_error_velocity_kosst"),
        (normalised_error_velocity_hf, "normalised_error_velocity_hf"),
        (normalised_error_velocity_current, "normalised_error_velocity_current"),
        (normalised_error_vorticity_kosst, "normalised_error_vorticity_kosst"),
        (normalised_error_vorticity_hf, "normalised_error_vorticity_hf"),
        (normalised_error_vorticity_current, "normalised_error_vorticity_current"),
    ]
    
    interpolated_fields = {}
    for field_data, field_name in fields_to_interpolate:
        interpolated_fields[field_name] = griddata(
            points_2d, field_data, (yi, zi), method='linear'
        )

    # Create figure with 3 rows (field types) and 3 columns (models) (3 rows × 3 cols)
    # Rows: velocity, velocity error, vorticity error
    # Columns: kOmegaSST, kOmegaSSTPDA, HF
    fig, ax = plt.subplots(
        figsize=get_figure_size(3, 3),
        ncols=3,
        nrows=3,
        sharey='row',
        sharex='col',
    )

    # Color limits and normalization
    vmin = -1
    vmax = 1
    delta = 0.5

    norm1 = matplotlib.colors.Normalize(vmin=0, vmax=1.5)
    norm2 = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    norm3 = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    # Model configuration: (name, column_index, velocity_field, error_velocity_field, error_vorticity_field, u_y, u_z)
    models = [
        (r"$k-\omega \ \text{SST}$", 0, "normalised_velocity_kosst", "normalised_error_velocity_kosst", 
         "normalised_error_vorticity_kosst", u_y_kosst, u_z_kosst),
        (r"$k-\omega \ \text{SST-PDA}$", 1, "normalised_velocity_current", "normalised_error_velocity_current",
         "normalised_error_vorticity_current", u_y_current, u_z_current),
        ("High Fidelity", 2, "normalised_velocity_hf", "normalised_error_velocity_hf",
         "normalised_error_vorticity_hf", u_y_hf, u_z_hf),
    ]

    # Pre-compute surface velocities and find HF maximum for normalization
    surface_velocity_plots = {}
    hf_surface_velocity_max = None
    
    for model_name, col_idx, field_vel, field_err_vel, field_err_vort, u_y, u_z in models:
        if u_y is not None and u_z is not None:
            surface_velocity = np.sqrt(u_y ** 2 + u_z ** 2)
            surface_velocity_plot = griddata(
                points_2d, surface_velocity, (yi, zi), method='nearest'
            )
            surface_velocity_plot = surface_velocity_plot.reshape(
                resolution, resolution
            )
            surface_velocity_plots[col_idx] = {
                'u_y': u_y,
                'u_z': u_z,
                'surface_velocity': surface_velocity_plot,
            }
            
            # Find HF maximum (HF is in column 2)
            if col_idx == 2:
                hf_surface_velocity_max = np.nanmax(surface_velocity_plot)
    
    # If HF data is not available, fall back to using the maximum of all models
    if hf_surface_velocity_max is None or hf_surface_velocity_max == 0:
        all_maxima = [
            np.nanmax(data['surface_velocity'])
            for data in surface_velocity_plots.values()
        ]
        if all_maxima:
            hf_surface_velocity_max = max(all_maxima)
        else:
            hf_surface_velocity_max = 1.0  # Fallback to prevent division by zero

    # Row 0: Normalized velocity with streamlines
    row_idx = 0
    cs = None
    for model_name, col_idx, field_vel, field_err_vel, field_err_vort, u_y, u_z in models:
        cs = ax[row_idx, col_idx].contourf(
            yi,
            zi,
            interpolated_fields[field_vel],
            levels=np.linspace(0, 1.5, 35),
            cmap='turbo',
            vmin=0,
            vmax=1.5,
            norm=norm1,
            extend='both',
        )

        # Add streamlines if available
        if col_idx in surface_velocity_plots:
            data = surface_velocity_plots[col_idx]
            gu = griddata(points_2d, -data['u_y'], (yi, zi), method='linear')
            gv = griddata(points_2d, data['u_z'], (yi, zi), method='linear')
            
            gu = gu.reshape(resolution, resolution)
            gv = gv.reshape(resolution, resolution)
            surface_velocity_plot = data['surface_velocity']

            # Normalize using HF maximum for consistent line thickness across all models
            lw = surface_velocity_plot / hf_surface_velocity_max
            ax[row_idx, col_idx].streamplot(
                yi,
                zi,
                gu,
                gv,
                density=2,
                linewidth=lw,
                color='k',
                arrowsize=0.5,
                arrowstyle='fancy',
            )

        # Set aspect ratio
        ax[row_idx, col_idx].set_aspect('equal', 'box')
        
        # Set column title (model name)
        ax[row_idx, col_idx].set_title(model_name, fontsize=10)

    # Row 1: Velocity error
    row_idx = 1
    csUst = None
    for model_name, col_idx, field_vel, field_err_vel, field_err_vort, u_y, u_z in models:
        csUst = ax[row_idx, col_idx].contourf(
            yi,
            zi,
            interpolated_fields[field_err_vel],
            levels=np.linspace(vmin, vmax, 35),
            cmap='RdBu_r',
            vmin=vmin,
            vmax=vmax,
            norm=norm2,
            extend='both',
        )

        # Set aspect ratio
        ax[row_idx, col_idx].set_aspect('equal', 'box')

    # Row 2: Vorticity error
    row_idx = 2
    csVort = None
    for model_name, col_idx, field_vel, field_err_vel, field_err_vort, u_y, u_z in models:
        csVort = ax[row_idx, col_idx].contourf(
            yi,
            zi,
            interpolated_fields[field_err_vort],
            levels=np.linspace(vmin, vmax, 35),
            cmap='PuOr_r',
            vmin=vmin,
            vmax=vmax,
            norm=norm3,
            extend='both',
        )

        # Set aspect ratio
        ax[row_idx, col_idx].set_aspect('equal', 'box')

    # Add colorbars (one per row, vertical, outside on the right side)
    # Anchor to all subplots in each row to match the row height exactly
    # Row 0: Velocity
    cbar1 = fig.colorbar(
        cs, cmap='turbo', pad=0.15, orientation='vertical',
        norm=norm1, ax=ax[0, :], aspect=20, shrink=1.0
    )
    # Row 1: Velocity error
    cbar2 = fig.colorbar(
        csUst, cmap='RdBu_r', pad=0.15, orientation='vertical',
        norm=norm2, ax=ax[1, :], aspect=20, shrink=1.0
    )
    # Row 2: Vorticity error
    cbar3 = fig.colorbar(
        csVort, cmap='PuOr_r', pad=0.15, orientation='vertical',
        norm=norm3, ax=ax[2, :], aspect=20, shrink=1.0
    )

    # Set colorbar labels
    cbar1.set_ticks([0, 0.5, 1.0, 1.5])
    cbar1.set_ticklabels([0, 0.5, 1.0, 1.5])
    cbar1.set_label(r'$\frac{\langle u \rangle}{U_b}$', rotation=90, labelpad=10)

    ticks = np.arange(vmin, vmax + delta, delta)
    ticklabels = np.round(np.arange(vmin, vmax + delta, delta), 2)

    cbar2.set_ticks(ticks)
    cbar2.set_ticklabels(ticklabels)
    cbar2.set_label(
        r'$\frac{\langle u \rangle - \langle u \rangle_{LES}}{\epsilon_{u}}$',
        rotation=90,
        labelpad=10,
    )

    cbar3.set_ticks(ticks)
    cbar3.set_ticklabels(ticklabels)
    cbar3.set_label(
        r'$\frac{\langle \omega_{x} \rangle - \langle \omega_{x} \rangle_{LES}}{\epsilon_{\omega_{x}}}$',
        rotation=90,
        labelpad=10,
    )

    # Set row titles (field types)
    ax[0, 0].set_ylabel(r'$\frac{\langle u \rangle}{U_b}$\n$z/h$', fontsize=10)
    ax[1, 0].set_ylabel(
        r'$\frac{\langle u \rangle - \langle u \rangle_{LES}}{\epsilon_{u}}$\n$z/h$',
        fontsize=10,
    )
    ax[2, 0].set_ylabel(
        r'$\frac{\langle \omega_{x} \rangle - \langle \omega_{x} \rangle_{LES}}{\epsilon_{\omega_{x}}}$\n$z/h$',
        fontsize=10,
    )

    # Set axis ticks and labels
    y_ticks = [y_norm.min(), (y_norm.min() + y_norm.max()) / 2, y_norm.max()]
    z_ticks = [z_norm.min(), (z_norm.min() + z_norm.max()) / 2, z_norm.max()]

    for col in range(3):
        ax[2, col].set_xticks(y_ticks)
        ax[2, col].set_xticklabels(['0', '0.5', '1'])
        ax[2, col].set_xlabel(r"$y/h$")

    for row in range(3):
        ax[row, 0].set_yticks(z_ticks)
        ax[row, 0].set_yticklabels(['0', '0.5', '1'])

    # Apply tight layout with space reserved on the right for colorbars
    # rect parameter: [left, bottom, right, top] in figure coordinates
    fig.tight_layout(pad=1.0, h_pad=0.3, w_pad=0.5, rect=[0, 0, 0.9, 1])

    return fig, ax


def create_matplotlib_contour_plot(
    sliced_mesh: pv.PolyData,
    normalised_velocity: np.ndarray,
    normalised_error_velocity: np.ndarray,
    normalised_error_vorticity: np.ndarray,
    norm_error_combined: float,
    norm_error_velocity: float,
    norm_error_vorticity: float,
    u_y: Optional[np.ndarray] = None,
    u_z: Optional[np.ndarray] = None,
) -> Tuple[plt.Figure, np.ndarray]:
    """
    Create matplotlib contour plots similar to postProcessMulti.py.

    Args:
        sliced_mesh: Sliced PyVista mesh
        normalised_velocity: Normalised velocity field
        normalised_error_velocity: Normalised velocity error field
        normalised_error_vorticity: Normalised vorticity error field
        norm_error_combined: Combined error metric (J)
        norm_error_velocity: Velocity error metric (j1)
        norm_error_vorticity: Vorticity error metric (j2)
        u_y: Optional y-component of velocity for streamlines
        u_z: Optional z-component of velocity for streamlines

    Returns:
        Matplotlib figure and axes array
    """
    # Get coordinates from sliced mesh
    points = sliced_mesh.points
    y_coords = points[:, 1]  # y coordinates
    z_coords = points[:, 2]  # z coordinates

    # Normalize coordinates to match postProcessMulti.py format
    y_norm = -y_coords / H
    z_norm = z_coords / H

    # Create regular grid for contour plotting
    resolution = 100
    y_grid = np.linspace(y_norm.min(), y_norm.max(), resolution)
    z_grid = np.linspace(z_norm.min(), z_norm.max(), resolution)
    yi, zi = np.meshgrid(y_grid, z_grid)

    # Interpolate fields onto regular grid
    points_2d = np.column_stack((y_norm, z_norm))
    normalised_velocity_grid = griddata(
        points_2d, normalised_velocity, (yi, zi), method='linear'
    )
    normalised_error_velocity_grid = griddata(
        points_2d, normalised_error_velocity, (yi, zi), method='linear'
    )
    normalised_error_vorticity_grid = griddata(
        points_2d, normalised_error_vorticity, (yi, zi), method='linear'
    )

    # Create figure with 3 subplots (1 row × 3 cols)
    fig, ax = plt.subplots(
        figsize=get_figure_size(3, 1), ncols=3, nrows=1, sharey=True
    )

    # Color limits and normalization
    vmin = -1
    vmax = 1
    delta = 0.5

    norm1 = matplotlib.colors.Normalize(vmin=0, vmax=1.5)
    norm2 = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    norm3 = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    # Create contour plots
    cs = ax[0].contourf(
        yi,
        zi,
        normalised_velocity_grid,
        levels=np.linspace(0, 1.5, 35),
        cmap='turbo',
        vmin=0,
        vmax=1.5,
        norm=norm1,
        extend='both',
    )
    csUst = ax[1].contourf(
        yi,
        zi,
        normalised_error_velocity_grid,
        levels=np.linspace(vmin, vmax, 35),
        cmap='RdBu_r',
        vmin=vmin,
        vmax=vmax,
        norm=norm2,
        extend='both',
    )
    csVort = ax[2].contourf(
        yi,
        zi,
        normalised_error_vorticity_grid,
        levels=np.linspace(vmin, vmax, 35),
        cmap='PuOr_r',
        vmin=vmin,
        vmax=vmax,
        norm=norm3,
        extend='both',
    )

    # Add colorbars
    cbar1 = fig.colorbar(
        cs, cmap='turbo', pad=0.15, orientation='horizontal', norm=norm1, ax=ax[0]
    )
    cbar2 = fig.colorbar(
        csUst, cmap='RdBu_r', pad=0.15, orientation='horizontal', norm=norm2, ax=ax[1]
    )
    cbar3 = fig.colorbar(
        csVort, cmap='PuOr_r', pad=0.15, orientation='horizontal', norm=norm3, ax=ax[2]
    )

    # Set colorbar ticks and labels
    ticks = np.arange(vmin, vmax + delta, delta)
    ticklabels = np.round(np.arange(vmin, vmax + delta, delta), 2)

    cbar1.set_ticks([0, 0.5, 1.0, 1.5])
    cbar1.set_ticklabels([0, 0.5, 1.0, 1.5])
    cbar1.set_label(r'$\frac{\langle u \rangle}{U_b}$', rotation=0, loc='right')

    cbar2.set_ticks(ticks)
    cbar2.set_ticklabels(ticklabels)
    cbar2.set_label(
        r'$\frac{\langle u \rangle - \langle u \rangle_{LES}}{\epsilon_{u}}$',
        rotation=0,
        loc='right',
    )

    cbar3.set_ticks(ticks)
    cbar3.set_ticklabels(ticklabels)
    cbar3.set_label(
        r'$\frac{\langle \omega_{x} \rangle - \langle \omega_{x} \rangle_{LES}}{\epsilon_{\omega_{x}}}$',
        rotation=0,
        loc='right',
    )

    # Add streamlines to first plot if velocity components are available
    if u_y is not None and u_z is not None:
        surface_velocity = np.sqrt(u_y ** 2 + u_z ** 2)
        gu = griddata(points_2d, -u_y, (yi, zi), method='linear')
        gv = griddata(points_2d, u_z, (yi, zi), method='linear')
        surface_velocity_plot = griddata(
            points_2d, surface_velocity, (yi, zi), method='nearest'
        )

        gu = gu.reshape(resolution, resolution)
        gv = gv.reshape(resolution, resolution)
        surface_velocity_plot = surface_velocity_plot.reshape(
            resolution, resolution
        )

        lw = surface_velocity_plot / surface_velocity_plot.max()
        ax[0].streamplot(
            yi,
            zi,
            gu,
            gv,
            density=2,
            linewidth=lw,
            color='k',
            arrowsize=0.5,
            arrowstyle='fancy',
        )

    # Set aspect ratio
    ax[0].set_aspect('equal', 'box')
    ax[1].set_aspect('equal', 'box')
    ax[2].set_aspect('equal', 'box')

    # Set axis ticks and labels
    y_ticks = [y_norm.min(), (y_norm.min() + y_norm.max()) / 2, y_norm.max()]
    z_ticks = [z_norm.min(), (z_norm.min() + z_norm.max()) / 2, z_norm.max()]
    
    ax[0].set_xticks(y_ticks)
    ax[0].set_xticklabels(['0', '0.5', '1'])
    ax[1].set_xticks(y_ticks)
    ax[1].set_xticklabels(['0', '0.5', '1'])
    ax[2].set_xticks(y_ticks)
    ax[2].set_xticklabels(['0', '0.5', '1'])

    ax[0].set_yticks(z_ticks)
    ax[0].set_yticklabels(['0', '0.5', '1'])

    ax[0].set_xlabel(r"$y/h$")
    ax[1].set_xlabel(r"$y/h$")
    ax[2].set_xlabel(r"$y/h$")
    ax[0].set_ylabel(r"$z/h$")

    # Set titles
    ax[0].set_title(
        r'$J = ' + str(np.round(norm_error_combined, 3)) + '$',
        loc='left',
        fontsize=7,
    )
    ax[1].set_title(
        r'$j_{1} = ' + str(np.round(norm_error_velocity, 3)) + '$',
        loc='left',
        fontsize=7,
    )
    ax[2].set_title(
        r'$j_{2} = ' + str(np.round(norm_error_vorticity, 3)) + '$',
        loc='left',
        fontsize=7,
    )

    # Apply tight layout with improved spacing
    fig.tight_layout(pad=1.0, h_pad=0.5, w_pad=0.5)

    return fig, ax


@click.command()
@click.option(
    '--plot-multi',
    is_flag=True,
    default=True,
    help='Create multi plot (velocity, velocity error, vorticity error) [default: True]'
)
@click.option(
    '--plot-comparison',
    is_flag=True,
    default=False,
    help='Create comparison contour plot (kOmegaSST, HF, kOmegaSSTPDA) [default: True]'
)
@click.option(
    '--plot-barycentric',
    is_flag=True,
    default=True,
    help='Create barycentric triangle plots [default: True]'
)
@click.option(
    '--plot-anisotropy',
    is_flag=True,
    default=False,
    help='Create anisotropy contour plot [default: True]'
)
@click.option(
    '--plot-velocity-profiles',
    is_flag=True,
    default=True,
    help='Create velocity profiles plot [default: True]'
)
@click.option(
    '--plot-rij-profiles',
    is_flag=True,
    default=False,
    help='Create Rij component profiles plot [default: True]'
)
@click.option(
    '--plot-diagonal-profile',
    is_flag=True,
    default=True,
    help='Create diagonal velocity profile plot (u₂/u_τ vs y/H along y=z) [default: True]'
)
@click.option(
    '--max-iterations',
    type=int,
    default=4000,
    help='Maximum number of iterations. If latest time folder equals this value, objectives are set to fail_value and marked as not converged [default: 10000]'
)
@click.option(
    '--fail-value',
    type=float,
    default=2.0,
    help='Objective value to set when max iterations are reached or simulation fails (only used if --use-fail-value is set) [default: 2.0]'
)
@click.option(
    '--use-fail-value',
    is_flag=True,
    default=False,
    help='If set, override objectives with fail_value when max iterations reached or simulation fails. If not set, use calculated objectives but mark as not converged [default: False]'
)
def main(
    plot_multi: bool,
    plot_comparison: bool,
    plot_barycentric: bool,
    plot_anisotropy: bool,
    plot_velocity_profiles: bool,
    plot_rij_profiles: bool,
    plot_diagonal_profile: bool,
    max_iterations: int,
    fail_value: float,
    use_fail_value: bool,
) -> None:
    """
    Main execution function.
    
    The objective function calculation (RMSE) is always performed.
    Individual plots can be enabled/disabled using command-line flags.
    """
    work_dir = Path.cwd()
    logger.info(f"Working directory: {work_dir}")
    logger.info(f"Plot options - Multi: {plot_multi}, Comparison: {plot_comparison}, "
                f"Barycentric: {plot_barycentric}, Anisotropy: {plot_anisotropy}, "
                f"Velocity profiles: {plot_velocity_profiles}, Rij profiles: {plot_rij_profiles}, "
                f"Diagonal profile: {plot_diagonal_profile}")

    # Load current simulation mesh using PyVista
    latest_folder = find_latest_time_folder(work_dir)

    if latest_folder is None:
        logger.warning("No time folders found, using default error values")
        norm_error_velocity = 1000.0
        norm_error_vorticity = 1000.0
        norm_error_combined = 1000.0
        converged = False
        logger.warning("Cannot create visualisations without simulation data")
        # Save results with not converged flag
        post_processing_dir = work_dir / "postProcessing"
        post_processing_dir.mkdir(exist_ok=True)
        results_file = post_processing_dir / "objective_fun_results.json"
        save_objective_function_results(
            results_file,
            norm_error_velocity,
            norm_error_vorticity,
            norm_error_combined,
            converged=converged,
        )
        return

    logger.info(f"Processing time folder: {latest_folder}")
    
    # Check if simulation reached maximum iterations
    converged = True
    max_iterations_reached = (latest_folder == str(max_iterations))
    
    if max_iterations_reached:
        converged = False
        if use_fail_value:
            logger.warning(
                f"Simulation reached maximum iterations ({max_iterations}). "
                f"Setting all objectives to {fail_value} and marking as not converged."
            )
        else:
            logger.warning(
                f"Simulation reached maximum iterations ({max_iterations}). "
                f"Using calculated objectives but marking as not converged."
            )

    try:
        # Load current simulation mesh
        reader, mesh, internal_mesh = load_openfoam_mesh(work_dir)

        # Get cell volumes from PyVista mesh
        volume = get_cell_volumes(internal_mesh)
        logger.info(f"Retrieved {len(volume)} cell volumes from mesh")

        # Slice at x midpoint (automatically calculated)
        sliced_mesh = slice_mesh_at_x_plane(internal_mesh)

        # Get velocity and vorticity fields from sliced mesh
        if "U" in sliced_mesh.point_data:
            u_current = np.array(sliced_mesh.point_data["U"])
            u_x_current = u_current[:, 0]
        else:
            raise KeyError("Velocity field 'U' not found in mesh point_data")

        if "vorticity" in sliced_mesh.point_data:
            vort_current = np.array(sliced_mesh.point_data["vorticity"])
            vort_x_current = vort_current[:, 0]
        else:
            raise KeyError(
                "Vorticity field 'vorticity' not found in mesh point_data"
            )

        # Get full 3D fields for error computation
        if "U" in internal_mesh.point_data:
            u_full = np.array(internal_mesh.point_data["U"])
            u_x_full = u_full[:, 0]
        else:
            raise KeyError("Velocity field 'U' not found in internal mesh")

        if "vorticity" in internal_mesh.point_data:
            vort_full = np.array(internal_mesh.point_data["vorticity"])
            vort_x_full = vort_full[:, 0]
        else:
            raise KeyError("Vorticity field 'vorticity' not found in internal mesh")

        # Convert point data to cell data for volumetric averaging
        u_x_full_cell = point_to_cell_data(internal_mesh, u_x_full)
        vort_x_full_cell = point_to_cell_data(internal_mesh, vort_x_full)

        # Load reference data from the 0 time folder using PyVista
        # The reference fields (U_HF, U_kOmegaSST, vort_HF, vort_kOmegaSST) are in the same case
        logger.info("Loading reference data from time 0 using PyVista...")

        # Load mesh at time 0 to get reference fields
        ref_reader, ref_mesh, ref_internal_mesh = load_openfoam_mesh(
            work_dir, time_value="0"
        )

        # Get reference fields from the mesh
        # These should be available in point_data or cell_data
        logger.info("Extracting reference fields from mesh...")

        # Get LES (HF) fields
        try:
            if "U_HF" in ref_internal_mesh.point_data:
                u_les_full = np.array(ref_internal_mesh.point_data["U_HF"])
                u_x_les_full = u_les_full[:, 0]
            elif "U_HF" in ref_internal_mesh.cell_data:
                u_les_full = np.array(ref_internal_mesh.cell_data["U_HF"])
                u_x_les_full = u_les_full[:, 0]
            else:
                raise KeyError("U_HF field not found in reference mesh")

            if "vort_HF" in ref_internal_mesh.point_data:
                vort_les_full = np.array(ref_internal_mesh.point_data["vort_HF"])
                vort_x_les_full = vort_les_full[:, 0]
            elif "vorticity_HF" in ref_internal_mesh.point_data:
                vort_les_full = np.array(
                    ref_internal_mesh.point_data["vorticity_HF"]
                )
                vort_x_les_full = vort_les_full[:, 0]
            elif "vort_HF" in ref_internal_mesh.cell_data:
                vort_les_full = np.array(ref_internal_mesh.cell_data["vort_HF"])
                vort_x_les_full = vort_les_full[:, 0]
            else:
                raise KeyError("vort_HF or vorticity_HF field not found in reference mesh")

            # Get kOmegaSST fields
            if "U_KOSST" in ref_internal_mesh.point_data:
                u_kosst_full = np.array(ref_internal_mesh.point_data["U_KOSST"])
                u_x_kosst_full = u_kosst_full[:, 0]
            elif "U_KOSST" in ref_internal_mesh.cell_data:
                u_kosst_full = np.array(ref_internal_mesh.cell_data["U_KOSST"])
                u_x_kosst_full = u_kosst_full[:, 0]
            else:
                raise KeyError("U_KOSST field not found in reference mesh")

            if "vort_KOSST" in ref_internal_mesh.point_data:
                vort_kosst_full = np.array(
                    ref_internal_mesh.point_data["vort_KOSST"]
                )
                vort_x_kosst_full = vort_kosst_full[:, 0]
            elif "vorticity_KOSST" in ref_internal_mesh.point_data:
                vort_kosst_full = np.array(
                    ref_internal_mesh.point_data["vorticity_KOSST"]
                )
                vort_x_kosst_full = vort_kosst_full[:, 0]
            elif "vort_KOSST" in ref_internal_mesh.cell_data:
                vort_kosst_full = np.array(
                    ref_internal_mesh.cell_data["vort_KOSST"]
                )
                vort_x_kosst_full = vort_kosst_full[:, 0]
            else:
                raise KeyError(
                    "vort_KOSST or vorticity_KOSST field not found in reference mesh"
                )

        except KeyError as e:
            logger.error(f"Reference field not found: {e}")
            raise

        # Interpolate reference fields to current mesh cell centers for error computation
        logger.info("Interpolating reference fields to current mesh...")
        # Get cell centers
        current_cell_centers = internal_mesh.cell_centers().points
        ref_cell_centers = ref_internal_mesh.cell_centers().points

        # Use nearest neighbor interpolation
        distances = cdist(current_cell_centers, ref_cell_centers)
        nearest_indices = np.argmin(distances, axis=1)

        # Convert reference point data to cell data if needed
        if "U_HF" in ref_internal_mesh.point_data:
            u_x_les_ref_cell = point_to_cell_data(
                ref_internal_mesh, u_x_les_full
            )
        else:
            u_x_les_ref_cell = u_x_les_full

        if "vort_HF" in ref_internal_mesh.point_data or "vorticity_HF" in ref_internal_mesh.point_data:
            vort_x_les_ref_cell = point_to_cell_data(
                ref_internal_mesh, vort_x_les_full
            )
        else:
            vort_x_les_ref_cell = vort_x_les_full

        if "U_KOSST" in ref_internal_mesh.point_data:
            u_x_kosst_ref_cell = point_to_cell_data(
                ref_internal_mesh, u_x_kosst_full
            )
        else:
            u_x_kosst_ref_cell = u_x_kosst_full

        if "vort_KOSST" in ref_internal_mesh.point_data or "vorticity_KOSST" in ref_internal_mesh.point_data:
            vort_x_kosst_ref_cell = point_to_cell_data(
                ref_internal_mesh, vort_x_kosst_full
            )
        else:
            vort_x_kosst_ref_cell = vort_x_kosst_full

        # Interpolate to current mesh
        u_x_les_cell = u_x_les_ref_cell[nearest_indices]
        vort_x_les_cell = vort_x_les_ref_cell[nearest_indices]
        u_x_kosst_cell = u_x_kosst_ref_cell[nearest_indices]
        vort_x_kosst_cell = vort_x_kosst_ref_cell[nearest_indices]

        # Interpolate to slice for visualisation
        slice_points = sliced_mesh.points
        if "U_HF" in ref_internal_mesh.point_data:
            ref_points = ref_internal_mesh.points
        else:
            ref_points = ref_internal_mesh.cell_centers().points

        distances_slice = cdist(slice_points, ref_points)
        nearest_indices_slice = np.argmin(distances_slice, axis=1)

        u_x_les_slice = u_x_les_full[nearest_indices_slice]
        vort_x_les_slice = vort_x_les_full[nearest_indices_slice]

        # Compute baseline errors
        diff_vort_kosst_cell = np.abs(vort_x_kosst_cell - vort_x_les_cell)
        avg_diff_vort_kosst = volumetric_average(
            diff_vort_kosst_cell, volume
        )

        diff_u_kosst_cell = np.abs(u_x_kosst_cell - u_x_les_cell)
        avg_diff_u_kosst = volumetric_average(diff_u_kosst_cell, volume)

        # Compute error metrics (normalised vs kOmegaSST baseline)
        norm_error_velocity = compute_error_metrics(
            u_x_full_cell, u_x_les_cell, volume, avg_diff_u_kosst
        )

        norm_error_vorticity = compute_error_metrics(
            vort_x_full_cell, vort_x_les_cell, volume, avg_diff_vort_kosst
        )

        # Compute additional error metrics vs HF
        # Current model vs HF
        rmse_velocity = compute_rmse(u_x_full_cell, u_x_les_cell, volume)
        rmse_vorticity = compute_rmse(vort_x_full_cell, vort_x_les_cell, volume)
        mae_velocity = compute_mae(u_x_full_cell, u_x_les_cell, volume)
        mae_vorticity = compute_mae(vort_x_full_cell, vort_x_les_cell, volume)
        relative_rmse_velocity = compute_relative_rmse(u_x_full_cell, u_x_les_cell, volume)
        relative_rmse_vorticity = compute_relative_rmse(vort_x_full_cell, vort_x_les_cell, volume)
        relative_mae_velocity = compute_relative_mae(u_x_full_cell, u_x_les_cell, volume)
        relative_mae_vorticity = compute_relative_mae(vort_x_full_cell, vort_x_les_cell, volume)

        # kOmegaSST vs HF (baseline for comparison)
        rmse_velocity_kosst = compute_rmse(u_x_kosst_cell, u_x_les_cell, volume)
        rmse_vorticity_kosst = compute_rmse(vort_x_kosst_cell, vort_x_les_cell, volume)
        mae_velocity_kosst = compute_mae(u_x_kosst_cell, u_x_les_cell, volume)
        mae_vorticity_kosst = compute_mae(vort_x_kosst_cell, vort_x_les_cell, volume)

    except (FileNotFoundError, KeyError, RuntimeError) as e:
        logger.error(f"Error processing simulation data: {e}")
        converged = False
        if use_fail_value:
            norm_error_velocity = fail_value
            norm_error_vorticity = fail_value
            norm_error_combined = fail_value
            logger.info(f"Objectives set to {fail_value} due to processing error")
        else:
            norm_error_velocity = 1000.0
            norm_error_vorticity = 1000.0
            norm_error_combined = 1000.0
            logger.info("Using default error values (fail-value override disabled)")
        sliced_mesh = None
        u_x_current = None
        vort_x_current = None
        u_x_les_slice = None
        vort_x_les_slice = None
        avg_diff_vort_kosst = 1.0
        # Set additional error metrics to None if computation failed
        rmse_velocity = None
        rmse_vorticity = None
        mae_velocity = None
        mae_vorticity = None
        relative_rmse_velocity = None
        relative_rmse_vorticity = None
        relative_mae_velocity = None
        relative_mae_vorticity = None
        rmse_velocity_kosst = None
        rmse_vorticity_kosst = None
        mae_velocity_kosst = None
        mae_vorticity_kosst = None

    # Combined error metric
    norm_error_combined = (norm_error_velocity + norm_error_vorticity) / 2.0
    
    # Override objectives if max iterations reached and flag is set
    if max_iterations_reached and use_fail_value:
        norm_error_velocity = fail_value
        norm_error_vorticity = fail_value
        norm_error_combined = fail_value
        logger.info(f"Objectives overridden to {fail_value} due to max iterations reached")
    elif max_iterations_reached:
        logger.info("Using calculated objectives (fail-value override disabled)")
    else:
        logger.info(f"Error metrics - j1: {norm_error_velocity:.6f}, j2: {norm_error_vorticity:.6f}, J: {norm_error_combined:.6f}")

    # Print objective values to screen
    logger.info("=" * 60)
    logger.info("Objective Function Results:")
    logger.info(f"  j1 (velocity error):     {norm_error_velocity:.6f}")
    logger.info(f"  j2 (vorticity error):    {norm_error_vorticity:.6f}")
    logger.info(f"  J (combined):            {norm_error_combined:.6f}")
    logger.info(f"  Converged:               {converged if 'converged' in locals() else True}")
    logger.info("=" * 60)

    # Create postProcessing directory if it doesn't exist
    post_processing_dir = work_dir / "postProcessing"
    post_processing_dir.mkdir(exist_ok=True)
    logger.info(f"Post-processing directory: {post_processing_dir}")

    # Save objective function results
    results_file = post_processing_dir / "objective_fun_results.json"
    save_objective_function_results(
        results_file,
        norm_error_velocity,
        norm_error_vorticity,
        norm_error_combined,
        rmse_velocity=rmse_velocity if 'rmse_velocity' in locals() else None,
        rmse_vorticity=rmse_vorticity if 'rmse_vorticity' in locals() else None,
        mae_velocity=mae_velocity if 'mae_velocity' in locals() else None,
        mae_vorticity=mae_vorticity if 'mae_vorticity' in locals() else None,
        relative_rmse_velocity=relative_rmse_velocity if 'relative_rmse_velocity' in locals() else None,
        relative_rmse_vorticity=relative_rmse_vorticity if 'relative_rmse_vorticity' in locals() else None,
        relative_mae_velocity=relative_mae_velocity if 'relative_mae_velocity' in locals() else None,
        relative_mae_vorticity=relative_mae_vorticity if 'relative_mae_vorticity' in locals() else None,
        rmse_velocity_kosst=rmse_velocity_kosst if 'rmse_velocity_kosst' in locals() else None,
        rmse_vorticity_kosst=rmse_vorticity_kosst if 'rmse_vorticity_kosst' in locals() else None,
        mae_velocity_kosst=mae_velocity_kosst if 'mae_velocity_kosst' in locals() else None,
        mae_vorticity_kosst=mae_vorticity_kosst if 'mae_vorticity_kosst' in locals() else None,
        converged=converged if 'converged' in locals() else True,
    )

    # Create visualisations if data is available
    if (
        sliced_mesh is not None
        and u_x_current is not None
        and u_x_les_slice is not None
    ):
        logger.info("Creating matplotlib visualisations...")
        
        # Create multi plot if requested
        if plot_multi:
            logger.info("Creating multi plot...")
            # Get kOmegaSST data for normalisation (interpolate to slice)
            # Load reference mesh again for slice interpolation
            ref_reader, ref_mesh, ref_internal_mesh = load_openfoam_mesh(
                work_dir, time_value="0"
            )

            # Get U_kOmegaSST for slice
            if "U_KOSST" in ref_internal_mesh.point_data:
                u_kosst_slice_full = np.array(ref_internal_mesh.point_data["U_KOSST"])
                u_x_kosst_slice_full = u_kosst_slice_full[:, 0]
                ref_points_slice = ref_internal_mesh.points
            elif "U_KOSST" in ref_internal_mesh.cell_data:
                u_kosst_slice_full = np.array(ref_internal_mesh.cell_data["U_KOSST"])
                u_x_kosst_slice_full = u_kosst_slice_full[:, 0]
                ref_points_slice = ref_internal_mesh.cell_centers().points
            else:
                raise KeyError("U_KOSST not found for slice interpolation")

            # Get slice points from sliced mesh
            slice_points_viz = sliced_mesh.points
            distances_slice_kosst = cdist(slice_points_viz, ref_points_slice)
            nearest_indices_slice_kosst = np.argmin(distances_slice_kosst, axis=1)
            u_x_kosst_slice = u_x_kosst_slice_full[nearest_indices_slice_kosst]

            u_x_kosst_mean = np.mean(np.abs(u_x_kosst_slice - u_x_les_slice))

            # Compute normalised fields for visualisation
            normalised_velocity = u_x_current / UB
            normalised_error_velocity = (u_x_current - u_x_les_slice) / u_x_kosst_mean
            normalised_error_vorticity = (
                vort_x_current - vort_x_les_slice
            ) / avg_diff_vort_kosst

            # Get velocity components for streamlines if available
            u_y_slice = None
            u_z_slice = None
            if "U" in sliced_mesh.point_data:
                u_slice_full = np.array(sliced_mesh.point_data["U"])
                u_y_slice = u_slice_full[:, 1]
                u_z_slice = u_slice_full[:, 2]

            # Create matplotlib contour plots
            fig, ax = create_matplotlib_contour_plot(
                sliced_mesh,
                normalised_velocity,
                normalised_error_velocity,
                normalised_error_vorticity,
                norm_error_combined,
                norm_error_velocity,
                norm_error_vorticity,
                u_y_slice,
                u_z_slice,
            )

            # Save plot to postProcessing directory
            output_name = post_processing_dir / f"{work_dir.name}_multi_.png"
            fig.savefig(output_name, dpi=300, transparent=False, bbox_inches='tight')
            plt.close(fig)
            logger.info(f"Figure saved to: {output_name}")
        else:
            logger.info("Skipping multi plot (--plot-multi flag not set)")

        # Create comparison plot with kOmegaSST, HF, and kOmegaSSTPDA
        if plot_comparison:
            logger.info("Creating comparison contour plot (kOmegaSST, HF, kOmegaSSTPDA)...")
            try:
                # Get kOmegaSST velocity and vorticity on slice
                if "U_KOSST" in ref_internal_mesh.point_data:
                    u_kosst_slice_full = np.array(ref_internal_mesh.point_data["U_KOSST"])
                    u_x_kosst_slice = u_kosst_slice_full[nearest_indices_slice_kosst, 0]
                    u_y_kosst_slice = u_kosst_slice_full[nearest_indices_slice_kosst, 1]
                    u_z_kosst_slice = u_kosst_slice_full[nearest_indices_slice_kosst, 2]
                elif "U_KOSST" in ref_internal_mesh.cell_data:
                    u_kosst_slice_full = np.array(ref_internal_mesh.cell_data["U_KOSST"])
                    u_x_kosst_slice = u_kosst_slice_full[nearest_indices_slice_kosst, 0]
                    u_y_kosst_slice = u_kosst_slice_full[nearest_indices_slice_kosst, 1]
                    u_z_kosst_slice = u_kosst_slice_full[nearest_indices_slice_kosst, 2]
                else:
                    u_x_kosst_slice = None
                    u_y_kosst_slice = None
                    u_z_kosst_slice = None

                # Get kOmegaSST vorticity on slice
                if "vort_KOSST" in ref_internal_mesh.point_data:
                    vort_kosst_slice_full = np.array(ref_internal_mesh.point_data["vort_KOSST"])
                    vort_x_kosst_slice = vort_kosst_slice_full[nearest_indices_slice_kosst, 0]
                elif "vorticity_KOSST" in ref_internal_mesh.point_data:
                    vort_kosst_slice_full = np.array(ref_internal_mesh.point_data["vorticity_KOSST"])
                    vort_x_kosst_slice = vort_kosst_slice_full[nearest_indices_slice_kosst, 0]
                elif "vort_KOSST" in ref_internal_mesh.cell_data:
                    vort_kosst_slice_full = np.array(ref_internal_mesh.cell_data["vort_KOSST"])
                    vort_x_kosst_slice = vort_kosst_slice_full[nearest_indices_slice_kosst, 0]
                else:
                    vort_x_kosst_slice = None

                # Get HF velocity and vorticity on slice (use same slice_points as LES)
                if "U_HF" in ref_internal_mesh.point_data:
                    u_hf_slice_full = np.array(ref_internal_mesh.point_data["U_HF"])
                    u_x_hf_slice = u_hf_slice_full[nearest_indices_slice, 0]
                    u_y_hf_slice = u_hf_slice_full[nearest_indices_slice, 1]
                    u_z_hf_slice = u_hf_slice_full[nearest_indices_slice, 2]
                elif "U_HF" in ref_internal_mesh.cell_data:
                    # Interpolate from cell data to slice points
                    u_hf_cell = np.array(ref_internal_mesh.cell_data["U_HF"])
                    u_x_hf_slice = u_hf_cell[nearest_indices_slice, 0]
                    u_y_hf_slice = u_hf_cell[nearest_indices_slice, 1]
                    u_z_hf_slice = u_hf_cell[nearest_indices_slice, 2]
                else:
                    u_x_hf_slice = None
                    u_y_hf_slice = None
                    u_z_hf_slice = None

                # Get HF vorticity on slice
                if "vort_HF" in ref_internal_mesh.point_data:
                    vort_hf_slice_full = np.array(ref_internal_mesh.point_data["vort_HF"])
                    vort_x_hf_slice = vort_hf_slice_full[nearest_indices_slice, 0]
                elif "vorticity_HF" in ref_internal_mesh.point_data:
                    vort_hf_slice_full = np.array(ref_internal_mesh.point_data["vorticity_HF"])
                    vort_x_hf_slice = vort_hf_slice_full[nearest_indices_slice, 0]
                elif "vort_HF" in ref_internal_mesh.cell_data:
                    vort_hf_cell = np.array(ref_internal_mesh.cell_data["vort_HF"])
                    vort_x_hf_slice = vort_hf_cell[nearest_indices_slice, 0]
                else:
                    vort_x_hf_slice = None

                # Compute normalized fields for all three models
                if u_x_kosst_slice is not None:
                    normalised_velocity_kosst = u_x_kosst_slice / UB
                    normalised_error_velocity_kosst = (u_x_kosst_slice - u_x_les_slice) / u_x_kosst_mean
                else:
                    normalised_velocity_kosst = np.zeros_like(u_x_current)
                    normalised_error_velocity_kosst = np.zeros_like(u_x_current)

                if vort_x_kosst_slice is not None:
                    normalised_error_vorticity_kosst = (
                        vort_x_kosst_slice - vort_x_les_slice
                    ) / avg_diff_vort_kosst
                else:
                    normalised_error_vorticity_kosst = np.zeros_like(vort_x_current)

                if u_x_hf_slice is not None:
                    normalised_velocity_hf = u_x_hf_slice / UB
                    normalised_error_velocity_hf = (u_x_hf_slice - u_x_les_slice) / u_x_kosst_mean
                else:
                    normalised_velocity_hf = np.zeros_like(u_x_current)
                    normalised_error_velocity_hf = np.zeros_like(u_x_current)

                if vort_x_hf_slice is not None:
                    normalised_error_vorticity_hf = (
                        vort_x_hf_slice - vort_x_les_slice
                    ) / avg_diff_vort_kosst
                else:
                    normalised_error_vorticity_hf = np.zeros_like(vort_x_current)

                # Create comparison plot
                fig_comp, ax_comp = create_comparison_contour_plot(
                    sliced_mesh,
                    normalised_velocity_kosst,
                    normalised_velocity_hf,
                    normalised_velocity,
                    normalised_error_velocity_kosst,
                    normalised_error_velocity_hf,
                    normalised_error_velocity,
                    normalised_error_vorticity_kosst,
                    normalised_error_vorticity_hf,
                    normalised_error_vorticity,
                    u_y_kosst_slice,
                    u_z_kosst_slice,
                    u_y_hf_slice,
                    u_z_hf_slice,
                    u_y_slice,
                    u_z_slice,
                )

                # Apply tight layout with improved spacing before saving
                fig_comp.tight_layout(pad=1.0, h_pad=0.3, w_pad=0.5)
                
                # Save comparison plot
                comparison_output = post_processing_dir / f"{work_dir.name}_comparison.png"
                fig_comp.savefig(
                    comparison_output, dpi=300, transparent=False, bbox_inches='tight'
                )
                plt.close(fig_comp)
                logger.info(f"Comparison plot saved to: {comparison_output}")

            except (KeyError, ValueError, RuntimeError) as e:
                logger.warning(f"Could not create comparison plot: {e}")
        else:
            logger.info("Skipping comparison plot (--plot-comparison flag not set)")

        # Create barycentric plots
        if plot_barycentric:
            logger.info("Creating barycentric triangle plots...")
            try:
                # Load k field (needed for all calculations)
                if "k" in internal_mesh.cell_data:
                    k = np.array(internal_mesh.cell_data["k"])
                elif "k" in internal_mesh.point_data:
                    k = np.array(internal_mesh.point_data["k"])
                    k = point_to_cell_data(internal_mesh, k)
                else:
                    raise KeyError("k not found in mesh")

                # For current model (kOmegaSSTPDA): Load turbulenceProperties:R directly
                # This is the FULL Reynolds stress tensor written by OpenFOAM, not a linear part
                rij_current = None
                if "turbulenceProperties:R" in internal_mesh.cell_data:
                    rij_current = np.array(internal_mesh.cell_data["turbulenceProperties:R"])
                    logger.info("Using turbulenceProperties:R field directly for kOmegaSSTPDA")
                elif "turbulenceProperties:R" in internal_mesh.point_data:
                    rij_current = np.array(internal_mesh.point_data["turbulenceProperties:R"])
                    rij_current = point_to_cell_data(internal_mesh, rij_current)
                    logger.info("Using turbulenceProperties:R field directly for kOmegaSSTPDA")
                elif "Rij" in internal_mesh.cell_data:
                    rij_current = np.array(internal_mesh.cell_data["Rij"])
                    logger.info("Using Rij field directly for kOmegaSSTPDA")
                elif "Rij" in internal_mesh.point_data:
                    rij_current = np.array(internal_mesh.point_data["Rij"])
                    rij_current = point_to_cell_data(internal_mesh, rij_current)
                    logger.info("Using Rij field directly for kOmegaSSTPDA")
                else:
                    raise KeyError("Neither turbulenceProperties:R nor Rij found in mesh")

                # Reorder from PyVista format to OpenFOAM format
                rij_current = reorder_pyvista_tensor_to_openfoam(rij_current)

                # Ensure rij_current is in (n_cells, 6) format
                if rij_current.ndim == 2 and rij_current.shape[1] == 9:
                    # Convert from 9-component to 6-component format
                    rij_current = np.column_stack([
                        rij_current[:, 0], rij_current[:, 1], rij_current[:, 2],
                        rij_current[:, 4], rij_current[:, 5], rij_current[:, 8],
                    ])

                # Calculate barycentric coordinates for current model (kOmegaSSTPDA)
                # Rij is in (n_points, 6) format
                logger.info(f"kOmegaSSTPDA: Rij shape {rij_current.shape}")
                eta_current, zeta_current, colours_current = barycentric_calculation_kOmegaSST(rij_current)

                # Load kOmegaSST Rij from reference mesh
                # Prefer Rij_KOSST (pre-computed reference data) over generic Rij
                rij_kosst = None
                if "Rij_KOSST" in ref_internal_mesh.cell_data:
                    rij_kosst = np.array(ref_internal_mesh.cell_data["Rij_KOSST"])
                    logger.info("Using Rij_KOSST field for kOmegaSST")
                elif "Rij_KOSST" in ref_internal_mesh.point_data:
                    rij_kosst = np.array(ref_internal_mesh.point_data["Rij_KOSST"])
                    rij_kosst = point_to_cell_data(ref_internal_mesh, rij_kosst)
                    logger.info("Using Rij_KOSST field for kOmegaSST")
                elif "Rij" in ref_internal_mesh.cell_data:
                    rij_kosst = np.array(ref_internal_mesh.cell_data["Rij"])
                    logger.info("Using Rij field for kOmegaSST")
                elif "Rij" in ref_internal_mesh.point_data:
                    rij_kosst = np.array(ref_internal_mesh.point_data["Rij"])
                    rij_kosst = point_to_cell_data(ref_internal_mesh, rij_kosst)
                    logger.info("Using Rij field for kOmegaSST")
                else:
                    logger.warning("Rij_KOSST or Rij not found in reference mesh, skipping kOmegaSST barycentric plot")
                    rij_kosst = None

                # Reorder from PyVista format to OpenFOAM format
                if rij_kosst is not None:
                    rij_kosst = reorder_pyvista_tensor_to_openfoam(rij_kosst)

                # Ensure rij_kosst is in (n_cells, 6) format if it exists
                if rij_kosst is not None and rij_kosst.ndim == 2 and rij_kosst.shape[1] == 9:
                    # Convert from 9-component to 6-component format
                    rij_kosst = np.column_stack([
                        rij_kosst[:, 0], rij_kosst[:, 1], rij_kosst[:, 2],
                        rij_kosst[:, 4], rij_kosst[:, 5], rij_kosst[:, 8],
                    ])

                # Load HF Rij: prefer Rij_HF directly, fall back to computing from Aij_HF + 2/3*k*I
                rij_hf_cc = None  # Store for Rij profiles plot
                if "Rij_HF" in ref_internal_mesh.cell_data:
                    rij_hf_cc = np.array(ref_internal_mesh.cell_data["Rij_HF"])
                    logger.info("Using Rij_HF field directly for High Fidelity")
                elif "Rij_HF" in ref_internal_mesh.point_data:
                    rij_hf_cc = np.array(ref_internal_mesh.point_data["Rij_HF"])
                    rij_hf_cc = point_to_cell_data(ref_internal_mesh, rij_hf_cc)
                    logger.info("Using Rij_HF field directly for High Fidelity")
                elif "Aij_HF" in ref_internal_mesh.cell_data and "k_HF" in ref_internal_mesh.cell_data:
                    # Compute Rij for HF: Rij_HF = Aij_HF + 2/3 * k_HF * I
                    aij_hf = np.array(ref_internal_mesh.cell_data["Aij_HF"])
                    k_hf = np.array(ref_internal_mesh.cell_data["k_HF"])
                    # Reorder Aij_HF from PyVista format to OpenFOAM format
                    aij_hf = reorder_pyvista_tensor_to_openfoam(aij_hf)
                    identity = np.array([[1, 0, 0, 1, 0, 1]])  # OpenFOAM format: [Rxx, Rxy, Rxz, Ryy, Ryz, Rzz]
                    rij_hf_cc = aij_hf + 2.0 / 3.0 * identity * np.reshape(k_hf, (k_hf.shape[0], 1))
                    logger.info("Computed Rij_HF from Aij_HF and k_HF using: Rij_HF = Aij_HF + 2/3 * I * k_HF")
                elif "Aij_HF" in ref_internal_mesh.point_data and "k_HF" in ref_internal_mesh.point_data:
                    aij_hf = np.array(ref_internal_mesh.point_data["Aij_HF"])
                    k_hf = np.array(ref_internal_mesh.point_data["k_HF"])
                    aij_hf = point_to_cell_data(ref_internal_mesh, aij_hf)
                    k_hf = point_to_cell_data(ref_internal_mesh, k_hf)
                    # Reorder Aij_HF from PyVista format to OpenFOAM format
                    aij_hf = reorder_pyvista_tensor_to_openfoam(aij_hf)
                    identity = np.array([[1, 0, 0, 1, 0, 1]])  # OpenFOAM format: [Rxx, Rxy, Rxz, Ryy, Ryz, Rzz]
                    rij_hf_cc = aij_hf + 2.0 / 3.0 * identity * np.reshape(k_hf, (k_hf.shape[0], 1))
                    logger.info("Computed Rij_HF from Aij_HF and k_HF using: Rij_HF = Aij_HF + 2/3 * I * k_HF")
                else:
                    logger.warning("Rij_HF, Aij_HF or k_HF not found in reference mesh, skipping HF barycentric plot")

                # Reorder from PyVista format to OpenFOAM format (only if loaded directly from PyVista)
                # If computed from Aij_HF + k_HF, Aij_HF was already reordered above
                if rij_hf_cc is not None:
                    if "Rij_HF" in ref_internal_mesh.cell_data or "Rij_HF" in ref_internal_mesh.point_data:
                        rij_hf_cc = reorder_pyvista_tensor_to_openfoam(rij_hf_cc)

                # Ensure rij_hf_cc is in (n_cells, 6) format
                if rij_hf_cc is not None and rij_hf_cc.ndim == 2 and rij_hf_cc.shape[1] == 9:
                    # Convert from 9-component to 6-component format
                    rij_hf_cc = np.column_stack([
                        rij_hf_cc[:, 0], rij_hf_cc[:, 1], rij_hf_cc[:, 2],
                        rij_hf_cc[:, 4], rij_hf_cc[:, 5], rij_hf_cc[:, 8],
                    ])

                # Calculate barycentric coordinates for HF
                if rij_hf_cc is not None:
                    # Convert to (6, n_points) format for HF calculation
                    rij_hf_transposed = rij_hf_cc.T  # Shape: (6, n_points)
                    logger.info(f"HF: Rij shape {rij_hf_transposed.shape}")
                    eta_hf, zeta_hf, colours_hf = barycentric_calculation_hf(rij_hf_transposed)
                else:
                    eta_hf = None
                    zeta_hf = None
                    colours_hf = None

                # Calculate barycentric coordinates for kOmegaSST
                # Rij is in (n_points, 6) format
                if rij_kosst is not None:
                    logger.info(f"kOmegaSST: Rij shape {rij_kosst.shape}")
                    eta_kosst, zeta_kosst, colours_kosst = barycentric_calculation_kOmegaSST(rij_kosst)
                else:
                    eta_kosst = None
                    zeta_kosst = None
                    colours_kosst = None

                # Create barycentric plot if we have at least current model data
                if eta_current is not None and zeta_current is not None:
                    # Use empty arrays if kOmegaSST or HF data is missing
                    if eta_kosst is None:
                        eta_kosst = np.array([])
                        zeta_kosst = np.array([])
                        colours_kosst = np.empty((0, 3))
                    if eta_hf is None:
                        eta_hf = np.array([])
                        zeta_hf = np.array([])
                        colours_hf = np.empty((0, 3))

                    barycentric_output = post_processing_dir / f"{work_dir.name}_barycentric.png"
                    create_barycentric_plot(
                        eta_kosst,
                        zeta_kosst,
                        colours_kosst,
                        eta_current,
                        zeta_current,
                        colours_current,
                        eta_hf,
                        zeta_hf,
                        colours_hf,
                        barycentric_output,
                    )
                else:
                    logger.warning("Cannot create barycentric plot: missing required fields")
            except (KeyError, ValueError, RuntimeError) as e:
                logger.warning(f"Could not create barycentric plot: {e}")
        else:
            logger.info("Skipping barycentric plot (--plot-barycentric flag not set)")
            # Initialize Rij variables to None if barycentric is disabled
            # They will be loaded below if needed for other plots
            rij_current = None
            rij_kosst = None
            rij_hf_cc = None

        # Load Rij data if needed for other plots (anisotropy, velocity profiles, Rij profiles)
        # but barycentric plot is disabled
        if (plot_anisotropy or plot_velocity_profiles or plot_rij_profiles) and not plot_barycentric:
            logger.info("Loading Rij data for other plots (barycentric plot disabled)...")
            try:
                # Load k field
                if "k" in internal_mesh.cell_data:
                    k = np.array(internal_mesh.cell_data["k"])
                elif "k" in internal_mesh.point_data:
                    k = np.array(internal_mesh.point_data["k"])
                    k = point_to_cell_data(internal_mesh, k)
                else:
                    raise KeyError("k not found in mesh")

                # Load current model Rij
                if "turbulenceProperties:R" in internal_mesh.cell_data:
                    rij_current = np.array(internal_mesh.cell_data["turbulenceProperties:R"])
                elif "turbulenceProperties:R" in internal_mesh.point_data:
                    rij_current = np.array(internal_mesh.point_data["turbulenceProperties:R"])
                    rij_current = point_to_cell_data(internal_mesh, rij_current)
                elif "Rij" in internal_mesh.cell_data:
                    rij_current = np.array(internal_mesh.cell_data["Rij"])
                elif "Rij" in internal_mesh.point_data:
                    rij_current = np.array(internal_mesh.point_data["Rij"])
                    rij_current = point_to_cell_data(internal_mesh, rij_current)
                else:
                    raise KeyError("Neither turbulenceProperties:R nor Rij found in mesh")

                rij_current = reorder_pyvista_tensor_to_openfoam(rij_current)
                if rij_current.ndim == 2 and rij_current.shape[1] == 9:
                    rij_current = np.column_stack([
                        rij_current[:, 0], rij_current[:, 1], rij_current[:, 2],
                        rij_current[:, 4], rij_current[:, 5], rij_current[:, 8],
                    ])

                # Load kOmegaSST Rij
                if "Rij_KOSST" in ref_internal_mesh.cell_data:
                    rij_kosst = np.array(ref_internal_mesh.cell_data["Rij_KOSST"])
                elif "Rij_KOSST" in ref_internal_mesh.point_data:
                    rij_kosst = np.array(ref_internal_mesh.point_data["Rij_KOSST"])
                    rij_kosst = point_to_cell_data(ref_internal_mesh, rij_kosst)
                elif "Rij" in ref_internal_mesh.cell_data:
                    rij_kosst = np.array(ref_internal_mesh.cell_data["Rij"])
                elif "Rij" in ref_internal_mesh.point_data:
                    rij_kosst = np.array(ref_internal_mesh.point_data["Rij"])
                    rij_kosst = point_to_cell_data(ref_internal_mesh, rij_kosst)
                else:
                    rij_kosst = None

                if rij_kosst is not None:
                    rij_kosst = reorder_pyvista_tensor_to_openfoam(rij_kosst)
                    if rij_kosst.ndim == 2 and rij_kosst.shape[1] == 9:
                        rij_kosst = np.column_stack([
                            rij_kosst[:, 0], rij_kosst[:, 1], rij_kosst[:, 2],
                            rij_kosst[:, 4], rij_kosst[:, 5], rij_kosst[:, 8],
                        ])

                # Load HF Rij
                if "Rij_HF" in ref_internal_mesh.cell_data:
                    rij_hf_cc = np.array(ref_internal_mesh.cell_data["Rij_HF"])
                elif "Rij_HF" in ref_internal_mesh.point_data:
                    rij_hf_cc = np.array(ref_internal_mesh.point_data["Rij_HF"])
                    rij_hf_cc = point_to_cell_data(ref_internal_mesh, rij_hf_cc)
                elif "Aij_HF" in ref_internal_mesh.cell_data and "k_HF" in ref_internal_mesh.cell_data:
                    aij_hf = np.array(ref_internal_mesh.cell_data["Aij_HF"])
                    k_hf = np.array(ref_internal_mesh.cell_data["k_HF"])
                    aij_hf = reorder_pyvista_tensor_to_openfoam(aij_hf)
                    identity = np.array([[1, 0, 0, 1, 0, 1]])
                    rij_hf_cc = aij_hf + 2.0 / 3.0 * identity * np.reshape(k_hf, (k_hf.shape[0], 1))
                elif "Aij_HF" in ref_internal_mesh.point_data and "k_HF" in ref_internal_mesh.point_data:
                    aij_hf = np.array(ref_internal_mesh.point_data["Aij_HF"])
                    k_hf = np.array(ref_internal_mesh.point_data["k_HF"])
                    aij_hf = point_to_cell_data(ref_internal_mesh, aij_hf)
                    k_hf = point_to_cell_data(ref_internal_mesh, k_hf)
                    aij_hf = reorder_pyvista_tensor_to_openfoam(aij_hf)
                    identity = np.array([[1, 0, 0, 1, 0, 1]])
                    rij_hf_cc = aij_hf + 2.0 / 3.0 * identity * np.reshape(k_hf, (k_hf.shape[0], 1))
                else:
                    rij_hf_cc = None

                if rij_hf_cc is not None:
                    if "Rij_HF" in ref_internal_mesh.cell_data or "Rij_HF" in ref_internal_mesh.point_data:
                        rij_hf_cc = reorder_pyvista_tensor_to_openfoam(rij_hf_cc)
                    if rij_hf_cc.ndim == 2 and rij_hf_cc.shape[1] == 9:
                        rij_hf_cc = np.column_stack([
                            rij_hf_cc[:, 0], rij_hf_cc[:, 1], rij_hf_cc[:, 2],
                            rij_hf_cc[:, 4], rij_hf_cc[:, 5], rij_hf_cc[:, 8],
                        ])

                logger.info("Rij data loaded for other plots")
            except (KeyError, ValueError, RuntimeError) as e:
                logger.warning(f"Could not load Rij data for other plots: {e}")
                # Set to None so dependent plots will skip
                rij_current = None
                rij_kosst = None
                rij_hf_cc = None

        # Create anisotropy contour plot
        if plot_anisotropy:
            logger.info("Creating anisotropy contour plot...")
            try:
                # Check if Rij data is available
                if 'rij_current' not in locals() or rij_current is None:
                    raise KeyError("rij_current not available for anisotropy plot")
                
                # Interpolate Rij data from full mesh to slice for each model
                slice_points = sliced_mesh.points

                # Current model (kOmegaSSTPDA) - interpolate from internal mesh
                current_cell_centers = internal_mesh.cell_centers().points
                distances_current = cdist(slice_points, current_cell_centers)
                nearest_current = np.argmin(distances_current, axis=1)
                rij_current_slice = rij_current[nearest_current]

                # kOmegaSST - interpolate from reference mesh
                if rij_kosst is not None:
                    ref_cell_centers = ref_internal_mesh.cell_centers().points
                    distances_kosst = cdist(slice_points, ref_cell_centers)
                    nearest_kosst = np.argmin(distances_kosst, axis=1)
                    rij_kosst_slice = rij_kosst[nearest_kosst]
                else:
                    rij_kosst_slice = None

                # HF - interpolate from reference mesh
                if rij_hf_cc is not None:
                    # rij_hf_cc is (n_cells, 6), interpolate to slice
                    distances_hf = cdist(slice_points, ref_cell_centers)
                    nearest_hf = np.argmin(distances_hf, axis=1)
                    rij_hf_slice = rij_hf_cc[nearest_hf]
                else:
                    rij_hf_slice = None

                anisotropy_output = post_processing_dir / f"{work_dir.name}_anisotropy_contour.png"
                create_anisotropy_contour_plot(
                    sliced_mesh,
                    rij_kosst_slice,
                    rij_current_slice,
                    rij_hf_slice,
                    anisotropy_output,
                )
            except (KeyError, ValueError, RuntimeError) as e:
                logger.warning(f"Could not create anisotropy contour plot: {e}")
        else:
            logger.info("Skipping anisotropy contour plot (--plot-anisotropy flag not set)")

        # Create velocity profiles plot
        if plot_velocity_profiles:
            logger.info("Creating velocity profiles plot...")
            try:
                # Get velocity fields at cell centers for all models
                # Current model
                if "U" in internal_mesh.point_data:
                    u_current_cc = point_to_cell_data(internal_mesh, u_full)
                else:
                    u_current_cc = np.array(internal_mesh.cell_data["U"])

                # kOmegaSST
                if "U_KOSST" in ref_internal_mesh.point_data:
                    u_kosst_cc = point_to_cell_data(ref_internal_mesh, u_kosst_full)
                else:
                    u_kosst_cc = np.array(ref_internal_mesh.cell_data["U_KOSST"])

                # HF
                if "U_HF" in ref_internal_mesh.point_data:
                    u_hf_cc = point_to_cell_data(ref_internal_mesh, u_les_full)
                else:
                    u_hf_cc = np.array(ref_internal_mesh.cell_data["U_HF"])

                profiles_output = post_processing_dir / f"{work_dir.name}_velocity_profiles.png"
                create_velocity_profiles_plot(
                    internal_mesh,
                    ref_internal_mesh,
                    u_current_cc,
                    u_kosst_cc,
                    u_hf_cc,
                    profiles_output,
                )
            except (KeyError, ValueError, RuntimeError) as e:
                logger.warning(f"Could not create velocity profiles plot: {e}")
        else:
            logger.info("Skipping velocity profiles plot (--plot-velocity-profiles flag not set)")

        # Create Rij component profiles plot
        if plot_rij_profiles:
            logger.info("Creating Rij component profiles plot...")
            try:
                # Reuse Rij data from barycentric section if available
                # These should be computed in the barycentric try block above
                if 'rij_current' in locals() and rij_current is not None:
                    rij_current_cc = rij_current
                else:
                    raise KeyError("rij_current not available from barycentric section")

                if 'rij_kosst' in locals() and rij_kosst is not None:
                    rij_kosst_cc = rij_kosst
                else:
                    raise KeyError("rij_kosst not available from barycentric section")

                if 'rij_hf_cc' in locals() and rij_hf_cc is not None:
                    # Already in cell format
                    pass
                else:
                    raise KeyError("rij_hf_cc not available from barycentric section")

                rij_profiles_output = post_processing_dir / f"{work_dir.name}_rij_profiles.png"
                create_rij_profiles_plot(
                    internal_mesh,
                    ref_internal_mesh,
                    rij_current_cc,
                    rij_kosst_cc,
                    rij_hf_cc,
                    rij_profiles_output,
                )
            except (KeyError, ValueError, RuntimeError) as e:
                logger.warning(f"Could not create Rij profiles plot: {e}")
        else:
            logger.info("Skipping Rij profiles plot (--plot-rij-profiles flag not set)")

        # Create diagonal velocity profile plot
        if plot_diagonal_profile:
            logger.info("Creating diagonal velocity profile plot...")
            try:
                # Get velocity fields at cell centers for all models
                # Current model
                if "U" in internal_mesh.point_data:
                    u_current_cc = point_to_cell_data(internal_mesh, u_full)
                else:
                    u_current_cc = np.array(internal_mesh.cell_data["U"])

                # kOmegaSST
                if "U_KOSST" in ref_internal_mesh.point_data:
                    u_kosst_cc = point_to_cell_data(ref_internal_mesh, u_kosst_full)
                else:
                    u_kosst_cc = np.array(ref_internal_mesh.cell_data["U_KOSST"])

                # HF
                if "U_HF" in ref_internal_mesh.point_data:
                    u_hf_cc = point_to_cell_data(ref_internal_mesh, u_les_full)
                else:
                    u_hf_cc = np.array(ref_internal_mesh.cell_data["U_HF"])

                diagonal_output = post_processing_dir / f"{work_dir.name}_diagonal_profile.png"
                create_diagonal_velocity_profile_plot(
                    internal_mesh,
                    ref_internal_mesh,
                    u_current_cc,
                    u_kosst_cc,
                    u_hf_cc,
                    diagonal_output,
                )
            except (KeyError, ValueError, RuntimeError) as e:
                logger.warning(f"Could not create diagonal velocity profile plot: {e}")
        else:
            logger.info("Skipping diagonal velocity profile plot (--plot-diagonal-profile flag not set)")

        # Note: Barycentric data loading is needed for anisotropy, velocity profiles, and Rij profiles
        # So we need to load it even if barycentric plot is disabled
        if (plot_anisotropy or plot_velocity_profiles or plot_rij_profiles) and not plot_barycentric:
            logger.info("Loading Rij data for other plots (barycentric plot disabled)...")
            try:
                # Load k field
                if "k" in internal_mesh.cell_data:
                    k = np.array(internal_mesh.cell_data["k"])
                elif "k" in internal_mesh.point_data:
                    k = np.array(internal_mesh.point_data["k"])
                    k = point_to_cell_data(internal_mesh, k)
                else:
                    raise KeyError("k not found in mesh")

                # Load current model Rij
                rij_current = None
                if "turbulenceProperties:R" in internal_mesh.cell_data:
                    rij_current = np.array(internal_mesh.cell_data["turbulenceProperties:R"])
                elif "turbulenceProperties:R" in internal_mesh.point_data:
                    rij_current = np.array(internal_mesh.point_data["turbulenceProperties:R"])
                    rij_current = point_to_cell_data(internal_mesh, rij_current)
                elif "Rij" in internal_mesh.cell_data:
                    rij_current = np.array(internal_mesh.cell_data["Rij"])
                elif "Rij" in internal_mesh.point_data:
                    rij_current = np.array(internal_mesh.point_data["Rij"])
                    rij_current = point_to_cell_data(internal_mesh, rij_current)
                else:
                    raise KeyError("Neither turbulenceProperties:R nor Rij found in mesh")

                rij_current = reorder_pyvista_tensor_to_openfoam(rij_current)
                if rij_current.ndim == 2 and rij_current.shape[1] == 9:
                    rij_current = np.column_stack([
                        rij_current[:, 0], rij_current[:, 1], rij_current[:, 2],
                        rij_current[:, 4], rij_current[:, 5], rij_current[:, 8],
                    ])

                # Load kOmegaSST Rij
                rij_kosst = None
                if "Rij_KOSST" in ref_internal_mesh.cell_data:
                    rij_kosst = np.array(ref_internal_mesh.cell_data["Rij_KOSST"])
                elif "Rij_KOSST" in ref_internal_mesh.point_data:
                    rij_kosst = np.array(ref_internal_mesh.point_data["Rij_KOSST"])
                    rij_kosst = point_to_cell_data(ref_internal_mesh, rij_kosst)
                elif "Rij" in ref_internal_mesh.cell_data:
                    rij_kosst = np.array(ref_internal_mesh.cell_data["Rij"])
                elif "Rij" in ref_internal_mesh.point_data:
                    rij_kosst = np.array(ref_internal_mesh.point_data["Rij"])
                    rij_kosst = point_to_cell_data(ref_internal_mesh, rij_kosst)

                if rij_kosst is not None:
                    rij_kosst = reorder_pyvista_tensor_to_openfoam(rij_kosst)
                    if rij_kosst.ndim == 2 and rij_kosst.shape[1] == 9:
                        rij_kosst = np.column_stack([
                            rij_kosst[:, 0], rij_kosst[:, 1], rij_kosst[:, 2],
                            rij_kosst[:, 4], rij_kosst[:, 5], rij_kosst[:, 8],
                        ])

                # Load HF Rij
                rij_hf_cc = None
                if "Rij_HF" in ref_internal_mesh.cell_data:
                    rij_hf_cc = np.array(ref_internal_mesh.cell_data["Rij_HF"])
                elif "Rij_HF" in ref_internal_mesh.point_data:
                    rij_hf_cc = np.array(ref_internal_mesh.point_data["Rij_HF"])
                    rij_hf_cc = point_to_cell_data(ref_internal_mesh, rij_hf_cc)
                elif "Aij_HF" in ref_internal_mesh.cell_data and "k_HF" in ref_internal_mesh.cell_data:
                    aij_hf = np.array(ref_internal_mesh.cell_data["Aij_HF"])
                    k_hf = np.array(ref_internal_mesh.cell_data["k_HF"])
                    aij_hf = reorder_pyvista_tensor_to_openfoam(aij_hf)
                    identity = np.array([[1, 0, 0, 1, 0, 1]])
                    rij_hf_cc = aij_hf + 2.0 / 3.0 * identity * np.reshape(k_hf, (k_hf.shape[0], 1))
                elif "Aij_HF" in ref_internal_mesh.point_data and "k_HF" in ref_internal_mesh.point_data:
                    aij_hf = np.array(ref_internal_mesh.point_data["Aij_HF"])
                    k_hf = np.array(ref_internal_mesh.point_data["k_HF"])
                    aij_hf = point_to_cell_data(ref_internal_mesh, aij_hf)
                    k_hf = point_to_cell_data(ref_internal_mesh, k_hf)
                    aij_hf = reorder_pyvista_tensor_to_openfoam(aij_hf)
                    identity = np.array([[1, 0, 0, 1, 0, 1]])
                    rij_hf_cc = aij_hf + 2.0 / 3.0 * identity * np.reshape(k_hf, (k_hf.shape[0], 1))

                if rij_hf_cc is not None:
                    if "Rij_HF" in ref_internal_mesh.cell_data or "Rij_HF" in ref_internal_mesh.point_data:
                        rij_hf_cc = reorder_pyvista_tensor_to_openfoam(rij_hf_cc)
                    if rij_hf_cc.ndim == 2 and rij_hf_cc.shape[1] == 9:
                        rij_hf_cc = np.column_stack([
                            rij_hf_cc[:, 0], rij_hf_cc[:, 1], rij_hf_cc[:, 2],
                            rij_hf_cc[:, 4], rij_hf_cc[:, 5], rij_hf_cc[:, 8],
                        ])

                logger.info("Rij data loaded for other plots")
            except (KeyError, ValueError, RuntimeError) as e:
                logger.warning(f"Could not load Rij data for other plots: {e}")
                # Set to None so dependent plots will skip
                rij_current = None
                rij_kosst = None
                rij_hf_cc = None

    else:
        logger.warning("Cannot create visualisations: mesh data not available")


if __name__ == "__main__":
    main()
