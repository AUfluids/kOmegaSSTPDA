#!/usr/bin/env python3
"""
plot_residuals.py

Description:    Plot OpenFOAM residuals including linear residuals and continuity residuals
                for Channel Flow test case at Re_tau = 5200.

Author:         Mario Javier Rincon, PhD
Affiliation:    Aarhus University/Kamstrup A/S
Contact:        mjrp@mpe.au.dk/mjrp@kamstrup.com
Version:        2.0.0
Last Updated:   2024-01-01
"""

import logging
from pathlib import Path
from typing import List, Optional

import click
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from rich.logging import RichHandler

# Configure logging using RichHandler for coloured output
logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)],
)
logger = logging.getLogger(__name__)

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

# Figure size constants
cm = 1 / 2.54  # inches to cm
textWidth = 9  # cm
figSize = (textWidth * cm, textWidth * cm * 3 / 4)
figSizeMedium = (14 * cm, 14 * cm * 3 / 4)
figSizeFull = (19 * cm, 19 * cm)

# Colour palette
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


def load_residual_data(
    log_dir: Path,
    filename: str,
    headers: List[str],
) -> Optional[pd.DataFrame]:
    """
    Load residual data from a log file.

    Args:
        log_dir: Directory containing log files
        filename: Name of the log file
        headers: Column headers for the data

    Returns:
        DataFrame containing the data, or None if loading failed
    """
    file_path = log_dir / filename

    if not file_path.exists():
        logger.warning(f"File not found: {file_path}")
        return None

    if not file_path.is_file():
        logger.warning(f"Path is not a file: {file_path}")
        return None

    try:
        df = pd.read_csv(
            file_path,
            delimiter=r'\s+',
            names=headers,
            comment='#',
        )

        if df.empty:
            logger.warning(f"Empty file: {file_path}")
            return None

        return df

    except pd.errors.EmptyDataError:
        logger.warning(f"Empty data file: {file_path}")
        return None
    except Exception as e:
        logger.error(f"Error loading {file_path}: {e}")
        return None


def plot_linear_residuals(
    ax: plt.Axes,
    log_dir: Path,
    variable_names: List[str],
    legend_names: List[str],
) -> bool:
    """
    Plot linear residuals on the given axes.

    Args:
        ax: Matplotlib axes to plot on
        log_dir: Directory containing log files
        variable_names: List of variable file names
        legend_names: List of LaTeX-formatted legend names

    Returns:
        True if at least one variable was plotted successfully
    """
    plotted = False

    for var_name, legend_name in zip(variable_names, legend_names):
        headers = ['Time', var_name]
        df = load_residual_data(log_dir, var_name, headers)

        if df is not None:
            try:
                ax.plot(
                    df['Time'],
                    df[var_name],
                    label=legend_name,
                    linewidth=0.75,
                )
                plotted = True
            except KeyError as e:
                logger.warning(f"Missing column in {var_name}: {e}")
            except Exception as e:
                logger.error(f"Error plotting {var_name}: {e}")

    if plotted:
        ax.set_ylabel("Initial residual")
        ax.set_xlabel("Iteration")
        ax.set_yscale('log')
        ax.legend(
            loc='center left',
            bbox_to_anchor=(1.0, 0.5),
            frameon=True,
            fancybox=True,
        )
        ax.grid(True, alpha=0.3)

    return plotted


def plot_continuity_residuals(
    ax: plt.Axes,
    log_dir: Path,
) -> bool:
    """
    Plot continuity residuals on the given axes.

    Args:
        ax: Matplotlib axes to plot on
        log_dir: Directory containing log files

    Returns:
        True if continuity residuals were plotted successfully
    """
    continuity_vars = ['contGlobal_0', 'contCumulative_0', 'contLocal_0']
    continuity_labels = ['Global', 'Cumulative', 'Local']

    # Load continuity data
    continuity_data = {}
    for var_name in continuity_vars:
        headers = ['Time', var_name]
        df = load_residual_data(log_dir, var_name, headers)
        if df is not None:
            continuity_data[var_name] = df

    if not continuity_data:
        logger.warning("No continuity residual data found")
        return False

    # Plot cumulative and global on primary y-axis
    if 'contCumulative_0' in continuity_data:
        df_cum = continuity_data['contCumulative_0']
        ax.plot(
            df_cum['Time'],
            df_cum['contCumulative_0'],
            color=colours[3],  # green
            label=continuity_labels[1],
            linewidth=0.75,
            zorder=2,
        )

    if 'contGlobal_0' in continuity_data:
        df_global = continuity_data['contGlobal_0']
        ax.plot(
            df_global['Time'],
            df_global['contGlobal_0'],
            color=colours[1],  # red
            label=continuity_labels[0],
            linewidth=0.75,
            linestyle='--',
            zorder=2,
        )

        # Add zero reference line
        if len(df_global) > 0:
            ax.axhline(
                y=0,
                color='k',
                linestyle=':',
                linewidth=0.5,
                zorder=10,
            )

    ax.set_xlabel("Iteration")
    ax.set_ylabel("Continuity residual", color=colours[1])
    ax.tick_params(axis='y', labelcolor=colours[1])
    # Format y-axis in scientific notation with exponents
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0), useMathText=True)
    ax.grid(True, alpha=0.3)

    # Create secondary y-axis for local continuity if available
    ax2 = None
    if 'contLocal_0' in continuity_data:
        ax2 = ax.twinx()
        df_local = continuity_data['contLocal_0']
        ax2.plot(
            df_local['Time'],
            df_local['contLocal_0'],
            color=colours[0],  # navy blue
            label=continuity_labels[2],
            linewidth=0.75,
            zorder=2,
        )
        ax2.set_ylabel(continuity_labels[2], color=colours[0])
        ax2.tick_params(axis='y', labelcolor=colours[0])
        # Format secondary y-axis in scientific notation with exponents
        ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0), useMathText=True)

    # Place legend inside the plot in the best position
    # Combine legends from both axes if ax2 exists
    if ax2 is not None:
        # Get handles and labels from both axes
        handles1, labels1 = ax.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        # Combine them
        handles = handles1 + handles2
        labels = labels1 + labels2
        # Create combined legend inside plot (upper left to avoid secondary axis label)
        ax.legend(
            handles,
            labels,
            loc='upper left',
            frameon=True,
            fancybox=True,
        )
    else:
        ax.legend(
            loc='best',
            frameon=True,
            fancybox=True,
        )

    return True


@click.command()
@click.option(
    '--plot-continuity',
    is_flag=True,
    default=True,
    help='Plot continuity residuals [default: True]'
)
def main(plot_continuity: bool) -> None:
    """
    Main function to plot OpenFOAM residuals.
    
    Creates plots for linear residuals (velocity components, pressure, k, omega)
    and optionally continuity residuals (global, cumulative, local).
    """
    work_dir = Path.cwd()
    log_dir = work_dir / "logs"
    post_processing_dir = work_dir / "postProcessing"

    logger.info(f"Working directory: {work_dir}")
    logger.info(f"Plot options - Continuity residuals: {plot_continuity}")

    # Check if logs directory exists
    if not log_dir.exists():
        logger.error(f"Logs directory not found: {log_dir}")
        return

    if not log_dir.is_dir():
        logger.error(f"Path is not a directory: {log_dir}")
        return

    # Create postProcessing directory if it doesn't exist
    try:
        post_processing_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Error creating postProcessing directory: {e}")
        return

    # Define variables to plot (matching CF_5200: Ux, Uy, p, k, omega)
    variable_names = ['Ux_0', 'Uy_0', 'p_0', 'k_0', 'omega_0']
    legend_names = ['$u$', '$v$', '$p$', '$k$', r'$\omega$']

    # Determine number of subplots
    n_plots = 1 if not plot_continuity else 2
    height_ratios = [2, 1] if plot_continuity else [1]

    # Create figure with subplots
    fig, axes = plt.subplots(
        nrows=n_plots,
        ncols=1,
        figsize=figSizeMedium,
        sharex=True,
        height_ratios=height_ratios if plot_continuity else None,
    )

    # Ensure axes is always a list
    if n_plots == 1:
        axes = [axes]

    # Plot linear residuals
    logger.info("Plotting linear residuals...")
    linear_plotted = plot_linear_residuals(
        axes[0],
        log_dir,
        variable_names,
        legend_names,
    )

    # Plot continuity residuals if requested
    continuity_plotted = False
    if plot_continuity:
        logger.info("Plotting continuity residuals...")
        continuity_plotted = plot_continuity_residuals(axes[1], log_dir)
    else:
        logger.info("Skipping continuity residuals (--plot-continuity flag not set)")

    if not linear_plotted and not continuity_plotted:
        logger.error("No residual data could be loaded or plotted")
        plt.close(fig)
        return

    # Adjust layout
    fig.tight_layout()

    # Save figure
    output_path = post_processing_dir / "residuals.png"
    try:
        fig.savefig(
            output_path,
            dpi=300,
            bbox_inches='tight',
            pad_inches=0.1,
            transparent=False,
        )
        logger.info(f"Residuals plot saved to: {output_path}")
    except Exception as e:
        logger.error(f"Error saving plot: {e}")
        plt.close(fig)
        return

    plt.close(fig)


if __name__ == "__main__":
    main()
