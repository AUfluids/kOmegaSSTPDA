# kOmegaSSTPDA
A generalisable data-augmented k-omega SST model with progressive and interpretable corrections for OpenFOAM. Developed by Rincón et al. 2025 and the Fluid Physics & Turbulence research group at Aarhus University.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Model Configuration](#model-configuration)
- [Validation](#validation)
- [License](#license)
- [Citations](#how-to-cite)
- [Authors](#authors)
- [References](#references)
- [Disclaimer](#disclaimer)

## Overview
The k-omega-SST-PDA model is a progressive data-augmented turbulence model that combines Bayesian optimisation with physics-guided corrections to improve predictions of anisotropy-induced secondary flows and flow separation simultaneously. The model features:

### Key Features
- **Progressive Data Augmentation (PDA)**: A framework that systematically embeds interpretable modifications through Bayesian optimisation
- **Separation Flow Enhancement**: 
  - Activation-based separation correction in the omega-equation
  - Optimised power-law function for local turbulent viscosity adjustment
  - Improved predictions in adverse pressure gradient regions
- **Anisotropy-induced Secondary Flow Prediction**:
  - Non-linear Reynolds stress anisotropy correction
  - Enhanced prediction of Prandtl's second kind of secondary flows
  - Improved corner flow and streamwise vorticity predictions

### Implementation Details
- Optimised for 2D flows
- Verified in 3D flows
- Configurable through user-modifiable coefficients
- Default optimised values provided
- Compatible with incompressible and compressible flows

## Requirements
- OpenFOAM-v2412 or previous ESI versions
- C++11 or later

## Installation
1. Clone the repository:
     ```
     git clone https://github.com/AUfluids/KOSSTPDA.git
     ```

2. Make the installation script executable:
     ```
     cd KOSSTPDA
     chmod a+x Allwmake
     ```

3. Compile the model:
     ```
     ./Allwmake
     ```

## Usage
1. Add the required library to `controlDict`:
   - For incompressible flow:
     ```
     libs ( "libPDAIncompressibleTurbulenceModels" );
     ```
   - For compressible flow:
     ```
     libs ( "libPDAcompressibleTurbulenceModels" );
     ```

2. Specify in `turbulentProperties`:
   ```
   RASModel kOmegaSSTPDA;
   ```

3. Add to `system/fvSchemes`:
   ```
   divSchemes
   {
       div(dev(((2*k)*bijDelta)))    Gauss linear;
   }
   ```

> **Note**: If the solution does not converge at first, it is recommended to first run your case with standard kOmegaSST before switching to kOmegaSSTPDA.

## Model Configuration
### Mode Selection
By default, the model activates both secondary and separation effects. If desired, one can change the models as follows:
   ```
   separationCorrection true; // default: true - off: false
   secondaryCorrection true; // default: true - off: false
   ```

If you use false in both corrections, the PDA model is deactivated, and the standard kOmegaSST is used.

### Optional Stability Settings
In case of stability and convergence issues, we also suggest the following setting for the new model to be tested.
Otherwise, these coefficients are automatically assigned values corresponding to the optimised model.

   ```
   // Separation Flow coefficients
   separationCorrection true;
   C0 -1;
   C1 0;
   C2 0;
   lambda1 1;
   lambda2 1;
   // Anisotropy Secondary Flow coefficients
   anisotropyCorrection true;
   A0 -1;
   A1 0;
   A2 0;
   anisotropyRelaxation 0.6;  // Relaxation factor for more stable simulations
   ```

## Validation
The model has been validated across multiple test cases:

### Separation Flow Cases
- Periodic hills
- Curved backward-facing step (Reb = 13700)
- Converging-diverging channel
- Parametric bumps

### Anisotropy-induced Secondary Flow Cases
- Duct flow (AR = 1, Reb = 3500)
- Duct flow (AR = 3, Reb = 2600)
- Roughness-induced atmospheric boundary layer flows

### Results
#### Law-of-the-wall
Results for channel flow (Re_tau = 590):
![Law-of-the-wall](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/CF_590/CF_u_CF_590.png)

#### Separation
Results for periodic hill (Reb = 2800):
![Velocity Contours](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/PH_2800/contours_comparison_PH_2800.png)
![Velocity Profiles Comparison](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/PH_2800/profiles_comparison_PH_2800.png)

#### Anisotropy-induced Secondary Flow
Results for duct flow (AR = 1, Reb = 3500):
![Anisotropy Secondary Flow](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/SD_ReB3500_AR1/contours_comparison_SD_ReB3500_AR1.png)
![Velocity Profiles](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/SD_ReB3500_AR1/profiles_comparison_SD_ReB3500_AR1.png)
![Reynolds Stress Profiles](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/SD_ReB3500_AR1/profiles_comparison_Rij_SD_ReB3500_AR1.png)

## Target platform
The code is known to work with OpenFOAM-v2406 and previous ESI versions.

## Authors
Mario Javier Rincón <mjrp@mpe.au.dk>

Ali Amarloo <amarloo@mpe.au.dk>

## References
For more details about the model development and validation, refer to:
- [A generalisable data-augmented turbulence model with progressive and interpretable corrections](https://arxiv.org/abs/2503.18568)
- [Progressive augmentation of turbulence models for flow separation by multi-case computational fluid dynamics driven surrogate optimization](https://doi.org/10.1063/5.0174470)
- [Progressive augmentation of Reynolds stress tensor models for secondary flow prediction by computational fluid dynamics driven surrogate optimisation](https://doi.org/10.1016/j.ijheatfluidflow.2023.109242)

## How to cite
Please cite this library using the following publications:

Rincón et al. (2025)
```
@article{rincon2025generalisable,
  title={A generalisable data-augmented turbulence model with progressive and interpretable corrections},
  author={Rinc{\'o}n, Mario J and Reclari, Martino and Yang, Xiang IA and Abkar, Mahdi},
  journal={arXiv preprint arXiv:2503.18568},
  year={2025}
}
```

Amarloo and Rincón (2023)
```
@article{amarloo2023progressive,
  title={Progressive augmentation of turbulence models for flow separation by multi-case computational fluid dynamics driven surrogate optimization},
  author={Amarloo, Ali and Rinc{\'o}n, Mario Javier and Reclari, Martino and Abkar, Mahdi},
  journal={Physics of Fluids},
  volume={35},
  number={12},
  year={2023},
  publisher={AIP Publishing}
}
```

Rincón and Amarloo (2023)
```
@article{rincon2023progressive,
  title={Progressive augmentation of Reynolds stress tensor models for secondary flow prediction by computational fluid dynamics driven surrogate optimisation},
  journal={International Journal of Heat and Fluid Flow},
  volume={104},
  pages={109242},
  year={2023},
  issn={0142-727X},
  doi={https://doi.org/10.1016/j.ijheatfluidflow.2023.109242},
  author={Mario Javier Rincón and Ali Amarloo and Martino Reclari and Xiang I.A. Yang and Mahdi Abkar}
}
```

## Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, the producer of the OpenFOAM software and owner of the OPENFOAM® and OpenCFD® trade marks.

Detailed information on the OpenFOAM trademark can be found at:
- http://www.openfoam.com/legal/trademark-policy.php
- http://www.openfoam.com/legal/trademark-guidelines.php

For further information on OpenCFD and OpenFOAM, please refer to:
- http://www.openfoam.com
