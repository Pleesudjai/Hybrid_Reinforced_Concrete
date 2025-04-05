# Hybrid Reinforced Concrete (HRC) Analysis Program

## Overview

This program analyzes the bending behavior of Hybrid Reinforced Concrete (HRC) beams, which combine fiber reinforced concrete (FRC) with conventional steel reinforcement. The program calculates moment-curvature relationships, load-deflection curves, and efficiency factors for different material combinations.

## Main Features

- Trilinear tension model and bilinear compression model for concrete
- Steel reinforcement modeling (both tension and compression)
- Detailed moment-curvature analysis with efficiency factors
- Calculation of key transition points in the moment-curvature relationship
- Load-deflection prediction for 3-point or 4-point bending tests
- Comprehensive visualization of results with multiple plots
- Different zone calculations based on material states (elastic/plastic)

## Files and Structure

### Main Program

- `MainProgram_HRC_Analysis.m` - The main program that controls the overall analysis

### Transition Point Calculator

- `Intersection_points_Tri_HRC.m` - Calculates key transition points in moment-curvature curve

### Zone Analysis Functions

For different combinations of concrete and steel behavior:

- `zone1_2024.m` - Elastic concrete in compression, elastic concrete in tension, elastic steel
- `zone21_2024.m` - Elastic concrete in compression, concrete in 2nd slope tension, elastic steel
- `zone22_2024.m` - Elastic concrete in compression, concrete in 2nd slope tension, yielded steel
- `zone31_2024.m` - Plastic concrete in compression, concrete in 2nd slope tension, elastic steel
- `zone32_2024.m` - Plastic concrete in compression, concrete in 2nd slope tension, yielded steel
- `zone41_2024.m` - Elastic concrete in compression, concrete in 3rd slope tension, elastic steel
- `zone42_2024.m` - Elastic concrete in compression, concrete in 3rd slope tension, yielded steel
- `zone51_2024.m` - Plastic concrete in compression, concrete in 3rd slope tension, elastic steel
- `zone52_2024.m` - Plastic concrete in compression, concrete in 3rd slope tension, yielded steel

### Utility Functions

- `Left_xbar.m` - Calculates the centroid x-coordinate for deflection calculations

## How to Use

1. Set the input parameters in the "Part A" and "Part B" sections of the main program
2. Run the main program `MainProgram_HRC_Analysis.m`
3. View the generated plots and output files

## Input Parameters

### Geometric Parameters

- `b` - Beam width
- `h` - Beam total depth
- `L` - Clear span
- `alpha` - Depth of steel to total depth ratio

### Test Configuration

- `pointBend` - 3 for 3-point bending, 4 for 4-point bending
- `S2` - Middle spacing for 4-point bending
- `Lp` - Plastic length for localized zone
- `c` - Localized length/spacing ratio

### Material Parameters

#### Concrete Parameters

- `E` - Young's modulus
- `epsilon_cr` - Cracking strain
- `xi` - Compressive Young's modulus parameter
- `omega` - Compressive yield strain parameter
- `lambda_cu` - Ultimate compressive strain parameter
- `mu` - Residual tension strength parameter
- `tau` - Transition zone tensile strain parameter
- `eta` - Calculated tension parameter

#### Steel Parameters

- `rho` - Steel area ratio
- `kappa` - Rebar yield strain parameter
- `n` - Rebar Young's modulus parameter
- `zeta` - Compression/tension steel area ratio

## Output Files

- `Analysis_output.dat` - Main output with all calculated data
- `Efficiency_output.dat` - Efficiency factors at key transition points
- `Intersection_output.dat` - Data for key transition points

## Technical Background

The program handles several different zones in the moment-curvature relationship:

- **Zone 1**: Both concrete and steel are in the elastic range
- **Zone 2**: Concrete in tension follows 2nd slope (post-cracking)
- **Zone 3**: Concrete in compression reaches plastic region
- **Zone 4**: Concrete in tension follows 3rd slope (residual strength)
- **Zone 5**: Both concrete in compression and tension are in their final states

Each zone is further divided based on whether the steel has yielded or not, creating the complete range of possible material states.

## Visualization Outputs

The program generates various plots to help interpret the analysis results:

1. Material models (compression, tension, and rebar)
2. Moment-curvature diagrams
3. Curvature distribution along the beam
4. Moment distribution along the beam
5. Load-deflection curve
6. Equivalent stress-deflection curve
7. Efficiency factors visualization
8. 3D deflection profile

## References
This program implements theories for analyzing fiber reinforced concrete with steel reinforcement
This program implements theories based on the recrnt manuscript for analyzing fiber reinforced concrete with steel reinforcement.
[1] C. Pleesudjai, D. D. Patel, and B. Mobasher, “Generalized Nonlinear Moment–Curvature Model for Serviceability-Based Design of Hybrid Reinforced Concrete,” J. Struct. Eng., vol. 149, no. 12, Dec. 2023
doi: 10.1061/JSENDH.STENG-12235.
[2] Y. Yao, K. Aswani, X. Wang, and B. Mobasher, “Analytical displacement solutions for statically determinate beams based on a trilinear moment–curvature model,” Struct. Concr., vol. 19, no. 6, pp. 1619–1632, Dec. 2018, doi: 10.1002/suco.201700150.
[3] Y. Yao, X. Wang, K. Aswani, and B. Mobasher, “Analytical procedures for design of strain softening and hardening cement composites,” Int. J. Adv. Eng. Sci. Appl. Math., vol. 9, no. 3, pp. 181–194, Sep. 2017, doi: 10.1007/s12572-017-0187-4.
[4] B. Mobasher, Y. Yao, and C. Soranakom, “Analytical solutions for flexural design of hybrid steel fiber reinforced concrete beams,” Eng. Struct., vol. 100, pp. 164–177, Oct. 2015, doi: 10.1016/j.engstruct.2015.06.006.
