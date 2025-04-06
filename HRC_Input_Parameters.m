function [params] = HRC_Input_Parameters()
% HRC_Input_Parameters - Function to define all input parameters for HRC analysis
% Returns a structure containing all necessary parameters for the analysis
%
% Function creates and returns a structure with all the input parameters
% needed for the Hybrid Reinforcement Program for designing FRC with Rebar
%
% Returns:
%   params: Structure containing all input parameters

% Initialize the parameter structure
params = struct();

%##########################################################################
% Part A: Output File Configuration
%##########################################################################
params.progVer = 'MainProg_LD3P4P_HRC_Double_Reinforced_Chidchanok_11072024';  % Program version
params.fname   = 'Main_output.dat';           % Main output file
params.fname2  = 'Efficiency_output.dat';     % Efficiency factors output file
params.fname3  = 'Intersection_output.dat';   % Intersection points output file

% Experimental Tension Data
params.hasTens       = 0;                     % Flag for tension data availability (0 for No, 1 for Yes)
params.fnameTens     = '.txt';                % Tension data file name
params.startRowTens  = 1;                     % Starting row for data
params.xColTens      = 1;                     % X data column
params.yColTens      = 2;                     % Y data column
params.bk0tb1cm2Tens = 0;                     % Data delimiter (0 for space, 1 for tab, 2 for comma)

% Experimental Flexural Data
params.hasFlex       = 0;                     % Flag for flexural data availability
params.fnameFlex     = 'FHWA MCV.dat';        % Flexural data file name
params.startRowFlex  = 1;                     % Starting row for data
params.xColFlex      = 1;                     % X data column
params.yColFlex      = 2;                     % Y data column
params.bk0tb1cm2Flex = 0;                     % Data delimiter (0 for space, 1 for tab, 2 for comma)

% Testing Type and Localization Parameters
params.pointBend     = 4;                     % 3 for 3-point bending test, 4 for 4-point bending test
params.S2            = 100/3;                 % Middle spacing S2, S1=(L-S2)/2; S1+S2+S1 = L (for 4PB only)
params.Lp            = 100/3;                 % Plastic length for localized zone (required for 3PB, ignored in 4PB)
params.c             = 0.5;                   % Localized length/spacing ratio (≤ 0.5 for 4PB, ignored for 3PB)
params.unLoadFactor1 = 1;                     % Unload modulus for damaged non-localized zone (1st descending)
params.unLoadFactor2 = 1;                     % Unload modulus for damaged non-localized zone (2nd descending)

% Moment-Curvature Animation Settings
params.mmovie        = 0;                     % 1 for animation, 0 for no animation
params.storemovie    = 0;                     % Store movie

%##########################################################################
% Part B: Analysis Parameters
%##########################################################################
params.Type = 1; % Model Type (0 for FRC Model, 1 for HRC Model)

% Beam Geometry
params.b       = 12;                          % Width
params.h       = 24;                          % Total depth
params.L       = 100;                         % Clear span
params.alpha   = 21.795 / params.h;           % Depth of steel to total depth ratio
params.d       = params.h;                    % Effective depth

% Material Parameters
params.Beta_design = 30;
params.E           = 6933 * 1000;             % Tensile Young's Modulus in psi or MPa

% Tension Model Parameters
params.epsilon_cr  = 147 * 10^(-6);           % First-cracking strain
params.mu          = 0.5;                     % Normalized residual tensile strength
params.beta_tu     = 200;                     % Ultimate tensile strain
params.tau         = 20.4;                    % Transition zone normalized tensile strain
% Auto-calculated eta
params.eta         = (params.mu * params.epsilon_cr - params.epsilon_cr) / ((params.tau - 1) * params.epsilon_cr);

% Compression Model Parameters
params.omega       = 18.4;                    % Compressive yield strain parameter
params.xi          = 0.99;                    % Compressive Young's Modulus parameter
params.lambda_cu   = 20;                      % Ultimate compressive strain parameter
params.fc          = params.xi * params.E * params.epsilon_cr * params.omega;

% Steel Properties
params.rho         = 0.005;                   % Steel area ratio
params.kappa       = 14.1;                    % Rebar yield strain parameter
params.n           = 4.2;                     % Rebar Young's Modulus parameter
params.zeta        = 0.0001;                  % Compression/Tension Steel Area ratio

% Tolerance for Load-Deflection Algorithm
params.tor         = 10^-4;                   % Tolerance for cracking moment check

% Discretization Parameters
params.subMC       = [50, 50, 40];            % Subdivision count in compressive regions for M-C diagram
params.nSeg        = [40 * 2, 20 * 2];        % Segments count for spacing (L-Lp)/2 and Lp/2 for 3PB, L/3 and L/6 for 4PB

% Combine Compression and Tension Models for Plotting
params.strain      = [-params.lambda_cu, -params.omega, 0, 1, params.tau, params.tau, params.beta_tu] * params.epsilon_cr;
params.stress      = [-params.omega * params.xi, -params.omega * params.xi, 0, 1, params.eta * (params.tau - 1) + 1, params.mu, params.mu] * params.E * params.epsilon_cr;

params.strain_st   = [0, params.kappa, params.kappa * 100] * params.epsilon_cr;
params.stress_st   = [0, params.kappa * params.n, params.kappa * params.n] * params.E * params.epsilon_cr;

% Derived parameters - common calculations
params.EIcr = (1/12) * params.E * params.b * params.h^3;    % First cracking Stiffness
params.Mcr = (1/6) * params.epsilon_cr * params.E * params.b * params.h^2;  % First cracking Moment
params.phicr = 2 * params.epsilon_cr / params.h;            % First cracking Curvature

end