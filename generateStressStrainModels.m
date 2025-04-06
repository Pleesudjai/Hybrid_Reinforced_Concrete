function [models] = generateStressStrainModels(params)
% generateStressStrainModels - Creates stress-strain models for analysis
%
% This function generates the stress-strain relationships for concrete in
% compression and tension, as well as for steel reinforcement.
% It returns a structure containing all the model data points.
%
% Parameters:
%   params: Structure containing model parameters
%
% Returns:
%   models: Structure with stress-strain model data

% Initialize the output structure
models = struct();

% Generate Concrete Compression-Tension Model points
models.strain = [-params.lambda_cu, -params.omega, 0, 1, params.tau, params.tau, params.beta_tu] * params.epsilon_cr;
models.stress = [-params.omega * params.xi, -params.omega * params.xi, 0, 1, params.eta * (params.tau - 1) + 1, params.mu, params.mu] * params.E * params.epsilon_cr;

% Generate Steel Rebar Model points
models.strain_st = [0, params.kappa, params.kappa * 100] * params.epsilon_cr;
models.stress_st = [0, params.kappa * params.n, params.kappa * params.n] * params.E * params.epsilon_cr;

% Plot Compression Model
fig1 = figure(1);
plot(-models.strain(1:3), -models.stress(1:3), '-or', 'LineWidth', 1.5);
grid on;
title('Compression Model');
xlabel('Strain');
ylabel('Stress');
saveas(fig1, 'Compression_Model.png');

% Plot Tension Model
fig2 = figure(2);
plot(models.strain(3:7), models.stress(3:7), '-o', 'LineWidth', 1.5);
grid on;
title('Tension Model');
xlabel('Strain');
ylabel('Stress');
saveas(fig2, 'Tension_Model.png');

% Plot Steel Rebar Model
fig3 = figure(3);
plot(models.strain_st(1:3), models.stress_st(1:3), '-ok', 'LineWidth', 1.5);
grid on;
title('Rebar Model');
xlabel('Strain');
ylabel('Stress');
saveas(fig3, 'Rebar_Model.png');

end