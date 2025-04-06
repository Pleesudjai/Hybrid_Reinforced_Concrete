%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainProgram_HRC_Analysis.m
% Hybrid Reinforced Concrete (HRC) Analysis Program
%
% This program analyzes HRC beams with fiber reinforced concrete (FRC) 
% and steel reinforcement. It calculates moment-curvature relationships,
% efficiency factors, and load-deflection responses for 3-point or 4-point
% bending tests.
%
% Features:
% - Trilinear tension model and bilinear compression model for concrete
% - Steel reinforcement modeling
% - Moment-curvature analysis with efficiency factors
% - Load-deflection prediction for 3-point and 4-point bending
% - Visualization of results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

%##########################################################################
% Load all input parameters from the separate function
%##########################################################################
params = HRC_Input_Parameters();

% Extract key parameters to local variables for cleaner code
b = params.b;
h = params.h;
L = params.L;
alpha = params.alpha;
d = params.d;
E = params.E;
epsilon_cr = params.epsilon_cr;
mu = params.mu;
beta_tu = params.beta_tu;
tau = params.tau;
eta = params.eta;
omega = params.omega;
xi = params.xi;
lambda_cu = params.lambda_cu;
fc = params.fc;
rho = params.rho;
kappa = params.kappa;
n = params.n;
zeta = params.zeta;
tor = params.tor;
subMC = params.subMC;
nSeg = params.nSeg;
strain = params.strain;
stress = params.stress;
strain_st = params.strain_st;
stress_st = params.stress_st;
Type = params.Type;
pointBend = params.pointBend;
S2 = params.S2;
Lp = params.Lp;
c = params.c;
unLoadFactor1 = params.unLoadFactor1;
unLoadFactor2 = params.unLoadFactor2;
mmovie = params.mmovie;
storemovie = params.storemovie;
hasTens = params.hasTens;
fnameTens = params.fnameTens;
startRowTens = params.startRowTens;
xColTens = params.xColTens;
yColTens = params.yColTens;
bk0tb1cm2Tens = params.bk0tb1cm2Tens;
hasFlex = params.hasFlex;
fnameFlex = params.fnameFlex;
startRowFlex = params.startRowFlex;
xColFlex = params.xColFlex;
yColFlex = params.yColFlex;
bk0tb1cm2Flex = params.bk0tb1cm2Flex;
fname = params.fname;
fname2 = params.fname2;
fname3 = params.fname3;
progVer = params.progVer;

% Get derived parameters
EIcr = params.EIcr;
Mcr = params.Mcr;
phicr = params.phicr;

%##########################################################################
% Part B: Define the material stress-strain model
%##########################################################################
% Define stress-strain points for plotting
strain      = [-lambda_cu, -omega, 0, 1, tau, tau, beta_tu] * epsilon_cr;
stress      = [-omega * xi, -omega * xi, 0, 1, eta * (tau - 1) + 1, mu, mu] * E * epsilon_cr;

strain_st   = [0, kappa, kappa * 100] * epsilon_cr;
stress_st   = [0, kappa * n, kappa * n] * E * epsilon_cr;

% Plot compression model
fig1 = figure(1);
plot(-strain(1:3), -stress(1:3), '-or', 'LineWidth', 1.5);
grid on;
title('Compression Model');
xlabel('Strain');
ylabel('Stress');
saveas(fig1, 'Compression_Model.png');

% Plot tension model
fig2 = figure(2);
plot(strain(3:7), stress(3:7), '-o', 'LineWidth', 1.5);
grid on;
title('Tension Model');
xlabel('Strain');
ylabel('Stress');
saveas(fig2, 'Tension_Model.png');

% Plot steel rebar model
fig3 = figure(3);
plot(strain_st(1:3), stress_st(1:3), '-ok', 'LineWidth', 1.5);
grid on;
title('Rebar Model');
xlabel('Strain');
ylabel('Stress');
saveas(fig3, 'Rebar_Model.png');

%##########################################################################
% Part C: Calculate Moment-Curvature relationship
%##########################################################################
% Calculate cracking values
EIcr = (1/12) * E * b * h^3;             % First cracking stiffness
Mcr = (1/6) * epsilon_cr * E * b * h^2;  % First cracking moment
phicr = 2 * epsilon_cr / h;              % First cracking curvature

% Calculate moment-curvature points for each zone
% Zone 1: Elastic concrete and elastic steel
for i = 1:(10*beta_tu)
    [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency_1] = zone1_2024(i,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha); 
    output1(i,1) = beta;                  % Beta with increment by 0.01
    output1(i,2) = lamb;                  % Lambda corresponding to beta
    output1(i,3) = kap;                   % Kappa at bottom steel corresponding to beta
    output1(i,4) = k;                     % Parameter of neutral axis
    output1(i,5) = nM;                    % Normalized Moment (nM = M/Mcr)
    output1(i,6) = M;                     % Actual Moment
    output1(i,7) = nphi;                  % Normalized Curvature (nphi = phi/phicr)
    output1(i,8) = phi;                   % Actual Curvature                 
    output1(i,9) = nstiff;                % Normalized Stiffness (nstiff = stiff/EIcr)
    output1(i,10) = stiff;                % Actual Stiffness
    output1(i,11) = netf;                 % Net force along section
    output1(i,12) = kap_com;              % Kappa at top steel corresponding to beta
    output1(i,13) = eps_bot;              % Bottom Strain
    output1(i,14) = eps_top;              % Top Strain
    output1(i,15) = eps_st_bot;           % Tension rebar strain
    output1(i,16) = eps_st_top;           % Compression rebar strain
    output1(i,17) = Matrix_sum_com;       % Summary matrix compression force
    output1(i,18) = Matrix_sum_ten;       % Summary matrix tension force
    output1(i,19) = Rebar_com;            % Rebar compression force
    output1(i,20) = Rebar_ten;            % Rebar tension force

    output1(i,21) = Efficiency_1(1);      % Efficiency factor: Concrete Compression
    output1(i,22) = Efficiency_1(2);      % Efficiency factor: Rebar Compression
    output1(i,23) = Efficiency_1(3);      % Efficiency factor: Concrete Tension
    output1(i,24) = Efficiency_1(4);      % Efficiency factor: Rebar Tension
end

% Zone 2.1: Elastic concrete in compression, concrete in tension follows 2nd slope, elastic steel
for i = 1:(10*beta_tu)
    [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency_21] = zone21_2024(i,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta); 
    output21(i,1) = beta;
    output21(i,2) = lamb;
    output21(i,3) = kap;
    output21(i,4) = k;
    output21(i,5) = nM;
    output21(i,6) = M;
    output21(i,7) = nphi;
    output21(i,8) = phi;
    output21(i,9) = nstiff; 
    output21(i,10) = stiff; 
    output21(i,11) = netf; 
    output21(i,12) = kap_com; 
    output21(i,13) = eps_bot; 
    output21(i,14) = eps_top; 
    output21(i,15) = eps_st_bot; 
    output21(i,16) = eps_st_top;
    output21(i,17) = Matrix_sum_com;
    output21(i,18) = Matrix_sum_ten;
    output21(i,19) = Rebar_com;
    output21(i,20) = Rebar_ten;

    output21(i,21) = Efficiency_21(1);    % Efficiency factor: Concrete Compression
    output21(i,22) = Efficiency_21(2);    % Efficiency factor: Rebar Compression
    output21(i,23) = Efficiency_21(3);    % Efficiency factor: Concrete Tension
    output21(i,24) = Efficiency_21(4);    % Efficiency factor: Rebar Tension
end

% Zone 2.2: Elastic concrete in compression, concrete in tension follows 2nd slope, yielded steel
for i = 1:(10*beta_tu)
    [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency_22] = zone22_2024(i,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa); 
    output22(i,1) = beta;
    output22(i,2) = lamb;
    output22(i,3) = kap;
    output22(i,4) = k;
    output22(i,5) = nM;
    output22(i,6) = M;
    output22(i,7) = nphi;
    output22(i,8) = phi;
    output22(i,9) = nstiff; 
    output22(i,10) = stiff; 
    output22(i,11) = netf; 
    output22(i,12) = kap_com; 
    output22(i,13) = eps_bot; 
    output22(i,14) = eps_top; 
    output22(i,15) = eps_st_bot; 
    output22(i,16) = eps_st_top;
    output22(i,17) = Matrix_sum_com;
    output22(i,18) = Matrix_sum_ten;
    output22(i,19) = Rebar_com;
    output22(i,20) = Rebar_ten;

    output22(i,21) = Efficiency_22(1);    % Efficiency factor: Concrete Compression
    output22(i,22) = Efficiency_22(2);    % Efficiency factor: Rebar Compression
    output22(i,23) = Efficiency_22(3);    % Efficiency factor: Concrete Tension
    output22(i,24) = Efficiency_22(4);    % Efficiency factor: Rebar Tension
end

% Zone 3.1: Plastic concrete in compression, concrete in tension follows 2nd slope, elastic steel
for i = 1:(10*beta_tu)
    [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency_31] = zone31_2024(i,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa,omega);
    output31(i,1) = beta;
    output31(i,2) = lamb;
    output31(i,3) = kap;
    output31(i,4) = k;
    output31(i,5) = nM;
    output31(i,6) = M;
    output31(i,7) = nphi;
    output31(i,8) = phi;
    output31(i,9) = nstiff; 
    output31(i,10) = stiff; 
    output31(i,11) = netf; 
    output31(i,12) = kap_com; 
    output31(i,13) = eps_bot; 
    output31(i,14) = eps_top; 
    output31(i,15) = eps_st_bot; 
    output31(i,16) = eps_st_top;
    output31(i,17) = Matrix_sum_com;
    output31(i,18) = Matrix_sum_ten;
    output31(i,19) = Rebar_com;
    output31(i,20) = Rebar_ten;

    output31(i,21) = Efficiency_31(1);    % Efficiency factor: Concrete Compression
    output31(i,22) = Efficiency_31(2);    % Efficiency factor: Rebar Compression
    output31(i,23) = Efficiency_31(3);    % Efficiency factor: Concrete Tension
    output31(i,24) = Efficiency_31(4);    % Efficiency factor: Rebar Tension
end

% Zone 3.2: Plastic concrete in compression, concrete in tension follows 2nd slope, yielded steel
for i = 1:(10*beta_tu)
    [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency_32] = zone32_2024(i,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa,omega);
    output32(i,1) = beta;
    output32(i,2) = lamb;
    output32(i,3) = kap;
    output32(i,4) = k;
    output32(i,5) = nM;
    output32(i,6) = M;
    output32(i,7) = nphi;
    output32(i,8) = phi;
    output32(i,9) = nstiff; 
    output32(i,10) = stiff; 
    output32(i,11) = netf; 
    output32(i,12) = kap_com; 
    output32(i,13) = eps_bot; 
    output32(i,14) = eps_top; 
    output32(i,15) = eps_st_bot; 
    output32(i,16) = eps_st_top;
    output32(i,17) = Matrix_sum_com;
    output32(i,18) = Matrix_sum_ten;
    output32(i,19) = Rebar_com;
    output32(i,20) = Rebar_ten;

    output32(i,21) = Efficiency_32(1);    % Efficiency factor: Concrete Compression
    output32(i,22) = Efficiency_32(2);    % Efficiency factor: Rebar Compression
    output32(i,23) = Efficiency_32(3);    % Efficiency factor: Concrete Tension
    output32(i,24) = Efficiency_32(4);    % Efficiency factor: Rebar Tension
end

% Zone 4.1: Elastic concrete in compression, concrete in tension follows 3rd slope, elastic steel
for i = 1:(10*beta_tu)
    [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency_41] = zone41_2024(i,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa,omega,tau,mu);
    output41(i,1) = beta;
    output41(i,2) = lamb;
    output41(i,3) = kap;
    output41(i,4) = k;
    output41(i,5) = nM;
    output41(i,6) = M;
    output41(i,7) = nphi;
    output41(i,8) = phi;
    output41(i,9) = nstiff; 
    output41(i,10) = stiff; 
    output41(i,11) = netf; 
    output41(i,12) = kap_com; 
    output41(i,13) = eps_bot; 
    output41(i,14) = eps_top; 
    output41(i,15) = eps_st_bot; 
    output41(i,16) = eps_st_top;
    output41(i,17) = Matrix_sum_com;
    output41(i,18) = Matrix_sum_ten;
    output41(i,19) = Rebar_com;
    output41(i,20) = Rebar_ten;

    output41(i,21) = Efficiency_41(1);    % Efficiency factor: Concrete Compression
    output41(i,22) = Efficiency_41(2);    % Efficiency factor: Rebar Compression
    output41(i,23) = Efficiency_41(3);    % Efficiency factor: Concrete Tension
    output41(i,24) = Efficiency_41(4);    % Efficiency factor: Rebar Tension
end

% Zone 4.2: Elastic concrete in compression, concrete in tension follows 3rd slope, yielded steel
for i = 1:(10*beta_tu)
    [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency_42] = zone42_2024(i,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa,omega,tau,mu);
    output42(i,1) = beta;
    output42(i,2) = lamb;
    output42(i,3) = kap;
    output42(i,4) = k;
    output42(i,5) = nM;
    output42(i,6) = M;
    output42(i,7) = nphi;
    output42(i,8) = phi;
    output42(i,9) = nstiff; 
    output42(i,10) = stiff; 
    output42(i,11) = netf; 
    output42(i,12) = kap_com; 
    output42(i,13) = eps_bot; 
    output42(i,14) = eps_top; 
    output42(i,15) = eps_st_bot; 
    output42(i,16) = eps_st_top;
    output42(i,17) = Matrix_sum_com;
    output42(i,18) = Matrix_sum_ten;
    output42(i,19) = Rebar_com;
    output42(i,20) = Rebar_ten;

    output42(i,21) = Efficiency_42(1);    % Efficiency factor: Concrete Compression
    output42(i,22) = Efficiency_42(2);    % Efficiency factor: Rebar Compression
    output42(i,23) = Efficiency_42(3);    % Efficiency factor: Concrete Tension
    output42(i,24) = Efficiency_42(4);    % Efficiency factor: Rebar Tension
end

% Zone 5.1: Plastic concrete in compression, concrete in tension follows 3rd slope, elastic steel
for i = 1:(10*beta_tu)
    [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency_51] = zone51_2024(i,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa,omega,tau,mu);
    output51(i,1) = beta;
    output51(i,2) = lamb;
    output51(i,3) = kap;
    output51(i,4) = k;
    output51(i,5) = nM;
    output51(i,6) = M;
    output51(i,7) = nphi;
    output51(i,8) = phi;
    output51(i,9) = nstiff; 
    output51(i,10) = stiff; 
    output51(i,11) = netf; 
    output51(i,12) = kap_com; 
    output51(i,13) = eps_bot; 
    output51(i,14) = eps_top; 
    output51(i,15) = eps_st_bot; 
    output51(i,16) = eps_st_top;
    output51(i,17) = Matrix_sum_com;
    output51(i,18) = Matrix_sum_ten;
    output51(i,19) = Rebar_com;
    output51(i,20) = Rebar_ten;

    output51(i,21) = Efficiency_51(1);    % Efficiency factor: Concrete Compression
    output51(i,22) = Efficiency_51(2);    % Efficiency factor: Rebar Compression
    output51(i,23) = Efficiency_51(3);    % Efficiency factor: Concrete Tension
    output51(i,24) = Efficiency_51(4);    % Efficiency factor: Rebar Tension
end

% Zone 5.2: Plastic concrete in compression, concrete in tension follows 3rd slope, yielded steel
for i = 1:(10*beta_tu)
    [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency_52] = zone52_2024(i,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa,omega,tau,mu);
    output52(i,1) = beta;
    output52(i,2) = lamb;
    output52(i,3) = kap;
    output52(i,4) = k;
    output52(i,5) = nM;
    output52(i,6) = M;
    output52(i,7) = nphi;
    output52(i,8) = phi;
    output52(i,9) = nstiff; 
    output52(i,10) = stiff; 
    output52(i,11) = netf; 
    output52(i,12) = kap_com;
    output52(i,13) = eps_bot; 
    output52(i,14) = eps_top; 
    output52(i,15) = eps_st_bot; 
    output52(i,16) = eps_st_top;
    output52(i,17) = Matrix_sum_com;
    output52(i,18) = Matrix_sum_ten;
    output52(i,19) = Rebar_com;
    output52(i,20) = Rebar_ten;

    output52(i,21) = Efficiency_52(1);    % Efficiency factor: Concrete Compression
    output52(i,22) = Efficiency_52(2);    % Efficiency factor: Rebar Compression
    output52(i,23) = Efficiency_52(3);    % Efficiency factor: Concrete Tension
    output52(i,24) = Efficiency_52(4);    % Efficiency factor: Rebar Tension
end

%##########################################################################
% Part D: Generate envelope curve - select valid route for analysis
%##########################################################################
% Possible routes through the zones:
% (A) Stage 1 >> Stage 2.1 >> Stage 2.2 >> Stage 3.2 >> Stage 5.2
% (B) Stage 1 >> Stage 2.1 >> Stage 2.2 >> Stage 4.2 >> Stage 5.2
% (C) Stage 1 >> Stage 2.1 >> Stage 3.1 >> Stage 3.2 >> Stage 5.2
% (D) Stage 1 >> Stage 2.1 >> Stage 3.1 >> Stage 5.1 >> Stage 5.2
% (E) Stage 1 >> Stage 2.1 >> Stage 4.1 >> Stage 4.2 >> Stage 5.2 
% (F) Stage 1 >> Stage 2.1 >> Stage 4.1 >> Stage 5.1 >> Stage 5.2 

for j=1:24
    
    envelope(1:10,j) = output1(1:10,j);
    %------------------------------------------------------------------------
    % For Route (A) and (B)
    %------------------------------------------------------------------------
    if output21(tau*10,3) > kappa
        
       BB = find(output21(:,3) < kappa); 
       AA = find(output22(:,2) < omega); 
       
       % For (A) Stage 1 >> Stage 2.1 >> Stage 2.2 >> Stage 3.2 >> Stage 5.2 
        if length(AA(:,1)) < tau*10
            int22 = length(BB(:,1));
            envelope(11:int22,j) = output21(11:int22,j);
            dum32 = find(output22(:,2) < omega);
            int32 = length(dum32(:,1));
            envelope(int22+1:int32,j) = output22(int22+1:int32,j);
            envelope(int32+1:tau*10,j) = output32(int32+1:tau*10,j);
            envelope(tau*10+1:beta_tu*10,j) = output52(tau*10+1:beta_tu*10,j);
            
        % For (B) Stage 1 >> Stage 2.1 >> Stage 2.2 >> Stage 4.2 >> Stage 5.2
        else
            int22 =  length(BB(:,1));
            envelope(11:int22,j) = output21(11:int22,j);
            envelope(int22+1:tau*10,j) = output22(int22+1:tau*10,j);
            dum52 = find(output42(:,2) < omega);
            int52 = length(dum52(:,1));
            envelope(tau*10+1:int52,j) = output42(tau*10+1:int52,j);
            envelope(int52+1:beta_tu*10,j) = output52(int52+1:beta_tu*10,j);
       end
    
    %------------------------------------------------------------------------
    % For Route (C),(D),(E) and (F)
    %------------------------------------------------------------------------
    else

        % For Route (E) and (F)
        if output21(tau*10,2) < omega % When tensile experiences residual stress, compression is still in Elastic
            envelope(11:(tau*10),j) = output21(11:(tau*10),j);
            a = find(output41(:,2) < omega);
            B = find(output41(:,3) < kappa);
            lengtha = length(a(:,1));
            lengthb = length(B(:,1));

            if a(lengtha,1) < B(lengthb,1) % For Route (F)
                int41 = a(lengtha,1);
                C = find(output51(:,3) < kappa);
                lengthc = length(C(:,1));
                int51 = C(lengthc,1);
                envelope((tau*10):int41,j) = output41((tau*10):int41,j);
                envelope((int41+1):int51,j) = output51((int41+1):int51,j);
                envelope((int51+1):beta_tu*10,j) = output52((int51+1):beta_tu*10,j);

            else % For Route (E)
                int41 = B(lengthb,1);
                C = find(output42(:,2) < omega);
                lengthc = length(C(:,1));
                int42 = C(lengthc,1);
                envelope((tau*10):int41,j) = output41((tau*10):int41,j);
                envelope((int41+1):int42,j) = output42((int41+1):int42,j);
                envelope((int42+1):beta_tu*10,j) = output52((int42+1):beta_tu*10,j);
            end

        %------------------------------------------------------------------------
        % For Route (C) and (D)
        %------------------------------------------------------------------------
        else
            AA = find(output31(:,3) < kappa);
            lengthAA = length(AA(:,1));
            
            BB = find(output21(:,2) < omega);
            lengthBB = length(BB(:,1));

            % For Route (D) Stage 1 >> Stage 2.1 >> Stage 3.1 >> Stage 5.1 >> Stage 5.2 
            if tau*10 < AA(lengthAA,1)  % Tensile is in phase 3 before rebar yield
                int31 = BB(lengthBB,1);
                envelope(11:int31,j) = output21(11:int31,j);
                envelope(int31+1:tau*10,j) = output31(int31+1:tau*10,j);
                dum51 = find(output51(:,3) < kappa);
                int51 = length(dum51(:,1));
                envelope(tau*10+1:int51,j) = output51(tau*10+1:int51,j);
                envelope(int51+1:beta_tu*10,j) = output52(int51+1:beta_tu*10,j);
                
            % For Route (C) Stage 1 >> Stage 2.1 >> Stage 3.1 >> Stage 3.2 >> Stage 5.2 
            else 
                int31 = BB(lengthBB,1);
                envelope(11:int31,j) = output21(11:int31,j);
                dum32 = find(output31(:,3) < kappa);
                int32 = length(dum32(:,1));
                envelope(int31+1:int32,j) = output31(int31+1:int32,j);
                envelope(int32+1:tau*10,j) = output32(int32+1:tau*10,j);
                envelope(tau*10+1:beta_tu*10,j) = output52(tau*10+1:beta_tu*10,j);
            end
        end     
    end
end

% Calculate force contributions
Matrix_sum_com = envelope(1:end,17);
Matrix_sum_ten = envelope(1:end,18);
Rebar_com = envelope(1:end,19);
Rebar_ten = envelope(1:end,20);
Total_force = Matrix_sum_com + Matrix_sum_ten + Rebar_com + Rebar_ten;

Matrix_Residual_ten_over_total = Matrix_sum_ten./Total_force.*2;
Matrix_com_over_total = Matrix_sum_com./Total_force.*2;
Rebar_ten_over_total = Rebar_ten./Total_force.*2;
Rebar_com_over_total = Rebar_com./Total_force.*2;

%##########################################################################
% Part E: Load-Deflection Analysis
%##########################################################################
% Get curvature and moment from envelope
CV = envelope(:,8);                       % Curvature vector
MM = envelope(:,6);                       % Moment vector
EqualFlexStress = (envelope(:,6).*h/2)/(b*(h^3)/12);  % Equivalent flexural stress
EI = (MM(10)-MM(3))/(CV(10)-CV(3));       % Initial stiffness (before cracking)
[Mmax, rowMax] = max(MM);                 % Maximum moment
Cmax = CV(rowMax);                        % Curvature at maximum moment

% For debugging and analysis
k_fuction_phi(:,1) = envelope(:,8);       % Phi value, x axis
k_fuction_phi(1,1) = 0;                   % Zero initial value
k_fuction_phi(:,2) = envelope(:,4);       % k value, y axis
k_fuction_beta(:,1) = envelope(:,1);      % Beta, x-axis
k_fuction_beta(1,1) = 0;                  % Zero initial value
k_fuction_beta(:,2) = envelope(:,4);      % k value, y axis

% Break moment-curvature into segments for different phases
CIncCurve1 = [0, 0.5*10^-7, 1*10^-7];     % First ascending curve (C points)
MIncCurve1 = [0, 0.5, 1];                 % First ascending curve (M points)
CDecCurve1 = [0, 0.5*10^-7, 1*10^-7];     % First descending curve
MDecCurve1 = [0, 0.5, 1];                 % First descending curve
CIncCurve2 = [0, 0.5*10^-7, 1*10^-7];     % Second ascending curve
MIncCurve2 = [0, 0.5, 1];                 % Second ascending curve
CDecCurve2 = [0, 0.5*10^-7, 1*10^-7];     % Second descending curve
MDecCurve2 = [0, 0.5, 1];                 % Second descending curve

% Initialize
CIncCurve1(1) = CV(1);
MIncCurve1(1) = MM(1);
segment = 1;
N1 = 1;

% Identify segments in moment-curvature curve
for i = 2:length(MM)
    % Segment 1: First ascending
    if MM(i) - MM(i-1) > 0 && segment == 1
        N1 = N1 + 1;
        CIncCurve1(N1) = CV(i);
        CIncCurve1(1:1) = 0;
        MIncCurve1(N1) = MM(i); 
        MIncCurve1(1:1) = 0;
    elseif segment == 1
        segment = 2;
        N2 = 1;
        CDecCurve1(1) = CV(i-1);
        MDecCurve1(1) = MM(i-1);
    end
    
    % Segment 2: First descending
    if MM(i) - MM(i-1) < 0 && segment == 2
        N2 = N2 + 1;
        CDecCurve1(N2) = CV(i);
        MDecCurve1(N2) = MM(i);
    elseif segment == 2
        segment = 3;
        N3 = 1;
        CIncCurve2(1) = CV(i-1);
        MIncCurve2(1) = MM(i-1);
    end
    
    % Segment 3: Second ascending
    if MM(i) - MM(i-1) > 0 && segment == 3
        N3 = N3 + 1;
        CIncCurve2(N3) = CV(i);
        MIncCurve2(N3) = MM(i);
    elseif segment == 3
        segment = 4;
        N4 = 1;
        CDecCurve2(1) = CV(i);
        MDecCurve2(1) = MM(i);
    end
    
    % Segment 4: Second descending
    if MM(i) - MM(i-1) < 0 && segment == 4
        N4 = N4 + 1;
        CDecCurve2(N4) = CV(i);
        MDecCurve2(N4) = MM(i);
    elseif segment == 4
        segment = 4;
        N4 = N4 + 1;
        CDecCurve2(N4) = CV(i);
        MDecCurve2(N4) = MM(i);
    end 
end

% Plot moment-curvature diagram with segments
fig7 = figure(7);
plot(CV(1:end,1), MM(1:end,1), '--ro', CIncCurve1, MIncCurve1, 'b.', CDecCurve1, MDecCurve1, 'g.', CIncCurve2, MIncCurve2, 'k.', CDecCurve2, MDecCurve2, 'm.'), grid, 
title(['Moment Curvature Diagram at \mu = ', num2str(mu)]);
xlabel('Curvature'), ylabel('Moment'),
legend('MC', '1st ascending', 'descending', '2nd ascending', '2nd descending')
saveas(fig7, 'Moment_Curvature_Diagram.png');

% Calculate curvature distribution along the beam
if pointBend == 3
    S = L/2;                            % Equal spacing for three-point bending
    S1 = (L-Lp)/2;                      % First spacing (support to first load P/2)
    X1 = linspace(0, S1, nSeg(1));
    X2 = linspace(S1, L/2, nSeg(2)+1);
    X = [X1(1:end), X2(2:end)];
elseif pointBend == 4
    S = S2;               
    S1 = (L-S2)/2;                      % First spacing (support to first load P/2)
    nSeg22 = round(2*c*nSeg(2));        % Number of segments for localized zone
    nSeg21 = nSeg(2) - nSeg22;          % Number of segments for non-localized zone
    X1 = linspace(0, (L-S2)/2, nSeg(1));
    X21 = linspace(S1, S1+((S/2)-(2*c*S/2)), nSeg21+1);
    X22 = linspace(S1+((S/2)-(2*c*S/2)), L/2, nSeg22+1);
    if c == 0
        X = [X1(1:end), X21(2:end)];
    elseif c == 0.5
        X = [X1(1:end), X22(2:end)];
    else
        X = [X1(1:end), X21(2:end), X22(2:end)];
    end
end

% Calculate reaction force at each load step
if pointBend == 3
    R = (MM/S);                          % Reaction force for 3PB
    Rs = (MM/S1);                        % Support reaction force
elseif pointBend == 4
    R = MM/S1;                           % Reaction force for 4PB
    Rs = (MM/S1);                        % Support reaction force
end

% Calculate flexural stress
flexStrs = (MM)*6/(b*h^2);              % Equivalent flexural stress

% Define key load steps
stepY = 11;                              % Yield point (intersection between zone 1 and 2)
[maxM, stepM] = max(MM);                 % Maximum moment
stepF = length(MM);                      % Last step (failure)

% Calculate deflection along beam
totSeg = length(X);
totStep = length(R);
Mmt = zeros(totStep, totSeg);            % Moment distribution
Phi = zeros(totStep, totSeg);            % Curvature distribution
damageX = zeros(totStep, totSeg);        % Damage location
Region = zeros(totStep, totSeg);         % Region type
pass = zeros(totStep, totSeg);           % Pass flag
delta = zeros(totStep, totSeg);          % Deflection
segment = zeros(totStep, totSeg);        % Segment type
rotation = zeros(totStep, totSeg);       % Rotation
segment(1,:) = 1;
segment(:,1) = 1;

% Calculate deflection for each load step
for i = 2:totStep
    delta(i) = 0;
    for j = 2:totSeg
        % Determine moment at each location
        if pointBend == 3 && i <= stepY    % Elastic zone (before damage)
            Mmt(i,j) = R(i)*X(j);          
            Region(i,j) = 1;              % Non-localized damage zone  
        else
            % Non-localized damage zone (3PB and 4PB, outside plastic length)
            if X(j) <= S1                    
                Mmt(i,j) = Rs(i)*(X(j));
                Region(i,j) = 1;          
            % Localized damage zone (3PB, inside plastic length) 
            elseif (X(j) > S1) && pointBend == 3
                Mmt(i,j) = Rs(i)*S1;
                Region(i,j) = 2;          % Localized damage zone
            % Non-localized damage zone (4PB, X > S1)
            elseif (X(j) > S1 && X(j) <= (S1+((S/2)-(2*c*S/2)))) && pointBend == 4 
                Mmt(i,j) = Rs(i)*S1;
                Region(i,j) = 1;        
            % Localized damage zone (4PB, X > S1)
            elseif (X(j) > (S1+((S/2)-(2*c*S/2)))) && (X(j) <= L/2) && pointBend == 4
                Mmt(i,j) = Rs(i)*S1;
                Region(i,j) = 2;        
            end
        end
        
        % Determine segment type based on loading history
        % Segment 1: First ascending
        if segment(i-1,j) == 1 && Mmt(i,j) - Mmt(i-1,j) > 0
            segment(i,j) = 1;
        elseif segment(i-1,j) == 1
            segment(i,j) = 2;
        end

        % Segment 2: First descending
        if segment(i-1,j) == 2 && Mmt(i,j) - Mmt(i-1,j) < 0
            segment(i,j) = 2;
        elseif segment(i-1,j) == 2 
            segment(i,j) = 3;
        end

        % Segment 3: Second ascending
        if segment(i-1,j) == 3 && Mmt(i,j) - Mmt(i-1,j) > 0
            segment(i,j) = 3;
        elseif segment(i-1,j) == 3
            segment(i,j) = 4;
        end

        % Segment 4: Second descending
        if segment(i-1,j) == 4 && Mmt(i,j) - Mmt(i-1,j) < 0
            segment(i,j) = 4;
        elseif segment(i-1,j) == 4
            segment(i,j) = 4;
        end
        
        % Calculate curvature based on region and segment type
        if Region(i,j) == 1  % Non-localized zone                   
            if segment(i,j) == 1
                Phi(i,j) = interp1(MIncCurve1, CIncCurve1, Mmt(i,j));
                
            elseif segment(i,j) == 2
                Phi(i,j) = Phi(i-1,j) - unLoadFactor1*(Mmt(i-1,j) - Mmt(i,j))/EI;  
                
            elseif segment(i,j) == 3  
                Phi(i,j) = interp1(MIncCurve2, CIncCurve2, Mmt(i,j));   
                % Handle interpolation issues
                if isnan(Phi(i,j))
                   slope1 = (Mmt(i-2,j) - Mmt(i-1,j))/(Phi(i-2,j)-Phi(i-1,j));
                   Phi(i,j) = Phi(i-1,j) + (Mmt(i-2,j) - Mmt(i-1,j))/slope1; 
                end  
                
            elseif segment(i,j) == 4 
                Phi(i,j) = Phi(i-1,j) - unLoadFactor2*(Mmt(i-1,j) - Mmt(i,j))/EI; 
            end
            
        elseif Region(i,j) == 2  % Localized zone
            if segment(i,j) == 1 
                Phi(i,j) = interp1(MIncCurve1, CIncCurve1, Mmt(i,j));
            elseif segment(i,j) == 2 
                Phi(i,j) = interp1(MDecCurve1, CDecCurve1, Mmt(i,j));
            elseif segment(i,j) == 3 
                Phi(i,j) = interp1(MIncCurve2, CIncCurve2, Mmt(i,j));  
                % Handle interpolation issues
                if isnan(Phi(i,j))
                  slope2 = (Mmt(i-2,j) - Mmt(i-1,j))/(Phi(i-2,j)-Phi(i-1,j));
                  Phi(i,j) = Phi(i-1,j) + (Mmt(i-2,j) - Mmt(i-1,j))/slope2; 
                end        
            elseif segment(i,j) == 4  
                Phi(i,j) = interp1(MDecCurve2, CDecCurve2, Mmt(i,j));
                % Handle interpolation issues
                if isnan(Phi(i,j))
                  slope3 = (Mmt(i-2,j) - Mmt(i-1,j))/(Phi(i-2,j)-Phi(i-1,j));
                  Phi(i,j) = Phi(i-1,j) + (Mmt(i-2,j) - Mmt(i-1,j))/slope3; 
                end        
            end
        end
        
        % Calculate rotation area and accumulate rotation
        area(i,j) = (0.5*(X(j) - X(j-1))*(Phi(i,j) + Phi(i,j-1)));
        rotation(i,j) = rotation(i,j-1) + (0.5*(X(j) - X(j-1))*(Phi(i,j) + Phi(i,j-1)));
        
        % Calculate beta distribution for shear
        Beta(i,j) = (Phi(i,j)*(1-interp1(k_fuction_phi(:,1), k_fuction_phi(:,2), Phi(i,j)))*d)/epsilon_cr;
    end
end

% Apply boundary conditions
rotation_leftend(:,1) = -1 * rotation(:,end);
for i = 1:totStep
    rotation(i,:) = rotation(i,:) + rotation_leftend(i,1);
end

% Calculate deflection profile
delta = zeros(totStep, totSeg);
delta_full_length = zeros(totStep, 2*totSeg-1);
area = zeros(totStep, totSeg);
dist = zeros(totStep, totSeg);
tangent = zeros(totStep, totSeg);

for i = 2:totStep
    delta(i,j) = 0;
    for j = 2:totSeg
        delta(i,j) = delta(i,j-1) + (0.5*(X(j) - X(j-1))*(rotation(i,j) + rotation(i,j-1)));
    end
end

% Convert to positive values for plotting
delta(:,:) = -1 * delta(:,:);

% Create full-length deflection profile
delta_full_length(:,1:totSeg) = delta(:,:);
for seg = 1:totSeg-1
    delta_full_length(:,totSeg+seg) = delta(:,totSeg-seg);
    XX(1,1:totSeg) = X(1,1:totSeg);
    XX(1,totSeg+seg) = ((L/2)-X(1,totSeg-seg))+L/2;
end

t_delta_full_length = delta_full_length';
t_XX = XX';

for nRow = 1:totStep
    location(nRow,:) = XX(1,:);
end

for nRow = 1:length(delta_full_length(1,:))
    loadstep(:,nRow) = 2*R(:,1);
end

t_loadstep = loadstep';

% Create 3D surface plot of deflection profile
adjust = beta_tu*10;
fig105 = figure(105);
surf(location(1:adjust,:), loadstep(1:adjust,:), delta_full_length(1:adjust,:), 'EdgeColor', 'none');
hold on;
ylabel('Load, N'); zlabel('Deflection, mm');
zlim([-0.5 30]);
ylim([0 200000]);
title('Deflection Profile');
saveas(fig105, 'Deflection_Profile_3D.png');

% Extract strain and stress at key load steps
strain_stepY = [envelope(stepY,14), envelope(stepY,16), envelope(stepY,15), envelope(stepY,13)] * epsilon_cr;
stress_stepY = [envelope(stepY,2)*xi, envelope(stepY,12)*n, envelope(stepY,3)*n, envelope(stepY,1)*eta*(tau-1)+1] * E * epsilon_cr;

strain_stepF = [envelope(stepF,14), envelope(stepF,16), envelope(stepF,15), envelope(stepF,13)] * epsilon_cr;
stress_stepF = [envelope(stepF,2)*xi, envelope(stepF,12)*n, envelope(stepF,3)*n, envelope(stepF,1)*mu] * E * epsilon_cr;

% Initialize arrays for experimental data
expTensStrn = [0, 0];
expTensStrs = [0, 0];
expFlexDisp = [0, 0];
expFlexLoad = [0, 0];
expFlexStrs = [0, 0];

% Load experimental data if available
if hasTens == 1
    % Load tension test data
    if bk0tb1cm2Tens == 0
        tensData = dlmread(fnameTens, '', startRowTens-1, 0);
    elseif bk0tb1cm2Tens == 1
        tensData = dlmread(fnameTens, '\t', startRowTens-1, 0);
    else
        tensData = dlmread(fnameTens, ',', startRowTens-1, 0);
    end
    % Read data up to maximum tensile strain
    [xMax, lastIndexTens] = max(tensData(:, xColTens));
    expTensStrn = tensData(1:lastIndexTens, xColTens);
    expTensStrs = tensData(1:lastIndexTens, yColTens);
end

if hasFlex == 1
    % Load flexural test data
    if bk0tb1cm2Flex == 0
        flexData = dlmread(fnameFlex, '', startRowFlex-1, 0);
    elseif bk0tb1cm2Flex == 1
        flexData = dlmread(fnameFlex, '\t', startRowFlex-1, 0);
    else
        flexData = dlmread(fnameFlex, ',', startRowFlex-1, 0);
    end
    % Read data up to maximum flexural deflection
    [xMax, lastIndexFlex] = max(flexData(:, xColFlex));
    expFlexDisp = flexData(1:lastIndexFlex, xColFlex);
    expFlexLoad = flexData(1:lastIndexFlex, yColFlex);
    expFlexStrs = (0.5*expFlexLoad*S1)*6/(b*d^2);
end

%##########################################################################
% Part F: Code to check design limits
%##########################################################################
% Calculate limits for concrete crushing and steel yielding
MMcr = envelope(10,6);
MM_ultimate = (Mcr*(-12*rho*n*kappa*mu+12*rho*n*kappa*alpha*mu+6*mu*omega+3*omega^2+3*rho*n*kappa*sqrt((rho*n*kappa-omega)^2)-3*rho^2*n^2*kappa^2-3*omega*sqrt((rho*n*kappa-omega)^2)-6*rho*n*kappa*omega+12*rho*n*kappa*alpha*omega)/(2*omega+2*mu));

if Type == 1 % Hybrid Design
    % Find beta value when strain limits are reached
    epsilon_steel_tension = (envelope(:,1).*epsilon_cr).*(alpha-envelope(:,4))./(1-alpha);
    epsilon_con_compression = (envelope(:,1).*epsilon_cr.*envelope(:,4))./(1-envelope(:,4));
    Beta_controll_1 = find(epsilon_con_compression <= 0.003);
    Beta_controll_2 = find(epsilon_steel_tension <= 0.005);
    Beta_controll = min(Beta_controll_1(length(Beta_controll_1),1), Beta_controll_2(length(Beta_controll_2),1));

    if Beta_controll_1(length(Beta_controll_1),1) < Beta_controll_2(length(Beta_controll_2),1)
       Beta_controll_str = "compression control";
       epsilon_comp = (envelope(Beta_controll_1(length(Beta_controll_1),1),1).*epsilon_cr.*envelope(Beta_controll_1(length(Beta_controll_1),1),4))./(1-envelope(Beta_controll_1(length(Beta_controll_1),1),4));
       epsilon_tension = (envelope(Beta_controll_1(length(Beta_controll_1),1),1).*epsilon_cr).*(alpha-envelope(Beta_controll_1(length(Beta_controll_1),1),4))./(1-alpha);
       
    elseif Beta_controll_1(length(Beta_controll_1),1) > Beta_controll_2(length(Beta_controll_2),1)
       Beta_controll_str = "Tension control";
       epsilon_comp = (envelope(Beta_controll_2(length(Beta_controll_2),1),1).*epsilon_cr.*envelope(Beta_controll_2(length(Beta_controll_2),1),4))./(1-envelope(Beta_controll_2(length(Beta_controll_2),1),4));
       epsilon_tension = (envelope(Beta_controll_2(length(Beta_controll_2),1),1).*epsilon_cr).*(alpha-envelope(Beta_controll_2(length(Beta_controll_2),1),4))./(1-alpha);
    
    else 
        Beta_controll_str = "strain in compression and tension is lower than ACI318";
        epsilon_comp = 0.0;
        epsilon_tension = 0.0;
    end
    
    % Get maximum moment under tension control
    MMmax = max(envelope(1:Beta_controll(1,1),6));
    
else
    % For pure FRC (no rebars)
    Beta_controll_str = "No Rebar";
    epsilon_comp = 0;
    epsilon_tension = 0;
    MMmax = max(envelope(:,6));
end

%##########################################################################
% Part G: Intersection points calculation
%##########################################################################
% Calculate key transition points in moment-curvature diagram
[Output_intersection, test] = Intersection_points_Tri_HRC(Mcr, phicr, rho, n, zeta, xi, alpha, eta, kappa, omega, tau, mu, epsilon_cr, d, E, b);

% Find corresponding indices in the envelope curve
i = 1;
for i = 1:length(Output_intersection(1:end,2))
    [~, index(i,1)] = min(abs(envelope(1:end,1) - Output_intersection(i,2)));
end

%##########################################################################
% Part H: Create additional plots
%##########################################################################
% Plot moment vs beta
fig4 = figure(4);
plot(envelope(:,1), envelope(:,6), 'k', 'linewidth', 2);
hold on
starter(1,1) = 0; 
starter(2,1) = envelope(1,1);
starter(1,6) = 0;
starter(2,6) = envelope(1,6);
plot(starter(:,1), starter(:,6), 'k', 'linewidth', 2);
xlabel('Beta');
ylabel('Moment');
xlim([-5 inf]);
ylim([0 inf]);
title(['Moment vs Beta at \mu = ', num2str(mu)]);
grid on;
saveas(fig4, 'Moment_vs_Beta.png');

% Plot normalized moment vs beta for all zones
lanbda_limit = find(envelope(:,2) <= lambda_cu);
fig5 = figure(5);
a1 = plot(output1(1:lanbda_limit(end,1),1), output1(1:lanbda_limit(end,1),5), 'k', 'linewidth', 0.005); M1 = "Zone1";
hold on
a2 = plot(output21(1:lanbda_limit(end,1),1), output21(1:lanbda_limit(end,1),5), 'm', 'linewidth', 0.005); M2 = "Zone21";
hold on
a3 = plot(output22(1:lanbda_limit(end,1),1), output22(1:lanbda_limit(end,1),5), 'm', 'linewidth', 2); M3 = "Zone22";
hold on
a4 = plot(output31(1:lanbda_limit(end,1),1), output31(1:lanbda_limit(end,1),5), 'r', 'linewidth', 0.005); M4 = "Zone31"; 
hold on
a5 = plot(output32(1:lanbda_limit(end,1),1), output32(1:lanbda_limit(end,1),5), 'r', 'linewidth', 2); M5 = "Zone32";
hold on
a6 = plot(output41(1:lanbda_limit(end,1),1), output41(1:lanbda_limit(end,1),5), 'b', 'linewidth', 0.005); M6 = "Zone41";
hold on
a7 = plot(output42(1:lanbda_limit(end,1),1), output42(1:lanbda_limit(end,1),5), 'b', 'linewidth', 2); M7 = "Zone42";
hold on
a8 = plot(output51(1:lanbda_limit(end,1),1), output51(1:lanbda_limit(end,1),5), 'g', 'linewidth', 0.005); M8 = "Zone51";
hold on
a9 = plot(output52(1:lanbda_limit(end,1),1), output52(1:lanbda_limit(end,1),5), 'g', 'linewidth', 2); M9 = "Zone52";
hold on
starter(1,1) = 0; 
starter(2,1) = envelope(1,1);
starter(1,5) = 0;
starter(2,5) = envelope(1,5);
plot(starter(:,1), starter(:,5), 'ko', 'linewidth', 1.0);
hold on
a10 = plot(envelope(:,1), envelope(:,5), 'k--', 'linewidth', 2.0); M10 = "Envelope";
legend([a1; a2; a3; a4; a5; a6; a7; a8; a9; a10], [M1; M2; M3; M4; M5; M6; M7; M8; M9; M10]);
xlabel('Beta');
ylabel('Normalized Moment');
xlim([0 beta_tu/10]);
ylim([0 max(envelope(:,5))*1.1]);
title(['Moment vs Beta at \mu = ', num2str(mu)]);
grid on;
saveas(fig5, 'Moment_vs_Beta_in_variable_zone.png');

% Plot curvature distribution along beam at different load stages
fig8 = figure(8);
plot(X, Phi(stepY, :), 'r', X, Phi(stepM, :), 'b', X, Phi(stepF, :), 'm');
grid on; 
title(['Curvature along X Axis at \mu = ', num2str(mu)]); 
xlabel('X locations (mm)'); 
ylabel('Curvature (1/mm)');
legend('Yield', 'Max', 'Failure');
saveas(fig8, 'Curvature_along_X_Axis.png');

% Plot moment distribution along beam at different load stages
fig9 = figure(9);
plot(X, Mmt(stepY, :), 'r', X, Mmt(stepM, :), 'b', X, Mmt(stepF, :), 'm');
grid on; 
title(['Moment along X Axis at \mu = ', num2str(mu)]); 
xlabel('X locations (mm)'); 
ylabel('Moment (N-mm)');
legend('Yield', 'Max', 'Failure');
saveas(fig9, 'Moment_along_X_Axis.png');

% Plot load-deflection curve
fig10 = figure(10);
plot(expFlexDisp, expFlexLoad, '-bx', delta(:,end), 2*Rs, '-ro');
grid on;
title(['Load - Deflection at \mu = ', num2str(mu)]);
xlabel('Deflection (mm)');
ylabel('Load (N)');
legend('Experiment', 'Predicted');
saveas(fig10, 'Load_Deflection.png');

% Plot equivalent stress-deflection curve
fig11 = figure(11);
plot(expFlexDisp, expFlexStrs, '-bx', delta(:,end), flexStrs, '-ro');
grid on;
title(['Equivalent Stress - Deflection at \mu = ', num2str(mu)]);
xlabel('Deflection (mm)');
ylabel('Equivalent Stress (MPa)');
legend('Experiment', 'Predicted');
saveas(fig11, 'Equivalent_Stress_Deflection.png');

% Plot moment and curvature distributions
fig12 = figure(12);
subplot(2,1,1);
plot(X(1,:), Mmt(beta_tu*10,:), 'k', 'linewidth', 2);
xlabel('Discrete Section (half of length)');
ylabel('Moment Distribution');
title(['Moment Distribution Timelapse at \mu = ', num2str(mu)]);
grid on;

subplot(2,1,2);
plot(X(1,:), Phi(beta_tu*10,:), 'm', 'linewidth', 2);
xlabel('Discrete Section (half of length)');
ylabel('Curvature Distribution');
title(['Curvature Distribution Timelapse at \mu = ', num2str(mu)]);
grid on;
saveas(fig12, 'Moment-Curvature_Distribution.png');

% Plot neutral axis depth (k) vs beta
MM_g = NaN(length(MM(:,1))+1,1);
cv_g = NaN(length(MM(:,1))+1,1);
CV_g = NaN(length(MM(:,1))+1,1);
mm_g = NaN(length(MM(:,1))+1,1);
bt_g = NaN(length(MM(:,1))+1,1); % beta 
kd_g = NaN(length(MM(:,1))+1,1); % kd
ld_g = NaN(length(MM(:,1))+1,1); % lambda
mm_g = NaN(length(MM(:,1))+1,1); % Normalized Moment
kapp_ten_g = NaN(length(MM(:,1))+1,1); % kappa tension
kapp_com_g = NaN(length(MM(:,1))+1,1); % kappa compression
nf_g = NaN(length(MM(:,1))+1,1); % netforce

MM_g(2:end,1) = envelope(:,6);
cv_g(2:end,1) = envelope(:,7);
CV_g(2:end,1) = envelope(:,8);
bt_g(2:end,1) = envelope(:,1); % beta 
kd_g(2:end,1) = envelope(:,4); % kd
ld_g(2:end,1) = envelope(:,2); % lambda
mm_g(2:end,1) = envelope(:,5); % Normalized Moment
kapp_ten_g(2:end,1) = envelope(:,3); % kappa tension
kapp_com_g(2:end,1) = envelope(:,12); % kappa compression
nf_g(2:end,1) = envelope(:,11); % netforce

MM_g(1,1) = 0;
cv_g(1,1) = 0;
CV_g(1,1) = 0;
bt_g(1,1) = 0; % beta 
kd_g(1,1) = 0; % kd
ld_g(1,1) = 0; % lambda
mm_g(1,1) = 0; % Normalized Moment
kapp_ten_g(1,1) = 0; % kappa tension
kapp_com_g(1,1) = 0; % kappa compression
nf_g(1,1) = 0; % netforce

fig13 = figure(13);
plot(bt_g, kd_g, 'r');
title(['k - beta at \mu = ', num2str(mu)]);
xlabel('Beta');
ylabel('k (neutral axis depth ratio)');
AX = axis;
axis([AX(1), AX(2), 0, 0.6]);
saveas(fig13, 'k_beta.png');

% Plot net force vs beta
fig14 = figure(14);
plot(bt_g, nf_g, 'r');
title(['Net Force - beta at \mu = ', num2str(mu)]);
xlabel('Beta');
ylabel('Net Force');
saveas(fig14, 'Net_Force_beta.png');

% Plot moment-curvature diagram with efficiency factors
fig15 = figure(15);

% Find indices of key points
for i = 1:length(Output_intersection(1:end,2))
    [~, index(i,1)] = min(abs(envelope(1:end,1) - Output_intersection(i,2)));
end

% Plot efficiency factors on right y-axis
yyaxis right;
a1 = plot(CV(1:lanbda_limit(end,1),1), envelope(1:lanbda_limit(end,1),21), '-r', 'linewidth', 2); M1 = "Concrete Compression";
hold on;
a2 = plot(CV(1:lanbda_limit(end,1),1), envelope(1:lanbda_limit(end,1),23), 'r', 'linewidth', 2); M2 = "Concrete Tension";
hold on;
a3 = plot(CV(1:lanbda_limit(end,1),1), envelope(1:lanbda_limit(end,1),22), '-b', 'linewidth', 2); M3 = "Rebar Compression";
hold on;
a4 = plot(CV(1:lanbda_limit(end,1),1), envelope(1:lanbda_limit(end,1),24), 'b', 'linewidth', 2); M4 = "Rebar Tension";
hold on;

xlabel('Curvature, 1/mm');
ylabel('Efficiency Factors');
ylim([-0.05, 1.05]);
grid on;

% Plot moment-curvature on left y-axis
yyaxis left;
b1 = plot(CV(1:lanbda_limit(end,1),1), MM(1:lanbda_limit(end,1),1), 'k', 'LineWidth', 2); N1 = "Moment-Curvature";
hold on;
b2 = scatter(Output_intersection(2,7), Output_intersection(2,6), 'o', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 100); N2 = "First Crack";
hold on;
b3 = scatter(Output_intersection(3,7), Output_intersection(3,6), 's', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 100); N3 = "Tension Rebar Yield";
hold on;
b4 = scatter(Output_intersection(4,7), Output_intersection(4,6), '^', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 100); N4 = "UHPC tensile strain limit";
hold on;
b5 = scatter(Output_intersection(5,7), Output_intersection(5,6), '*', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 100); N5 = "Concrete Compression Plastic";
hold on;

legend([b1; b2; b3; b4; b5; a1; a2; a3; a4], [N1; N2; N3; N4; N5; M1; M2; M3; M4]);
ylabel('Bending Moment (N-mm)');
title(['Force Contribution in UHPC Reinforced Concrete at \mu = ', num2str(mu), ' and \rho = ', num2str(rho)]);
saveas(fig15, 'Force_Contribution_in_UHPC_Reinforced_Concrete.png');

% Plot efficiency factors in concrete tension for different zones
fig16 = figure(16);
a1 = plot(output1(1:lanbda_limit(end,1),1), output1(1:lanbda_limit(end,1),23), 'k', 'linewidth', 0.005); M1 = "Zone1";
hold on;
a2 = plot(output21(1:lanbda_limit(end,1),1), output21(1:lanbda_limit(end,1),23), 'm', 'linewidth', 0.005); M2 = "Zone21";
hold on;
a3 = plot(output22(1:lanbda_limit(end,1),1), output22(1:lanbda_limit(end,1),23), 'm', 'linewidth', 2); M3 = "Zone22";
hold on;
a4 = plot(output31(1:lanbda_limit(end,1),1), output31(1:lanbda_limit(end,1),23), 'r', 'linewidth', 0.005); M4 = "Zone31"; 
hold on;
a5 = plot(output32(1:lanbda_limit(end,1),1), output32(1:lanbda_limit(end,1),23), 'r', 'linewidth', 2); M5 = "Zone32";
hold on;
a6 = plot(output41(1:lanbda_limit(end,1),1), output41(1:lanbda_limit(end,1),23), 'b', 'linewidth', 0.005); M6 = "Zone41";
hold on;
a7 = plot(output42(1:lanbda_limit(end,1),1), output42(1:lanbda_limit(end,1),23), 'b', 'linewidth', 2); M7 = "Zone42";
hold on;
a8 = plot(output51(1:lanbda_limit(end,1),1), output51(1:lanbda_limit(end,1),23), 'g', 'linewidth', 0.005); M8 = "Zone51";
hold on;
a9 = plot(output52(1:lanbda_limit(end,1),1), output52(1:lanbda_limit(end,1),23), 'g', 'linewidth', 2); M9 = "Zone52";
hold on;
starter(1,1) = 0; 
starter(2,1) = envelope(1,1);
starter(1,5) = 0;
starter(2,5) = envelope(1,5);
plot(starter(:,1), starter(:,5), 'ko', 'linewidth', 1.0);
hold on;
a10 = plot(envelope(:,1), envelope(:,23), 'r--', 'linewidth', 3.0); M10 = "Envelope concrete tension";
legend([a1; a2; a3; a4; a5; a6; a7; a8; a9; a10], [M1; M2; M3; M4; M5; M6; M7; M8; M9; M10]);
xlabel('\beta');
ylabel('Efficiency Factors');
ylim([-0.05, 1.05]);
xlim([envelope(1,7), envelope(end,7)]);
title(['Efficiency factor in concrete tension at \mu = ', num2str(mu), ' and \rho = ', num2str(rho)]);
grid on;
saveas(fig16, 'Efficiency_factor_in_concrete_tension.png');

% Plot efficiency factors in rebar tension for different zones
fig17 = figure(17);
a1 = plot(output1(1:lanbda_limit(end,1),1), output1(1:lanbda_limit(end,1),24), 'k', 'linewidth', 0.005); M1 = "Zone1";
hold on;
a2 = plot(output21(1:lanbda_limit(end,1),1), output21(1:lanbda_limit(end,1),24), 'm', 'linewidth', 0.005); M2 = "Zone21";
hold on;
a3 = plot(output22(1:lanbda_limit(end,1),1), output22(1:lanbda_limit(end,1),24), 'm', 'linewidth', 2); M3 = "Zone22";
hold on;
a4 = plot(output31(1:lanbda_limit(end,1),1), output31(1:lanbda_limit(end,1),24), 'r', 'linewidth', 0.005); M4 = "Zone31"; 
hold on;
a5 = plot(output32(1:lanbda_limit(end,1),1), output32(1:lanbda_limit(end,1),24), 'r', 'linewidth', 2); M5 = "Zone32";
hold on;
a6 = plot(output41(1:lanbda_limit(end,1),1), output41(1:lanbda_limit(end,1),24), 'b', 'linewidth', 0.005); M6 = "Zone41";
hold on;
a7 = plot(output42(1:lanbda_limit(end,1),1), output42(1:lanbda_limit(end,1),24), 'b', 'linewidth', 2); M7 = "Zone42";
hold on;
a8 = plot(output51(1:lanbda_limit(end,1),1), output51(1:lanbda_limit(end,1),24), 'g', 'linewidth', 0.005); M8 = "Zone51";
hold on;
a9 = plot(output52(1:lanbda_limit(end,1),1), output52(1:lanbda_limit(end,1),24), 'g', 'linewidth', 2); M9 = "Zone52";
hold on;
starter(1,1) = 0; 
starter(2,1) = envelope(1,1);
starter(1,5) = 0;
starter(2,5) = envelope(1,5);
plot(starter(:,1), starter(:,5), 'ko', 'linewidth', 1.0);
hold on;
a10 = plot(envelope(:,1), envelope(:,24), 'b--', 'linewidth', 3.0); M10 = "Envelope rebar tension";
legend([a1; a2; a3; a4; a5; a6; a7; a8; a9; a10], [M1; M2; M3; M4; M5; M6; M7; M8; M9; M10]);
xlabel('\beta');
ylabel('Efficiency Factors');
ylim([-0.05, 1.05]);
xlim([envelope(1,7), envelope(end,7)]);
title(['Efficiency factor in rebar tension at \mu = ', num2str(mu), ' and \rho = ', num2str(rho)]);
grid on;
saveas(fig17, 'Efficiency_factor_in_rebar_tension.png');

% Plot beta distribution with key points
for i = 1:1
    [~, idx1(i,1)] = min(abs(Beta(index(2,1),:) - index(i+1,1)/10));
end

for i = 1:2
    [~, idx2(i,1)] = min(abs(Beta(index(3,1),:) - index(i+1,1)/10));
end

for i = 1:3
    [~, idx3(i,1)] = min(abs(Beta(index(4,1),:) - index(i+1,1)/10));
end

fig18 = figure(18);
subplot(1,2,1); % 1 row, 2 columns, subplot 1
% Plot Beta(x) data with key points
d1 = plot(X(1,1:end), Beta(index(2,1),1:end), '-b', 'LineWidth', 2); hold on; dM1 = ['Mcrack, max Beta = ', num2str(index(2,1)/10)];
d2 = plot(X(1,1:end), Beta(index(3,1),1:end), '-k', 'LineWidth', 2); dM2 = ['Msy, max Beta = ', num2str(index(3,1)/10)];
d3 = plot(X(1,1:end), Beta(index(4,1),1:end), '-g', 'LineWidth', 2); dM3 = ['Mloc, max Beta = ', num2str(index(4,1)/10)];

% Scatter plots for key points
d4 = scatter(X(1,idx1(1,1)), index(2,1)/10, 'o', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 100); dM4 = "Mcrack";
hold on;
d5 = scatter(X(1,idx2(1,1)), index(2,1)/10, 'o', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 100); dM5 = "Mcrack";
hold on;
d6 = scatter(X(1,idx2(2,1)), index(3,1)/10, 's', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 100); dM6 = "Msy";
hold on;
d7 = scatter(X(1,idx3(1,1)), index(2,1)/10, 'o', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 100); dM7 = "Mcrack";
hold on;
d8 = scatter(X(1,idx3(2,1)), index(3,1)/10, 's', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 100); dM8 = "Msy";
hold on;
d9 = scatter(X(1,idx3(3,1)), index(4,1)/10, '^', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 100); dM9 = "Mloc";
hold on;
xlabel('X (Half Beam)');
ylabel('Beta(x)');
title('Beta(x) vs. X');
legend([d1; d2; d3; d4; d5; d6; d7; d8; d9], [dM1; dM2; dM3; dM4; dM5; dM6; dM7; dM8; dM9], 'Location', 'best');
grid on;

% Add text box with parameters
textStr1 = sprintf(['Parameters' ...
    '\ncase = %.0f' ...
    '\nfc = %.0f' ...
    '\nE = %.0f' ...
    '\n\\epsilon_c_r = %.6f' ...
    '\n\\mu = %.2f' ...
    '\n\\beta_1 = %.2f' ...
    '\n\\eta = %.3f' ...
    '\n\\gamma = %.1f' ...
    '\n\\omega = %.2f' ...
    '\n\\kappa = %.2f' ...
    '\nn = %.2f' ...
    '\n\\rho = %.4f' ...
    '\n\\alpha = %.2f' ...
    '\n\\zeta = %.2f'], test, E*epsilon_cr*omega, E, epsilon_cr, mu, tau, eta, xi, omega, kappa, n, rho, alpha, zeta);

% Position the text box
xLimits = xlim;
yLimits = ylim;
xPos1 = xLimits(1) + 0.1* (xLimits(2) - xLimits(1));
yPos1 = yLimits(1) + 0.6 * (yLimits(2) - yLimits(1));
% Add text box to subplot 1
text(xPos1, yPos1, textStr1, 'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');

subplot(1,2,2); % 1 row, 2 columns, subplot 2
% Left Y-Axis: Moment-Curvature with key points
yyaxis left;
b1 = plot(envelope(1:lanbda_limit(end,1),7), envelope(1:lanbda_limit(end,1),5), 'k', 'LineWidth', 2); hold on; N1 = "Moment-Curvature";
b2 = scatter(Output_intersection(2,5), Output_intersection(2,4), 100, 'k', 'filled', 'o'); N2 = "First Crack";
b3 = scatter(Output_intersection(3,5), Output_intersection(3,4), 100, 'k', 'filled', 's'); N3 = "Tension Rebar Yield";
b4 = scatter(Output_intersection(4,5), Output_intersection(4,4), 100, 'k', 'filled', '^'); N4 = "UHPC strain localization";
b5 = scatter(Output_intersection(5,5), Output_intersection(5,4), 100, 'k', 'filled', '*'); N5 = "Concrete Compression Plastic";
ylabel('Normalized Bending Moment');

% Right Y-Axis: Efficiency Factors
yyaxis right;
a1 = plot(envelope(1:lanbda_limit(end,1),7), envelope(1:lanbda_limit(end,1),21), '-r', 'LineWidth', 2); M1 = "Concrete Compression";
hold on;
a2 = plot(envelope(1:lanbda_limit(end,1),7), envelope(1:lanbda_limit(end,1),23), 'r', 'LineWidth', 2); M2 = "Concrete Tension";
a3 = plot(envelope(1:lanbda_limit(end,1),7), envelope(1:lanbda_limit(end,1),22), '-b', 'LineWidth', 2); M3 = "Rebar Compression";
a4 = plot(envelope(1:lanbda_limit(end,1),7), envelope(1:lanbda_limit(end,1),24), 'b', 'LineWidth', 2); M4 = "Rebar Tension";
ylabel('Efficiency Factors');
ylim([-0.05, 1.05]);

xlabel('Normalized Curvature');
title(['Force Contribution in UHPC Reinforced Concrete at \mu = ', num2str(mu), ' and \rho = ', num2str(rho)]);
grid on;

% Add legend
legend([b1; b2; b3; b4; b5; a1; a2; a3; a4], ...
       [N1; N2; N3; N4; N5; M1; M2; M3; M4], ...
       'Location', 'best');
saveas(fig18, 'Beta-X_Combined_Analysis.png');

%##########################################################################
% Part I: Output results to files
%##########################################################################
% Open file for writing main output
fid = fopen(fname , 'w');

N1 = length(expTensStrn);
N2 = length(strain);
N3 = length(envelope(1:1:lanbda_limit(end,1),11));
N4 = length(expFlexDisp);
N5 = length(delta(:,end));
N6 = length(output1(1:1:lanbda_limit(end,1),1));
N7 = length(X);
NJ = [N1,N1,N2,N2,N3,N3,N3,N3,N3,N4,N4,N4,N5,N5,N5,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,...
     N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N7,N7,N7,N7,N7,N7,N7];

for i=1:N1
    %'expTensStrn','expTensStrs'
    numTable(i,1:2)   = [expTensStrn(i),expTensStrs(i)];
end

for i=1:N2
    %'strain','stress',
    numTable(i,3:4)   = [strain(i),stress(i)];
end

for i=1:N3
    %'Bottom strain','Bottom st strain','Top strain','Top st strain','Equivalent Stress'
    numTable(i,5:9)   = [envelope(i,13),envelope(i,15),envelope(i,14),envelope(i,16),EqualFlexStress(i)];
end

for i=1:N4
    %'expFlexDisp','expFlexLoad','expFlexStrs'
    numTable(i,10:12) = [expFlexDisp(i),expFlexLoad(i),expFlexStrs(i)];
end

for i=1:N5
    %'delta','totalload','flexStrs',
    numTable(i,13:15) = [delta(i,end),2*Rs(i),flexStrs(i)];
end

for i=1:N6
    %'bt','ld','kappa_ten','kappa_comp','kd','Net Force','CV envelope','MM envelope','cv envelope','mm envelope',
    numTable(i,16:52) = [envelope(i,1),envelope(i,2),envelope(i,3),envelope(i,12),envelope(i,4),envelope(i,11),envelope(i,8),envelope(i,6),envelope(i,7),envelope(i,5),output1(i,7),output1(i,5),output21(i,7),output21(i,5),output22(i,7),...
                         output22(i,5),output31(i,7),output31(i,5),output32(i,7),output32(i,5),output41(i,7),output41(i,5),output42(i,7),output42(i,5),output51(i,7),output51(i,5),output52(i,7),output52(i,5),Matrix_sum_com(i,1)*E*b*h,Matrix_sum_ten(i,1)*E*b*h,Rebar_com(i,1)*E*b*h,Rebar_ten(i,1)*E*b*h,Total_force(i,1)*E*b*h,envelope(i,21),envelope(i,23)...
                         ,envelope(i,22),envelope(i,24)];
end

for i=1:N7
    %'Location_X','M_stepY','M_stepM','M_stepF','CV_stepY'
    numTable(i,53:59) = [X(i),Mmt(stepY,i),Mmt(stepM,i),Mmt(stepF,i),Phi(stepY, i),Phi(stepM, i),Phi(stepF, i)];
end

fid = fopen(fname,'w');
fprintf (fid,'%1s\n\n','Load deflection response of three or four point bending test predicted by uniaxial stress strain model');

if hasTens == 1
    fprintf (fid,'%12s %1s\n','fnameTens =',fnameTens);
end
if hasFlex == 1
    fprintf (fid,'%12s %1s\n','fnameFlex =',fnameFlex);
end
fprintf (fid,'\n');
fprintf (fid,' If 3 point bending, Lp and unLoadFactor are used, c is ignored\n');
fprintf (fid,' If 4              , Lp is ignored,                unLoadFactor and c are used\n');
fprintf (fid,' *****************************************************************************************************************************************************************************************************************************************************\n\n');

fprintf (fid,'%12s,%12s,%12s,%12s,%12s\n','pointBend','Lp','c','unLoadFactor1','unLoadFactor2');
fprintf (fid,'%12.4g,%12.4g,%12.4g,%12.4g,%12.4g\n', pointBend,Lp,c,unLoadFactor1,unLoadFactor2);
fprintf (fid,' *****************************************************************************************************************************************************************************************************************************************************\n\n');

fprintf (fid,'%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s\n',...
             'b','h','L','ecr','E','mu','beta_tu','gamma','omega','lambda_cu','n','kappa','rho','zeta','tau','alpha','eta','MMcr','MMmax','tor');
fprintf (fid,'%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g\n',...
             b, h, L, epsilon_cr,E,mu,beta_tu,xi,omega,lambda_cu,n,kappa,rho,zeta,tau,alpha,eta,MMcr,MMmax,tor);
fprintf (fid,' *****************************************************************************************************************************************************************************************************************************************************\n\n');

fprintf (fid,'%12s,%12s,%12s,%12s,%12s\n','subMC(1)','subMC(2)','subMC(3)','nSeg(1)','nSeg(2)');
fprintf (fid,'%12.4g,%12.4g,%12.4g,%12.4g,%12.4g\n',subMC(1), subMC(2), subMC(3), nSeg(1), nSeg(2));
fprintf (fid,' *****************************************************************************************************************************************************************************************************************************************************\n\n');

% print table
fprintf (fid,'%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s\n',...
             'expTensStrn','expTensStrs','strain','stress','Bottom strain','Bottom st strain','Top strain','Top st strain','Equivalent Stress','expFlexDisp','expFlexLoad','expFlexStrs','delta','totalload','flexStrs','bt','ld','kappa_ten','kappa_comp','kd','Net Force',...
             'CV envelope','MM envelope','cv envelope','mm envelope','cv zone1','mm zone1','cv zone21','mm zone21','cv zone22','mm zone22','cv zone31',...
             'mm zone31','cv zone32','mm zone32','cv zone41','mm zone41','cv zone42','mm zone42',...
             'cv zone51','mm zone51','cv zone52','mm zone52','Matrix com force','Matrix ten force','Rebar com force','Rebar ten force','Total force',...
             'Matrix com froce ratio','Matrix ten force ratio','Rebar com force ratio','Rebar ten force ratio','Location_X','M_stepY','M_stepM','M_stepF','CV_stepY',...
             'CV_stepM','CV_stepF');

[nRows, nCols] = size(numTable); 
for i=1:nRows
    for j=1:nCols
        if (i > NJ(j)) & (numTable(i,j) == 0)
            fprintf (fid,'%12s,', '');
        else
            fprintf (fid,'%12.4g,', numTable(i,j));
        end
    end
    fprintf (fid,'\n');
end
fclose(fid);

%##########################################################################
% (b) Store Efficiency force at the intersection points 
%##########################################################################

% Open file for writing
fileID = fopen(fname2 , 'w');

% Write data to the file
fprintf(fileID, '#Input parameters \n');
fprintf(fileID, '#Compressive strength, fc = %.3f\n',fc);
fprintf(fileID, '#Concrete young modulus, E = %f\n',E);
fprintf(fileID, '#Normalized concrete compressive modulus, gamma = %.3f\n',xi);
fprintf(fileID, '#Transition tensile strain, beta_1 = %.3f\n',tau);
fprintf(fileID, '#Normalized depth of steel reinforcement (d/h), alpha = %.3f\n',alpha);
fprintf(fileID, '#Compression rebar ratio (As_prime/As),zeta = %.3f\n',zeta);
fprintf(fileID, '#Normalized concrete compressive yield strain, omega = %.3f\n',omega);
fprintf(fileID, '#Cracking strain, epsilon_cr = %f\n',epsilon_cr);
fprintf(fileID, '#Rebar Grade, 29,000,000 psi');
fprintf(fileID, '#Rebar modulus ratio, n = %.3f\n',n);
fprintf(fileID, '#Ultimate normalized compressive strain, lambda_cu= %.3f\n\n\n',lambda_cu);
N1 = length(index(:,1));
NJ = [N1,N1,N1,N1] ;

for i=1:N1
    NumTable(i,1:4)   = [envelope(index(i,1),21),envelope(index(i,1),22),envelope(index(i,1),23),envelope(index(i,1),24)];
end

% print table
fprintf (fileID,'%12s,%12s,%12s,%12s\n',...
             'ConcreteCompression','RebarCompression','ConcreteTension','Rebartension');

[nRows, nCols] = size(NumTable); 
for i=1:nRows
    for j=1:nCols
        if (i > NJ(j)) & (NumTable(i,j) == 0)
            fprintf (fileID,'%12s,', '');
        else
            fprintf (fileID,'%12.4g,', NumTable(i,j));
        end
    end
    fprintf (fileID,'\n');
end
fclose(fileID);
% Close the file
disp('Results have been written to output.txt');

%##########################################################################
% End of Main Program
%##########################################################################