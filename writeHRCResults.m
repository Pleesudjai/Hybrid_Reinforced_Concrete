function writeHRCResults(fname, fname2, envelope, output1, output21, output22, ...
                        output31, output32, output41, output42, output51, output52, ...
                        lanbda_limit, Matrix_sum_com, Matrix_sum_ten, Rebar_com, Rebar_ten, ...
                        Total_force, delta, Rs, flexStrs, X, Mmt, Phi, EqualFlexStress, ...
                        strain, stress, expTensStrn, expTensStrs, expFlexDisp, expFlexLoad, expFlexStrs, ...
                        b, h, L, epsilon_cr, E, mu, beta_tu, xi, omega, lambda_cu, n, kappa, ...
                        rho, zeta, tau, alpha, eta, MMcr, MMmax, tor, subMC, nSeg, ...
                        pointBend, Lp, c, unLoadFactor1, unLoadFactor2, hasTens, fnameTens, ...
                        hasFlex, fnameFlex, index, stepY, stepM, stepF,fc)
% writeHRCResults - Writes analysis results to output files
%
% This function handles writing all simulation data to output files.
% It creates two files: the main output data file and an efficiency factors file.
%
% Parameters:
%   fname: Name of the main output file
%   fname2: Name of the efficiency factors output file
%   envelope: Envelope curve data
%   output1, output21, etc.: Data from different analysis zones
%   lanbda_limit: Lambda limit index
%   Matrix_sum_com: Matrix compression force data
%   Matrix_sum_ten: Matrix tension force data
%   Rebar_com: Rebar compression force data
%   Rebar_ten: Rebar tension force data
%   Total_force: Total force data
%   delta: Deflection data
%   Rs: Support reaction force data
%   flexStrs: Flexural stress data
%   X: Location data along beam
%   Mmt: Moment distribution data
%   Phi: Curvature distribution data
%   EqualFlexStress: Equivalent flexural stress data
%   strain, stress: Material model data
%   expTensStrn, expTensStrs: Experimental tension data
%   expFlexDisp, expFlexLoad, expFlexStrs: Experimental flexural data
%   b, h, L: Beam geometry parameters
%   epsilon_cr, E, mu, beta_tu, etc.: Material model parameters
%   pointBend, Lp, c, etc.: Testing configuration parameters
%   hasTens, fnameTens: Tension test data parameters
%   hasFlex, fnameFlex: Flexural test data parameters
%   index: Indices of intersection points
%   stepY: Step index for yield point
%   stepM: Step index for maximum moment
%   stepF: Step index for failure
%   fc: Compressive strength of concrete

%##########################################################################
% (a) Create main output file with all analysis data
%##########################################################################
% Open file for writing main output
fid = fopen(fname, 'w');

% Get sizes for data arrays
N1 = length(expTensStrn);
N2 = length(strain);
N3 = length(envelope(1:1:lanbda_limit(end,1),11));
N4 = length(expFlexDisp);
N5 = length(delta(:,end));
N6 = length(output1(1:1:lanbda_limit(end,1),1));
N7 = length(X);

% Create array of row counts for each column
NJ = [N1,N1,N2,N2,N3,N3,N3,N3,N3,N4,N4,N4,N5,N5,N5,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,...
     N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N6,N7,N7,N7,N7,N7,N7,N7];

% Initialize the table with data
% Experimental tension data
for i = 1:N1
    numTable(i,1:2) = [expTensStrn(i), expTensStrs(i)];
end

% Material model data
for i = 1:N2
    numTable(i,3:4) = [strain(i), stress(i)];
end

% Strain and stress data from envelope
for i = 1:N3
    numTable(i,5:9) = [envelope(i,13), envelope(i,15), envelope(i,14), envelope(i,16), EqualFlexStress(i)];
end

% Experimental flexural data
for i = 1:N4
    numTable(i,10:12) = [expFlexDisp(i), expFlexLoad(i), expFlexStrs(i)];
end

% Load-deflection data
for i = 1:N5
    numTable(i,13:15) = [delta(i,end), 2*Rs(i), flexStrs(i)];
end

% Detailed analysis data from all zones
for i = 1:N6
    numTable(i,16:52) = [envelope(i,1), envelope(i,2), envelope(i,3), envelope(i,12), envelope(i,4), ...
                         envelope(i,11), envelope(i,8), envelope(i,6), envelope(i,7), envelope(i,5), ...
                         output1(i,7), output1(i,5), output21(i,7), output21(i,5), output22(i,7), ...
                         output22(i,5), output31(i,7), output31(i,5), output32(i,7), output32(i,5), ...
                         output41(i,7), output41(i,5), output42(i,7), output42(i,5), output51(i,7), ...
                         output51(i,5), output52(i,7), output52(i,5), Matrix_sum_com(i,1)*E*b*h, ...
                         Matrix_sum_ten(i,1)*E*b*h, Rebar_com(i,1)*E*b*h, Rebar_ten(i,1)*E*b*h, ...
                         Total_force(i,1)*E*b*h, envelope(i,21), envelope(i,23), envelope(i,22), envelope(i,24)];
end

% Location and distribution data
for i = 1:N7
    numTable(i,53:59) = [X(i), Mmt(stepY,i), Mmt(stepM,i), Mmt(stepF,i), Phi(stepY,i), Phi(stepM,i), Phi(stepF,i)];
end

% Write header and setup information
fprintf(fid, '%1s\n\n', 'Load deflection response of three or four point bending test predicted by uniaxial stress strain model');

if hasTens == 1
    fprintf(fid, '%12s %1s\n', 'fnameTens =', fnameTens);
end
if hasFlex == 1
    fprintf(fid, '%12s %1s\n', 'fnameFlex =', fnameFlex);
end
fprintf(fid, '\n');
fprintf(fid, ' If 3 point bending, Lp and unLoadFactor are used, c is ignored\n');
fprintf(fid, ' If 4              , Lp is ignored,                unLoadFactor and c are used\n');
fprintf(fid, ' *****************************************************************************************************************************************************************************************************************************************************\n\n');

% Write test configuration parameters
fprintf(fid, '%12s,%12s,%12s,%12s,%12s\n', 'pointBend', 'Lp', 'c', 'unLoadFactor1', 'unLoadFactor2');
fprintf(fid, '%12.4g,%12.4g,%12.4g,%12.4g,%12.4g\n', pointBend, Lp, c, unLoadFactor1, unLoadFactor2);
fprintf(fid, ' *****************************************************************************************************************************************************************************************************************************************************\n\n');

% Write material and geometry parameters
fprintf(fid, '%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s\n',...
         'b', 'h', 'L', 'ecr', 'E', 'mu', 'beta_tu', 'gamma', 'omega', 'lambda_cu', 'n', 'kappa', 'rho', 'zeta', 'tau', 'alpha', 'eta', 'MMcr', 'MMmax', 'tor');
fprintf(fid, '%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g,%12.4g\n',...
         b, h, L, epsilon_cr, E, mu, beta_tu, xi, omega, lambda_cu, n, kappa, rho, zeta, tau, alpha, eta, MMcr, MMmax, tor);
fprintf(fid, ' *****************************************************************************************************************************************************************************************************************************************************\n\n');

% Write discretization parameters
fprintf(fid, '%12s,%12s,%12s,%12s,%12s\n', 'subMC(1)', 'subMC(2)', 'subMC(3)', 'nSeg(1)', 'nSeg(2)');
fprintf(fid, '%12.4g,%12.4g,%12.4g,%12.4g,%12.4g\n', subMC(1), subMC(2), subMC(3), nSeg(1), nSeg(2));
fprintf(fid, ' *****************************************************************************************************************************************************************************************************************************************************\n\n');

% Write column headers for data table
fprintf(fid, '%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s\n',...
             'expTensStrn','expTensStrs','strain','stress','Bottom strain','Bottom st strain','Top strain','Top st strain','Equivalent Stress','expFlexDisp','expFlexLoad','expFlexStrs','delta','totalload','flexStrs','bt','ld','kappa_ten','kappa_comp','kd','Net Force',...
             'CV envelope','MM envelope','cv envelope','mm envelope','cv zone1','mm zone1','cv zone21','mm zone21','cv zone22','mm zone22','cv zone31',...
             'mm zone31','cv zone32','mm zone32','cv zone41','mm zone41','cv zone42','mm zone42',...
             'cv zone51','mm zone51','cv zone52','mm zone52','Matrix com force','Matrix ten force','Rebar com force','Rebar ten force','Total force',...
             'Matrix com froce ratio','Matrix ten force ratio','Rebar com force ratio','Rebar ten force ratio','Location_X','M_stepY','M_stepM','M_stepF','CV_stepY',...
             'CV_stepM','CV_stepF');

% Write data table
[nRows, nCols] = size(numTable); 
for i = 1:nRows
    for j = 1:nCols
        if (i > NJ(j)) && (numTable(i,j) == 0)
            fprintf(fid, '%12s,', '');
        else
            fprintf(fid, '%12.4g,', numTable(i,j));
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);

%##########################################################################
% (b) Create efficiency factors output file
%##########################################################################
% Open file for writing
fileID = fopen(fname2, 'w');

% Write header information with material parameters
fprintf(fileID, '#Input parameters \n');
fprintf(fileID, '#Compressive strength, fc = %.3f\n', fc);
fprintf(fileID, '#Concrete young modulus, E = %f\n', E);
fprintf(fileID, '#Normalized concrete compressive modulus, gamma = %.3f\n', xi);
fprintf(fileID, '#Transition tensile strain, beta_1 = %.3f\n', tau);
fprintf(fileID, '#Normalized depth of steel reinforcement (d/h), alpha = %.3f\n', alpha);
fprintf(fileID, '#Compression rebar ratio (As_prime/As), zeta = %.3f\n', zeta);
fprintf(fileID, '#Normalized concrete compressive yield strain, omega = %.3f\n', omega);
fprintf(fileID, '#Cracking strain, epsilon_cr = %f\n', epsilon_cr);
fprintf(fileID, '#Rebar Grade, 29,000,000 psi\n');
fprintf(fileID, '#Rebar modulus ratio, n = %.3f\n', n);
fprintf(fileID, '#Ultimate normalized compressive strain, lambda_cu = %.3f\n\n\n', lambda_cu);

% Prepare efficiency data
N1 = length(index(:,1));
NJ = [N1, N1, N1, N1];

for i = 1:N1
    NumTable(i,1:4) = [envelope(index(i,1),21), envelope(index(i,1),22), envelope(index(i,1),23), envelope(index(i,1),24)];
end

% Write column headers
fprintf(fileID, '%12s,%12s,%12s,%12s\n', 'ConcreteCompression', 'RebarCompression', 'ConcreteTension', 'RebarTension');

% Write efficiency data
[nRows, nCols] = size(NumTable); 
for i = 1:nRows
    for j = 1:nCols
        if (i > NJ(j)) && (NumTable(i,j) == 0)
            fprintf(fileID, '%12s,', '');
        else
            fprintf(fileID, '%12.4g,', NumTable(i,j));
        end
    end
    fprintf(fileID, '\n');
end
fclose(fileID);

% Display success message
disp('Results have been written to output files.');

end