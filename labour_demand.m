%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               MACROECONOMICS III - DYNAMIC PROGRAMMING
%%                 Labour demand with adjustment costs

%                          José Pedro Garcia

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script is meant to develop a code which solves the labour demand
% model shown in section 4.2 of the class notes. 



%--------------------------------------------------------------------------
% HOUSEKEEPING
%--------------------------------------------------------------------------

clear;
clc;

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Set directory here:
directory = '/Users/josegarcia/Desktop/Zé Pedro/EUI 2023-2024/Block 4/Teaching/Class notes';
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cd(directory);

% Set parameters of the model:
aalpha                                     = 0.73;                         % Curvature of the revenue function
bbeta                                      = 0.9;                          % Time discount factor

rrho                                       = 0.9;                          % Serial correlation of profitability
ssigma                                     = 0.6;                          % Std deviation of profitability innovations

w0                                         = 0.33;                         % Base wage
w1                                         = 0.01;                         % Overtime pay rate
eeta                                       = 1.9;                          % Elasticity of hourly wage wrt hours worked

ttheta                                     = ...
    (eeta * w1 / aalpha) ^ (1 / (aalpha - eeta));                          % Composite parameter (see notes)

GGamma_f                                   = 2.9;                          % Fixed cost of firing
ggamma_f                                   = 1.3;                          % Linear cost of firing
GGamma_h                                   = 1;                            % Fixed cost of hiring
ggamma_h                                   = 1.1;                          % Linear cost of hiring

% Set numerical parameters and grid for profitability:
nProf                                      = 20;                           % Number of points in profitability grid
mmu                                        = exp(0);                       % Mean of TFPR process
kkappa                                     = 2;                            % How many standard deviations should the grid be in length

[vProfGrid, mProfTransition, vProfErgodic] = ...
    compute_tauchen(rrho, ssigma, mmu, nProf, kkappa);                     % Tauchen method spits back the discrete grid for profitability

% Set numerical parameters and grid for labour:
nLabour                                    = 750;                          % Number of points in labour grid

minLabour                                  = 0;    
maxLabour                                  = log(500) / log(10);
vLabourGrid                                = ...
    logspace(minLabour, maxLabour, nLabour)';                              % Grid for labour
 
% Define grids for state and choice sets:
% CONVENTION:
% ROWS    = State variables (productivity + current capital)
% COLUMNS = Choice variables (future capital)
nState                                     = nProf * nLabour;              % Number of realisations of the state-space
[mProfGrid, mLabourGrid]                   = ndgrid(vProfGrid, vLabourGrid); 
mStateGrid                                 = [mProfGrid(:) mLabourGrid(:)];

nChoice                                    = nLabour;



%--------------------------------------------------------------------------
% SET UP PROFIT FUNCTION
%--------------------------------------------------------------------------

% Objects of interest:
mProfState    = repmat(mStateGrid(:,1), 1, nChoice);
mLabourState  = repmat(mStateGrid(:,2), 1, nChoice);
mLabourChoice = repmat(vLabourGrid', nState, 1);

% From the FOC wrt hours (h) we know that the optimal choice of h is a
% function of two arguments: profitability and new labour demand.
mHours        = ttheta * (exp(mProfState) .^ (1 / (eeta - aalpha))) .* ...
    (mLabourChoice .^ ((1 - aalpha) / (aalpha - eeta)));

% Define the revenue function:
mRevenue      = exp(mProfState) .* ...
    ((mLabourChoice .* mHours) .^ aalpha);

% Define labour adjustment costs:
mDeltaLabour  = mLabourChoice - mLabourState;

mHiringCosts  = (GGamma_h + ggamma_f * abs(mDeltaLabour)) .* (mDeltaLabour > 0);
mFiringCosts  = (GGamma_f + ggamma_h * abs(mDeltaLabour)) .* (mDeltaLabour < 0);
mAdjCosts     = mHiringCosts + mFiringCosts;

% Define the wage function:
mWage         = mLabourChoice .* (w0 + w1 * (mHours .^ eeta));

% Define the profit function:
mProfit       = mRevenue - mAdjCosts - mWage;



%--------------------------------------------------------------------------
% VALUE FUNCTION ITERATION
%--------------------------------------------------------------------------

% Pre-allocate necessary objects:
vValue        = zeros(nState, 1);                                          % Initial guess for the value function
dampening     = 0;                                                         % Weight placed on previous iteration
toleranceVFI  = 1e-8;                                                      % Tolerance allowed for convergence
iterationVFI  = 1;                                                         % To keep track of number of iterations
maxIterations = 1e4;                                                       % Maximum number of iterations allowed
difference    = 100;                                                       % To kick off iteration
vDifference   = zeros(maxIterations, 1);                                   % Record the convergence for easier debugging

tic;
while abs(difference) > toleranceVFI && iterationVFI <= maxIterations

    % Compute expected value:
    mEValue                    = repmat(mProfTransition * reshape(vValue, nProf, nLabour), ...
        nLabour, 1);

    % Compute objective function:
    mObjectiveFunction         = mProfit + bbeta * mEValue;

    % Apply the max operator. Each row is a realisation of the state
    % vector, so we just take the max for each row and record the column
    % number where the max is located:
    [vValueNew, vLabourIndex] = max(mObjectiveFunction, [], 2); 

    % Update:
    difference                 = max(abs(vValueNew - vValue));
    vDifference(iterationVFI)  = difference;
    vValue                     = dampening * vValue + ...
        (1 - dampening) * vValueNew;
    iterationVFI               = iterationVFI + 1;    

end
time = toc;

X = ['Solution found in ', num2str(iterationVFI - 1), ' iterations. ' ...
    'It took ', num2str(time), ' seconds.'];
disp(X);

% Retrieve value and capital policy functions:
% CONVENTION:
% ROWS    = State profitability
% COLUMNS = State labour
mValue         = reshape(vValue, nProf, nLabour);
mPolicyLabour  = reshape(vLabourGrid(vLabourIndex), nProf, nLabour);
mPolicyHours   = ttheta .* (exp(repmat(vProfGrid, 1, nLabour)) .^ (1 / (eeta - aalpha))) .* ...
    (mPolicyLabour .^ ((1 - aalpha) / (aalpha - eeta)));



%--------------------------------------------------------------------------
% REPORT
%--------------------------------------------------------------------------

% Store all output in the corresponding folder:
results = strcat(directory, '/LabourDemand');
cd(results);

% Labour grid:
f1 = figure;
scatter(1 : nLabour, vLabourGrid);
hold on
plot(1 : nLabour, 1 : nLabour, 'k--', 'LineWidth', 0.5);
hold off
xlabel('Ranking of grid points', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Grid levels', 'interpreter', 'latex', 'FontSize', 14);
legend('Labour grid', '45o line','location','Southeast', 'interpreter', 'latex')
title('\textbf{The choice of grid for labour}', 'interpreter', 'latex','fontsize',16);
savefig("labour_grid")
saveas(gcf,'labour_grid.jpg')


% Policy function for employment
f2 = figure;
plot(vLabourGrid, mPolicyLabour(9, :), 'LineWidth', 1.8);
hold on
plot(vLabourGrid, mPolicyLabour(16, :), 'LineWidth', 1.8);
hold on
plot(1 : 10 ^ maxLabour, 1 : 10 ^ maxLabour, 'k:', 'LineWidth', 0.5);
hold off
xlabel('Inherited labour: $n_{-1}$', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Labour demand choice: $n$', 'interpreter', 'latex', 'FontSize', 14);
legend('Labour demand: $n=\psi(A_{9},\cdot)$','Labour demand: $n=\psi(A_{16},\cdot)$','45o line','location','Southeast', 'interpreter', 'latex')
title('\textbf{Policy Function for employment}', 'interpreter', 'latex','fontsize',16);
savefig("pol_n")
saveas(gcf,'pol_n.jpg')


% Policy function for average hours per worker
f3 = figure;
plot(vLabourGrid, mPolicyHours(9, :), 'LineWidth', 1.8);
hold on
plot(vLabourGrid, mPolicyHours(16, :), 'LineWidth', 1.8);
hold off
xlabel('Inherited labour: $n_{-1}$', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Average hours of work: $h$', 'interpreter', 'latex', 'FontSize', 14);
legend('Avg hours: $h=\phi(A_{9},\cdot)$','Avg hours: $h=\phi(A_{16},\cdot)$', 'location','east', 'interpreter', 'latex')
title('\textbf{Policy Function for hours per worker}', 'interpreter', 'latex','fontsize',16);
savefig("pol_h")
saveas(gcf,'pol_h.jpg')