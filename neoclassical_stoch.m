%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               MACROECONOMICS III - DYNAMIC PROGRAMMING

%                          José Pedro Garcia

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script is meant to develop a code which solves the neoclassical
% growth model *with* uncertainty.


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
bbeta     = 0.9;                                                           % Time discount factor
aalpha    = 0.36;                                                          % Capital share
ddelta    = 0.025;                                                         % Depreciation rate

rrho      = 0.95;                                                          % Serial correlation of idiosyncratic productivity
ssigma    = 0.007;                                                         % Std deviation of innovations to idiosyncratic productivity

% The steady-state level of capital is given as follows (see notes):
ssCapital = ((aalpha * bbeta) / (1 - bbeta * (1 - ddelta))) ^ ...
    (1 / (1 - aalpha));

% Set numerical parameters and grid for productivity:
nProd     = 20;                                                            % Number of grid points for productivity
mmu       = exp(0);                                                        % Mean of TFP process
kkappa    = 2;                                                             % How many standard deviations should the grid be in length

[vProdGrid, mProdTransition, vProdErgodic] = ...
    compute_tauchen(rrho, ssigma, mmu, nProd, kkappa);                     % Tauchen method spits back the discrete grid for productivity

% Set numerical parameters and grid for capital:
nCapital                       = 500;
vCapitalGrid                   = linspace(0.5 * ssCapital, 1.5 * ssCapital, nCapital)';

% Define grids for state and choice sets:
% CONVENTION:
% ROWS    = State variables (productivity + current capital)
% COLUMNS = Choice variables (future capital)
nState                         = nProd * nCapital;                         % Number of realisations of the state-space
[mProdGrid, mCapitalGrid]      = ndgrid(vProdGrid, vCapitalGrid); 
mStateGrid                     = [mProdGrid(:) mCapitalGrid(:)];

nChoice                        = nCapital;

% Compute all possible consumption levels:
mProdState                     = repmat(mStateGrid(:,1), 1, nChoice);
mCapitalState                  = repmat(mStateGrid(:,2), 1, nChoice);
mCapitalChoice                 = repmat(vCapitalGrid', nState, 1);

mConsumption                   = (exp(mProdState) .* (mCapitalState .^ aalpha)) + mCapitalState * (1 - ddelta) - ...
    mCapitalChoice;
mConsumption(mConsumption < 0) = 0;                                        % Impose non-negativity constraint on consumption



%--------------------------------------------------------------------------
% VALUE FUNCTION ITERATION
%--------------------------------------------------------------------------

% Pre-allocate necessary objects:
vValue        = zeros(nState, 1);                                          % Initial guess for the value function
dampening     = 0;                                                         % Weight placed on previous iteration
toleranceVFI  = 1e-8;                                                     % Tolerance allowed for convergence
iterationVFI  = 1;                                                         % To keep track of number of iterations
maxIterations = 1e4;                                                       % Maximum number of iterations allowed
difference    = 100;                                                       % To kick off iteration
vDifference   = zeros(maxIterations, 1);                                   % Record the convergence for easier debugging

tic;
while abs(difference) > toleranceVFI && iterationVFI <= maxIterations

    % Compute expected value:
    mEValue                    = repmat(mProdTransition * reshape(vValue, nProd, nCapital), ...
        nCapital, 1);

    % Compute objective function:
    mObjectiveFunction         = log(mConsumption) + bbeta * mEValue;

    % Apply the max operator. Each row is a realisation of the state
    % vector, so we just take the max for each row and record the column
    % number where the max is located:
    [vValueNew, vCapitalIndex] = max(mObjectiveFunction, [], 2); 

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
% ROWS    = State productivity
% COLUMNS = State capital
mValue         = reshape(vValue, nProd, nCapital);

mPolicyCapital = reshape(vCapitalGrid(vCapitalIndex), nProd, nCapital);
mIndexCapital  = reshape(vCapitalIndex, nProd, nCapital);



%--------------------------------------------------------------------------
% SIMULATION OF CAPITAL ACCUMULATION THROUGHOUT TIME
%--------------------------------------------------------------------------

% To ensure replicability of results:
rng(1893);                             

% Pre-allocate vectors to fill in with simulation:
nTime = 1e4;

vCapitalPathIndex      = zeros(nTime, 1);
vCapitalPathIndex(1,1) = round(nCapital / 2);                              % Economy starts at the mean point of capital in the grid

vProdPathIndex         = zeros(nTime, 1);
vProdPathIndex(1,1)    = randsample(1:nProd, 1, true, vProdErgodic);       % Draw an initial value for productivity

for t = 2 : nTime

    % Iterate on productivity indices:
    vProdPathIndex(t)    = randsample(1:nProd, 1, true, ...
        mProdTransition(vProdPathIndex(t - 1), :));

    % Iterate on capital indices:
    vCapitalPathIndex(t) = mIndexCapital(vProdPathIndex(t), ...
        vCapitalPathIndex(t - 1));

end

% Obtain the levels of simulated productivity and capital:
vProdPath    = vProdGrid(vProdPathIndex);
vCapitalPath = vCapitalGrid(vCapitalPathIndex);



%--------------------------------------------------------------------------
% REPORT
%--------------------------------------------------------------------------

% Store all output in the corresponding folder:
results = strcat(directory, '/stochastic');
cd(results);

% Define indices for productivity:
low  = 1;
med  = median(0 : nProd);
high = nProd;

% Value function for different levels of productivity
f1 = figure;
plot(vCapitalGrid, mValue(low, :), 'LineWidth', 1.8);
hold on
plot(vCapitalGrid, mValue(med, :), 'LineWidth', 1.8);
plot(vCapitalGrid, mValue(high, :), 'LineWidth', 1.8);
hold off
xlabel('Current capital: $k$', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Value function: $V(k,z)$', 'interpreter', 'latex', 'FontSize', 14);
legend('Min productivity: $z=z_{1}$', 'Median productivity: $z=z_{10}$', 'Max productivity: $z=z_{20}$','location','Southeast', 'interpreter', 'latex')
title('\textbf{Value Function}', 'interpreter', 'latex','fontsize',16);
savefig("value_function")
saveas(gcf,'value_31.jpg')


% Policy function for different levels of productivity
f2 = figure;
plot(vCapitalGrid, mPolicyCapital(low, :), 'LineWidth', 1.8);
hold on
plot(vCapitalGrid, mPolicyCapital(med, :), 'LineWidth', 1.8);
hold on
plot(vCapitalGrid, mPolicyCapital(high, :), 'LineWidth', 1.8);
hold on
plot(vCapitalGrid, vCapitalGrid, 'k--', 'LineWidth', 0.5);
hold off
xlabel('Current capital: $k$', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Policy function: $k^{\prime}\equiv \phi(k,z)$', 'interpreter', 'latex', 'FontSize', 14);
legend('Min productivity: $z=z_{1}$', 'Median productivity: $z=z_{10}$', 'Max productivity: $z=z_{20}$','45o line','location','Southeast', 'interpreter', 'latex');
title('\textbf{Policy Function}', 'interpreter', 'latex','fontsize',16);
savefig("policy_function");
saveas(gcf,'policy_31.jpg');


% Simulation of productivity
f3 = figure;
plot(1 : nTime, exp(vProdPath), 'k--', 'LineWidth', 0.5);
xlabel('Time: $t$', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Simulated levels of productivity: $e^{z}$', 'interpreter', 'latex', 'FontSize', 14);
title('\textbf{Simulated path of productivity}', 'interpreter', 'latex','fontsize',16);
savefig("simulation_z")
saveas(gcf,'simulation_z.jpg')


% Histogram of capital
f4 = figure;
histogram(vCapitalPath, 25, 'Normalization', 'probability');
xlabel('Capital levels: $k$', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Probability', 'interpreter', 'latex', 'FontSize', 14);
title('\textbf{Histogram of capital levels}', 'interpreter', 'latex','fontsize',16);
savefig("histogram_k")
saveas(gcf,'histogram_k.jpg')
