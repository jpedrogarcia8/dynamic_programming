%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               MACROECONOMICS III - DYNAMIC PROGRAMMING
%%                           Problem Set 2

%                          José Pedro Garcia

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script is meant to develop a code which solves the neoclassical
% growth model *without* uncertainty: exercise 1 of the problem set.


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

% The steady-state level of capital is given as follows (see notes):
ssCapital = ((aalpha * bbeta) / (1 - bbeta * (1 - ddelta))) ^ ...
    (1 / (1 - aalpha));

% Set numerical parameters and grid for capital:
nCapital  = 500;                                                           % Grid points for capital
vCapital  = linspace(0.9 * ssCapital, 1.1 * ssCapital, nCapital)';         % Grid for capital

% Define grids for state and choice sets:
% CONVENTION:
% ROWS    = State capital 
% COLUMNS = Choice capital
nState         = nCapital;                                                 % Number of possible realisations of the state vector
nChoice        = nCapital;                                                 % Number of possible choices

mStateCapital  = repmat(vCapital, 1, nChoice);
mChoiceCapital = repmat(vCapital', nState, 1); 

% Compute all possible consumption levels:
mConsumption   = (mStateCapital .^ aalpha) + mStateCapital * (1 - ddelta) - ...
    mChoiceCapital;



%--------------------------------------------------------------------------
% VALUE FUNCTION ITERATION
%--------------------------------------------------------------------------

% Pre-allocate necessary objects:
vValue        = zeros(nState, 1);                                          % Initial guess for the value function
dampening     = 0;                                                         % Weight placed on previous iteration
toleranceVFI  = 1e-10;                                                     % Tolerance allowed for convergence
iterationVFI  = 1;                                                         % To keep track of number of iterations
maxIterations = 1e4;                                                       % Maximum number of iterations allowed
difference    = 100;                                                       % To kick off iteration
vDifference   = zeros(maxIterations, 1);                                   % Record the convergence for easier debugging

tic;
while abs(difference) > toleranceVFI && iterationVFI <= maxIterations

    % Compute continuation value:
    mContinuationValue         = repmat(vValue', nState, 1);

    % Compute objective function:
    mObjectiveFunction         = log(mConsumption) + bbeta * mContinuationValue;

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

% Retrieve capital policy function:
vPolicyCapital = vCapital(vCapitalIndex);



%--------------------------------------------------------------------------
% SIMULATION OF CAPITAL ACCUMULATION THROUGHOUT TIME
%--------------------------------------------------------------------------

% Pre-allocate necessary objects:
nTime                = 1e2;                                                % Number of time periods to be simulated
vCapitalPathIndex    = zeros(nTime, 1);                                    % Vector for simulation
vCapitalPathIndex(1) = 1;                                                  % Say the economy starts at the lowest value of capital (first in the grid)

for t = 2 : nTime
    vCapitalPathIndex(t) = vCapitalIndex(vCapitalPathIndex(t - 1));
end

% Obtain the levels of capital:
vCapitalPath = vPolicyCapital(vCapitalPathIndex);



%--------------------------------------------------------------------------
% REPORT
%--------------------------------------------------------------------------

% Store all output in the corresponding folder:
results = strcat(directory, '/deterministic');
cd(results);

% Value function as a function of current capital:
f1 = figure;
plot(vCapital, vValue,'LineWidth',1.5); 
title('$\textbf{Value function}$','Interpreter','latex','fontsize',16);
xlabel('Current capital: $k$','Interpreter','latex','fontsize',14);   
ylabel('Value function: $V(k)$','Interpreter','latex','fontsize',14);
savefig("value_function");
saveas(gcf,'value.jpg');

% Policy function for k'
f2 = figure;
plot(vCapital, vPolicyCapital, 'LineWidth', 1.8);
hold on
plot(vCapital, vCapital, 'k--', 'LineWidth', 0.5);
hold off
xlabel('Current capital: $k$', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Future capital: $k^{\prime}$', 'interpreter', 'latex', 'FontSize', 14);
legend('Policy function', '45º line','location','Southeast');
title('\textbf{Policy Function}', 'interpreter', 'latex','fontsize',16);
savefig("policy_function");
saveas(gcf,'policy.jpg')

% Policy function for savings; k' - k
vPolicySavings = vPolicyCapital - vCapital;

f3 = figure;
plot(vCapital, vPolicySavings, 'LineWidth', 1.8); 
yline(0, 'k--','LineWidth',0.8);
title('$\textbf{Savings function}$','Interpreter','latex','fontsize',16);
xlabel('Current capital: $k$','Interpreter','latex','fontsize',14);   
ylabel('Savings: $k{\prime} - k$','Interpreter','latex','fontsize',14);
savefig("savings_function");
saveas(gcf,'savings.jpg')

% Simulation of capital accumulation
f4 = figure;
plot(1 : nTime, vCapitalPath, 'LineWidth', 1.8) 
title('$\textbf{Capital accumulation}$','Interpreter','latex','fontsize',16)
xlabel('Time period: $t$','Interpreter','latex','fontsize',14)   
ylabel('Capital: $k^{\prime}$','Interpreter','latex','fontsize',14)
savefig("simulation")
saveas(gcf,'simulation.jpg')
