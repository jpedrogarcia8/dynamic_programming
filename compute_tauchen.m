%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%          REPLICATION OF THE HOPENHAYN & ROGERSON MODEL (1993 JPE)

%              José Pedro Garcia - European University Institute
%                              October, 2025

%                          IMPLEMENT TAUCHEN METHOD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function implements the Tauchen method which discretises a
% stochastic productivity process of the following type: 

% log(z_t) = \rho * log(z_t-1) + (1 - rho) * log(z0) + \sigma * epsilon_t,
% where epsilon ~ N(0, 1)

% INPUTS:
% 1) rrho    --> Serial correlation of productivity 
% 2) ssigma  --> Standard deviation of productivity innovations 
% 3) mmu0    --> Mean of level TFP of incumbent firms
% 4) nProd   --> Number of grid points 
% 5) kappa   --> Lower bound for the grid

% OUTPUT
% 1) vProdGrid       --> A vector of size "nProd" with the grid point for
% productivity
% 2) mProdTransition --> A transition matrix of dimension (nProd x nProd)
% and whose generic element (i,j) gives the probability of producitivity j
% tomorrow, given productivity i today. This means that each row should sum
% up to 1.
% 3) vProdErgodic    --> The ergodic/stationary distribution of
% productivity


function [vProdGrid, mProdTransition, vProdErgodic] = ...
    compute_tauchen(rrho, ssigma, mmu0, nProd, kappa)


% Define the standard deviation of the log TFP:
ssigmaTFP = ssigma / sqrt(1 - rrho ^ 2);

% Set up grid for productivity:
vProdRaw  = linspace(1, nProd, nProd)';
ddelta    = 2 * kappa * ssigmaTFP / (nProd - 1);                           % Distance b/w grid points

vProdGrid = log(mmu0) - (kappa * ssigmaTFP) + ddelta * (vProdRaw - 1);     % Grid for log(z)

 

% -------------------------------------------------------------------------
% Set up matrix for the transition probabilities
% -------------------------------------------------------------------------

% Set up matrices to deal with top and bottom values of the productivity
% grid:
vPlus             = zeros(1, nProd);
vPlus(1, nProd)   = 1e9;
vPlus(1, 1:end-1) = vProdGrid(2:end)';

vMinus            = zeros(1, nProd);
vMinus(1,1)       = -vPlus(1, nProd);
vMinus(1,2:end)   = vProdGrid(1:end-1)';

mPlusCutoff       = (0.5 * ones(nProd, 1) * (vProdGrid' + vPlus) - ...
    rrho * vProdGrid * ones(1, nProd) - (1 - rrho) * log(mmu0)) / ssigma;
mMinusCutoff      = (0.5 * ones(nProd, 1) * (vProdGrid' + vMinus) - ...
    rrho * vProdGrid * ones(1, nProd) - (1 - rrho) * log(mmu0)) / ssigma;

% The innovations follow a standard normal distribution, so we can just
% rely on the corresponding CDF:
mProdTransition   = normcdf(mPlusCutoff, 0, 1) - ...
    normcdf(mMinusCutoff, 0, 1);

% Sanity check: verify that each row sums up to 1:
vSanity           = (sum(mProdTransition') == ones(1, nProd))';



% -------------------------------------------------------------------------
% Compute ergodic distribution of productivity
% -------------------------------------------------------------------------

% The ergodic distribution of the process (call it \nu) solves the
% following problem: \nu' = \nu' * F. This amounts to solving the left
% eigenvectors of the probability transition matrix, F, associated with the
% unit eigenvalue.

[~, mEigenValues, mLeftEigenVectors] = eig(mProdTransition);

% Identify the column where the unit eigenvalue is located:
vEigenValues                         = round(mEigenValues(:), 1);
vEigenValues(vEigenValues == 0)      = [];
unitColumn                           = find(vEigenValues == 1);            % Number of unit eigenvalue

% Ergodic distribution corresponds to the eigenvector in column
% "unitColumn" of the mLeftEigenVectors matrix:
vProdErgodic                         = mLeftEigenVectors(:,unitColumn);
vProdErgodic                         = vProdErgodic / sum(vProdErgodic);   % Normalise such that it sums up to 1


end