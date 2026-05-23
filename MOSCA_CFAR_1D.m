function [detection, threshold, num_targets, Z, X_k, MeanW] = MOSCA_CFAR_1D(signal, numGuard, numTrain, P_fa, k_rank)
% MOSCA_CFAR_1D - Mean of Order Statistics and Cell Averaging CFAR
%
% Combines the advantages of OS-CFAR and CA-CFAR:
%   - OS-CFAR: Robust to interfering targets (uses k-th order statistic)
%   - CA-CFAR: Good performance in uniform clutter (uses averaging)
%
% Inputs:
%   signal        - Input signal (Nx1 vector, range profile)
%   numGuard      - Number of guard cells on each side
%   numTrain      - Number of training cells on each side (M = N = numTrain)
%   P_fa          - Desired false alarm rate (e.g., 1e-6)
%   k_rank        - Order statistic rank for leading window (optional)
%                   If not provided, uses 3/4 rule: k = 3M/4
%
% Outputs:
%   detection     - Binary detection mask (1 = target detected)
%   threshold     - CFAR threshold for each cell
%   num_targets   - Number of detected targets
%   Z             - Combined noise estimate (X_k + MeanW)
%   X_k           - k-th order statistic from leading window
%   MeanW         - Mean of lagging window (CA component)
%
% Reference: 
%   Yu et al., "Research and Implementation of Orderly Statistical 
%   Constant False Alarm Detector", IEEE 2019

%% Input validation and default parameters
if nargin < 5
    % Default: use 3/4 rule (75th percentile) for k
    M = numTrain;  % Leading window size
    k_rank = round(0.75 * M);
    fprintf('Using default k=%d (75%% of M=%d)\n', k_rank, M);
end

N = length(signal);
M = numTrain;  % Leading window cells
Nn = numTrain;  % Lagging window cells (named Nn to avoid conflict with N)

% Initialize outputs
detection = zeros(N, 1);
threshold = zeros(N, 1);
X_k = zeros(N, 1);      % OS component
MeanW = zeros(N, 1);    % CA component
Z = zeros(N, 1);        % Combined estimate

%% Calculate threshold multiplier T using Equation (3)
T = MOSCA_threshold_multiplier(M, Nn, k_rank, P_fa);

fprintf('\n=== MOSCA-CFAR Parameters ===\n');
fprintf('Leading window (M): %d cells\n', M);
fprintf('Lagging window (N): %d cells\n', Nn);
fprintf('Guard cells: %d on each side\n', numGuard);
fprintf('Order statistic: k = %d (%.1f%% of M)\n', k_rank, 100*k_rank/M);
fprintf('Threshold multiplier: T = %.4f\n', T);
fprintf('False alarm rate: P_fa = %.2e\n', P_fa);

% Calculate ADT (Average Decision Threshold) for performance assessment
ADT = T * (sum(1 ./ (M - (1:k_rank) + k_rank + 1)) + 1/(k_rank + 1));
fprintf('Average Decision Threshold: ADT = %.4f\n', ADT);
fprintf('============================\n\n');

%% Sliding window CFAR processing
for idx = M + numGuard + 1 : N - (Nn + numGuard)
    
    % ========== Leading Window (OS-CFAR Processing) ==========
    % Collect M training cells BEFORE the CUT
    leading_start = idx - M - numGuard;
    leading_end = idx - numGuard - 1;
    leading_cells = signal(leading_start:leading_end);
    
    % Sort and select k-th order statistic (k-th smallest value)
    sorted_leading = sort(leading_cells, 'ascend');
    X_k(idx) = sorted_leading(k_rank);
    
    % ========== Lagging Window (CA-CFAR Processing) ==========
    % Collect N training cells AFTER the CUT
    lagging_start = idx + numGuard + 1;
    lagging_end = idx + numGuard + Nn;
    lagging_cells = signal(lagging_start:lagging_end);
    
    % Calculate mean (average)
    MeanW(idx) = mean(lagging_cells);
    
    % ========== MOSCA Combination ==========
    % Z = X(k) + MeanW (selection logic from paper)
    Z(idx) = X_k(idx) + MeanW(idx);
    
    % Calculate threshold: S = T × Z
    threshold(idx) = T * Z(idx);
    
    % ========== Detection Test ==========
    if signal(idx) > threshold(idx)
        detection(idx) = 1;
    end
end
function T = MOSCA_threshold_multiplier(M, N, k, P_fa)
% Calculate threshold multiplier T for MOSCA-CFAR
% Based on Equation (3) from the paper:
%
% P_fa = (M choose k) × [Γ(M-k+1) × Γ(k+T)] / [Γ(M+T+1)] × (1 + 1/N)^T
%
% Inputs:
%   M     - Leading window cells
%   N     - Lagging window cells  
%   k     - Order statistic rank
%   P_fa  - Desired false alarm rate
%
% Output:
%   T     - Threshold multiplier

% This is a transcendental equation that requires numerical solution
% We'll use an approximation based on the paper's approach

% Method 1: Iterative solution (accurate but slower)
if M <= 32 && N <= 32
    T = solve_T_iterative(M, N, k, P_fa);
else
    % Method 2: Approximation (faster for large windows)
    T = solve_T_approximation(M, N, k, P_fa);
end

% Ensure T is positive and reasonable
T = max(T, 1);

end

function T = solve_T_iterative(M, N, k, P_fa)
% Iterative solution using Newton-Raphson or binary search

% Initial guess based on OS-CFAR formula
T_guess = M * (P_fa^(-1/k) - 1);

% Binary search for T
T_min = 0.1;
T_max = 100;
tolerance = 1e-6;
max_iter = 100;

for iter = 1:max_iter
    T_mid = (T_min + T_max) / 2;
    
    % Calculate P_fa for current T
    P_calc = calculate_MOSCA_Pfa(M, N, k, T_mid);
    
    if abs(P_calc - P_fa) < tolerance
        T = T_mid;
        return;
    end
    
    if P_calc > P_fa
        T_min = T_mid;  % Increase T to decrease P_fa
    else
        T_max = T_mid;  % Decrease T to increase P_fa
    end
end

T = (T_min + T_max) / 2;

end

function P_fa = calculate_MOSCA_Pfa(M, N, k, T)
% Calculate false alarm probability for given parameters
% Based on Equation (3)

try
    % Binomial coefficient: M choose k
    binom_coeff = nchoosek(M, k);
    
    % Gamma function terms
    gamma_term = (gamma(M - k + 1) * gamma(k + T)) / gamma(M + T + 1);
    
    % Lagging window term
    lag_term = (1 + 1/N)^T;
    
    % Combined
    P_fa = binom_coeff * gamma_term * lag_term;
    
catch
    % If numerical overflow, use log-space calculation
    log_binom = sum(log(1:M)) - sum(log(1:k)) - sum(log(1:(M-k)));
    log_gamma = gammaln(M - k + 1) + gammaln(k + T) - gammaln(M + T + 1);
    log_lag = T * log(1 + 1/N);
    
    P_fa = exp(log_binom + log_gamma + log_lag);
end

% Ensure valid probability
P_fa = min(max(P_fa, 0), 1);

end
%%
function T = solve_T_approximation(M, N, k, P_fa)
% Fast approximation for large window sizes
% Based on simplified analysis

% Approximate using weighted combination of OS and CA formulas
T_os = M * (P_fa^(-1/k) - 1);
T_ca = M * (P_fa^(-1/M) - 1);

% Weight factor (more weight on OS component)
weight = k / M;  % Higher k → more weight on OS

T = weight * T_os + (1 - weight) * T_ca;

% Apply correction factor from paper's empirical data
% From Table 1: T_mosca ≈ 0.7 × T_os for typical parameters
correction = 0.7;
T = correction * T;

end
%% Post-processing: Remove duplicate detections
detection = remove_duplicates_1D(detection, numGuard);
num_targets = sum(detection);

fprintf('Detection complete: %d targets found\n', num_targets);
function detection_clean = remove_duplicates_1D(detection, window_size)
% Remove duplicate detections within a local window
% Keep only the first detection in each cluster

detection_clean = detection;
detect_idx = find(detection);

if isempty(detect_idx)
    return;
end

% Group nearby detections
i = 1;
while i < length(detect_idx)
    if i < length(detect_idx) && (detect_idx(i+1) - detect_idx(i)) <= 2*window_size
        % Remove duplicate (keep first)
        detection_clean(detect_idx(i+1)) = 0;
        detect_idx(i+1) = [];
    else
        i = i + 1;
    end
end

end

end