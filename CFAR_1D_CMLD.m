function [detection, threshold, k_detect] = CFAR_1D_CMLD(signal, numGuard, numTrain, P_fa, k_censor)
    % CMLD-CFAR (Censored Mean Level Detector) 1D implementation
    % 
    % signal: Input signal (power or magnitude squared)
    % numGuard: Number of guard cells on each side
    % numTrain: Number of training cells on each side
    % P_fa: Desired probability of false alarm
    % k_censor: Number of largest cells to censor (remove)
    %           Typical values: k = 0.25*N to 0.5*N
    
    % Ensure signal is a row vector
    signal = signal(:)';
    N_sig = length(signal);
    
    % Total training cells
    N = 2 * numTrain;
    
    % Validate k_censor parameter
    if k_censor < 0 || k_censor >= N
        error('k_censor must be between 0 and %d (< 2*numTrain)', N-1);
    end
    
    % Calculate threshold factor T for CMLD-CFAR
    T = calculate_T_CMLD(N, k_censor, P_fa);
    
    % Initialize arrays
    threshold = zeros(1, N_sig);
    noise_level = zeros(1, N_sig);
    
    % Window parameters
    half_window = numTrain + numGuard;
    
    % Sliding window CMLD-CFAR
    for idx = 1:N_sig
        % Define training cell indices (left and right windows)
        left_start = max(1, idx - half_window);
        left_end = max(1, idx - numGuard - 1);
        right_start = min(N_sig, idx + numGuard + 1);
        right_end = min(N_sig, idx + half_window);
        
        % Collect training cells from both sides
        if left_end >= left_start && right_end >= right_start
            train_cells = [signal(left_start:left_end), signal(right_start:right_end)];
        elseif left_end >= left_start
            train_cells = signal(left_start:left_end);
        elseif right_end >= right_start
            train_cells = signal(right_start:right_end);
        else
            train_cells = [];
        end
        
        % CMLD Processing: Sort and censor k largest cells
        if ~isempty(train_cells)
            N_actual = length(train_cells);
            k_actual = min(k_censor, N_actual - 1); % Need at least 1 cell
            
            % Sort in ascending order
            sorted_cells = sort(train_cells, 'ascend');
            
            % Censor k largest cells (keep only smallest N-k cells)
            censored_cells = sorted_cells(1:end-k_actual);
            
            % Average the censored (remaining) cells
            noise_level(idx) = mean(censored_cells);
        else
            noise_level(idx) = 0;
        end
    end
    
    % Calculate threshold: threshold = T × noise_level
    threshold = T * noise_level;
    
    % Detection decision
    detections = signal > threshold;
    
    % Remove duplicate detections (keep strongest in local window)
    detection = remove_duplicates_1D(detections, signal, numGuard);
    k_detect = sum(detection);
end


function T = calculate_T_CMLD(N, k, P_fa)
    % Calculate threshold factor T for CMLD-CFAR
    % P_fa = C(N,k) * prod_{j=1}^{k} [T + (N-j+1)/(k-j+1)]^(-1)
    %
    % Where C(N,k) is the binomial coefficient
    
    % Binomial coefficient C(N, k)
    C_Nk = nchoosek(N, k);
    
    % Define the equation to solve: f(T) = 0
    % Rearrange: prod_{j=1}^{k} [T + (N-j+1)/(k-j+1)] = C(N,k) / P_fa
    
    target_prod = C_Nk / P_fa;
    
    equation = @(T) compute_product_CMLD(T, N, k) - target_prod;
    
    % Solve using fzero with reasonable bounds
    options = optimset('TolX', 1e-10, 'Display', 'off');
    
    try
        % Try to find T in reasonable range
        T = fzero(equation, [0.01, 1000], options);
    catch
        % If fzero fails, use approximation
        warning('fzero failed, using approximate T');
        % Approximation: T ≈ (N-k) * (C(N,k)/P_fa)^(1/k) - average term
        avg_term = 0;
        for j = 1:k
            avg_term = avg_term + (N - j + 1) / (k - j + 1);
        end
        avg_term = avg_term / k;
        T = (target_prod)^(1/k) - avg_term;
        T = max(T, 0.1); % Ensure positive
    end
end


function prod_val = compute_product_CMLD(T, N, k)
    % Compute the product term for CMLD-CFAR
    % prod_{j=1}^{k} [T + (N-j+1)/(k-j+1)]
    
    prod_val = 1;
    for j = 1:k
        term = T + (N - j + 1) / (k - j + 1);
        prod_val = prod_val * term;
    end
end


function detection_clean = remove_duplicates_1D(detection, signal, numGuard)
    % Remove duplicate detections within a local window
    % Keep only the local maximum (strongest signal)
    
    detection_clean = detection;
    detect_idx = find(detection);
    
    if isempty(detect_idx)
        return;
    end
    
    % Group nearby detections and keep strongest
    i = 1;
    while i <= length(detect_idx)
        % Find all detections within the window
        window_start = detect_idx(i);
        j = i + 1;
        
        % Collect all indices within 2*numGuard distance
        window_indices = detect_idx(i);
        while j <= length(detect_idx) && (detect_idx(j) - detect_idx(i)) <= 2*numGuard
            window_indices = [window_indices, detect_idx(j)];
            j = j + 1;
        end
        
        % If multiple detections in window, keep only the strongest
        if length(window_indices) > 1
            [~, max_local_idx] = max(signal(window_indices));
            strongest_idx = window_indices(max_local_idx);
            
            % Zero out all except the strongest
            for idx = window_indices
                if idx ~= strongest_idx
                    detection_clean(idx) = 0;
                end
            end
        end
        
        % Move to next group
        i = j;
    end
end