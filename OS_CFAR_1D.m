function [detection, threshold, k_detect] = OS_CFAR_1D(signal, numGuard, numTrain, P_fa, k)
    % OS-CFAR 1D implementation using efficient sliding window
    % signal: Input signal (power or magnitude squared)
    % numGuard: Number of guard cells on each side
    % numTrain: Number of training cells on each side
    % P_fa: Desired probability of false alarm
    % k: Order statistic index (typically 0.75*2*numTrain for robustness)
    
    % Ensure signal is a row vector
    signal = signal(:)';
    N = length(signal);
    total_train = 2 * numTrain;
    
    % Validate k parameter
    if k < 1 || k > total_train
        error('k must be between 1 and %d (2*numTrain)', total_train);
    end
    
    % Calculate alpha for OS-CFAR using the given formula
    % P_fa = prod_{i=0}^{k-1} (N-i)/(alpha+N-i)
    alpha = calculate_alpha_OS(total_train, k, P_fa);
    
    % Initialize threshold array
    threshold = zeros(1, N);
    noise_level = zeros(1, N);
    
    % Window parameters
    half_window = numTrain + numGuard;
    
    % Efficient vectorized approach using sliding window
    % Pre-allocate training cell matrix
    for idx = 1:N
        % Define training cell indices (left and right windows)
        left_start = max(1, idx - half_window);
        left_end = max(1, idx - numGuard - 1);
        right_start = min(N, idx + numGuard + 1);
        right_end = min(N, idx + half_window);
        
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
        
        % Sort and select k-th order statistic
        if ~isempty(train_cells)
            sorted_cells = sort(train_cells, 'ascend');
            k_actual = min(k, length(sorted_cells));
            noise_level(idx) = sorted_cells(k_actual);
        else
            noise_level(idx) = 0;
        end
    end
    
    % Calculate threshold
    threshold = alpha * noise_level;
    alpha
    % Detection decision
    detections = signal > threshold;
    
    % Remove duplicate detections (keep strongest in local window)
    detection = remove_duplicates_1D(detections, signal, numGuard);
    k_detect = sum(detection);
end


function alpha = calculate_alpha_OS(N, k, P_fa)
    % Calculate alpha for OS-CFAR using numerical solution
    % P_fa = prod_{i=0}^{k-1} (N-i)/(alpha+N-i)
    
    % Define the equation to solve: f(alpha) = 0
    equation = @(alpha) compute_Pfa(alpha, N, k) - P_fa;
    
    % Solve using fzero with reasonable bounds
    options = optimset('TolX', 1e-10, 'Display', 'off');
    try
        alpha = fzero(equation, [0.01, 10000], options);
    catch
        % If fzero fails
        warning('fzero failed, using approximate alpha');
        alpha = N * (P_fa^(-1/k) - 1); % Approximation
    end
end



function P_fa_calc = compute_Pfa(alpha, N, k)
    % Compute P_fa for given alpha
    % P_fa = prod_{i=0}^{k-1} (N-i)/(alpha+N-i)
    
    P_fa_calc = 1;
    for i = 0:(k-1)
        P_fa_calc = P_fa_calc * (N - i) / (alpha + N - i);
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


