function [detection, threshold, k] = CFAR_1D_CA_(signal, numGuard, numTrain, P_fa)
    % CA-CFAR 1D implementation using convolution
    
    N = length(signal);
    
    % Calculate CFAR parameter (alpha)
    total_train = 2 * numTrain;
    alpha = total_train * (P_fa^(-1/total_train) - 1)
    
    % Create CFAR kernel
    kernel_length = 2 * (numTrain + numGuard) + 1;
    cfar_kernel = ones(1, kernel_length) / total_train;
    
    % Zero out guard cells and CUT
    center = numTrain + numGuard + 1;
    cfar_kernel(center - numGuard : center + numGuard) = 0;
    
    % Compute noise level using convolution
    noise_floor = conv(signal, cfar_kernel, 'same');
    
    % Calculate threshold
    threshold = alpha * noise_floor;
    
    % Detection decision
    detections = signal > threshold;
    
    % Remove duplicate detections
    detection = remove_duplicates_1D(detections, numGuard);
    k = sum(detection);
    
end

%% remaining part is for target location estimation
function detection_clean = remove_duplicates_1D(detection, numGuard)
% Remove duplicate detections within a local window
% Keep only the local maximum
detection_clean = detection;
detect_idx = find(detection);
if isempty(detect_idx)
return;
end
% Group nearby detections
i = 1;
while i < length(detect_idx)
% Check if next detection is within window
    if (detect_idx(i+1) - detect_idx(i)) <= 2*numGuard
% Mark second detection as duplicate
        detection_clean(detect_idx(i+1)) = 0;
        detect_idx(i+1) = [];
    else
        i = i + 1;
    end
end
end
