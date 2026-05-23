
% 2Tx-4Rx processing 
% Most recently updated on - 20/08/2025
% 11:22
%%
clc;
close all;
clear all;

%% ---------------- Radar & Waveform Parameters ----------------
fc = 77.5e9;                  % Carrier frequency (Hz)
c = 3e8;                      % Speed of light (m/s)
lambda = c/fc;                % Wavelength (m)
BW = 2.04e9;                  % Bandwidth of chirp (Hz)
Tchirp = 45e-6;               % Chirp duration (s)
S = BW/Tchirp;                % Slope of chirp (Hz/s)
fs = 10e6;                    % ADC sampling rate (Hz)
Num_adc_samples = 256;        % Number of ADC samples per chirp
Num_chirps = 400;%48000;  %40;            % Total chirps 
Num_subframes = Num_chirps/2; % Number of subframes (each has 2 chirps due to BPM)
% -- Antenna geometry -- ---
numTx = 2;                    % number of Transmitting antennas
numRx = 4;                    % number of Receiving antennas
% TxPos = [0, 2*lambda];      % Tx1 at 0m, Tx2 at 2λ (2λ spacing)
% RxPos = (0:3)*(lambda/2);   % 4 Rx at 0, λ/2, λ, 3λ/2
%% Converting .bin file into matlab data file
% Reading data from .bin file

file_list = {% do label the files
    %'adc_data_5balls_v2.bin','5 balls'                                      %adc_data_Empty_NF_240chirps_100frames_v2.
    %'adc_data_2nuts_stacked_8metre_0Deg.bin', '8m @ 0°';                    %adc_data_Multichannel_Empty_NF_2.bin   
    % 'adc_data_2nuts_stacked_thicker_3metres_0Deg.bin', '3m @ 0°';                 %adc_data_Multichannel_Empty_NF_2.bin
    % 'adc_data_2nuts_4metre_3metre_0_30Deg.bin', '4m & 3m @ 0° & 30°'
    %'adc_data_3Reflectors_4_6_8metres.bin','3 reflectors 4,6,8 m'           %adc_data_Empty_NF_100chirps
    'adc_data_5balls_200chirps_1frame.bin','5 metal balls'                        %adc_data_Empty_NF_100chirps
};

%% clutter map with no target,
%raw_data_NF = readca1000('adc_data_Empty_NF_240chirps_100frames_v2.bin');
%raw_data_NF = readca1000('adc_data_Multichannel_Empty_NF_2.bin');
raw_data_NF = readca1000('adc_data_Empty_NF_100chirps.bin');
% Reshape into [samples x chirps x rx]
adc_data_NF = zeros(Num_adc_samples, Num_chirps, numRx);
for rx = 1:numRx
    adc_data_NF(:, :, rx) = reshape(raw_data_NF(rx, :), Num_adc_samples, Num_chirps);
end
%%  Load ADC data
for file_idx = 1:size(file_list, 1)
    filename = file_list{file_idx, 1};
    file_label = file_list{file_idx, 2};
% Reshape into [samples x chirps x rx]
     
    raw_data = readca1000(filename);
    adc_data = zeros(Num_adc_samples, Num_chirps, numRx);
    for rx = 1:numRx
        adc_data(:, :, rx) = reshape(raw_data(rx, :), Num_adc_samples, Num_chirps);
    end
    
    
    %% BPM Decoding
    %  BPM Code
    bpm = [1  1;               % Chirp 1: Tx1 + Tx2
           1 -1];              % Chirp 2: Tx1 - Tx2
    decoding_matrix = 0.5 * [1 1; 1 -1];  % Inverse of BPM matrix for decoding
    
    for rx = 1:numRx
        for subframe = 1:Num_subframes
            chirp_index = (subframe-1)*2 + (1:2);
            raw_data_temp = adc_data(:,chirp_index,rx);
            decoded = raw_data_temp * decoding_matrix.';
            raw_data_NF_temp = adc_data_NF(:,chirp_index,rx);
            decoded_NF = raw_data_NF_temp * decoding_matrix.'; % for background subtraction
      
            for tx=1:2
                decoded_block(:,subframe,tx,rx) = decoded(:,tx);
                decoded_block_NF(:,subframe,tx,rx) = decoded_NF(:,tx);
            end
        end
    end
    %% Background subtraction
    decoded_block_postsub = decoded_block - decoded_block_NF;
    % decoded_block_postsub_limited = decoded_block_postsub(1:30,:,:,:);
    %% Range FFT -> seems to be right .. :)
    % window = hamming(Num_adc_samples);
    % Taking 1D-FFT across 256- adc samples for 20 chirps, 2tx , 4rx 
    range_fft = fft(decoded_block_postsub, Num_adc_samples, 1);  % [range x chirp x tx x rx]
    % Taking Average of all 20x2x4 FFT values 
    range_fft_avg = mean(range_fft, [2 3 4]); 
    % Range axis 
    %range_Resolution=c/(2*BW)
    range_axis = (0:Num_adc_samples-1) * (c/(2*BW));
    % Plotting the average range FFT obtained -> peaks correspond to targets
    % fig_range = figure('Name', sprintf('File %d: Range Profile - %s', file_idx, file_label));
    % plot(range_axis, 10*log10(abs(range_fft_avg)), 'b');
    % xlabel('Range (m)');
    % ylabel('Amplitude (dB)');
    % title(sprintf('Range Profile - %s', file_label));
    % grid on;
    %% Doppler FFT
    % % Doppler FFT across subframes (slow time)
    % doppler_res = lambda/(2*Num_subframes*Tchirp);
    % % Taking FFT along slow-time axis 
    % doppler_fft = fftshift(fft(range_fft, Num_subframes, 2), 2);  % [range × doppler × tx × rx]
    % doppler_axis = (-Num_subframes/2 : Num_subframes/2 - 1) *doppler_res ;  % optional, Hz or m/s
    % % Taking Average of all  
    % % doppler_fft_avg = mean((fftshift(fft(range_fft, Num_subframes, 2), 2)),[1 3 4]); 
    % % Plotting
    % Plot Range-Doppler Map
    % fig_rd = figure('Name', sprintf('File %d: Range-Doppler - %s', file_idx, file_label));
    % imagesc(doppler_axis, range_axis, 10*log10(abs(doppler_fft(:,:,1,1))));
    % xlabel('Velocity (m/s)');
    % ylabel('Range (m)');
    % title(sprintf('Range-Doppler Heatmap - %s', file_label));
    % colormap jet; colorbar;
    % axis xy;
    
    
    %% 2D CA-CFAR
    % % CFAR parameters
    % numTrain = 8; % # no. of training cells
    % numGuard = 1; % # no. of guard/gap cells
    % P_far = 2*1e-1; % desired false alarm rate 
    % SNR_OFFSET = -5 ;% dB
    % %RDM_dB = 10*log10(abs(doppler_fft(:,:,1,2))/max(max(abs(doppler_fft(:,:,1,2))))); % Convert RDM to dB scale ( 2TX 3RX )
    % RDM = (mean(abs(doppler_fft), [3 4]));
    % %RDM = (abs(doppler_fft(:,:,1,2)));
    % RDM_dB=10*log10(abs(doppler_fft(:,:,2,3)));
    % [RDM_mask, cfar_ranges, cfar_dopps, k] = CA_CFAR(RDM, numGuard, numTrain, P_far, SNR_OFFSET);
    % disp(['no 0f targets idtenified:', num2str(k)]);
    % figure;
    % h=imagesc(doppler_axis,range_axis,RDM_mask);
    % xlabel('Velocity (m/s)')
    % ylabel('Range (m)')
    % title('CA-CFAR')
    % 
    
    % CFAR_1D_CA
    % CFAR parameters
%     numTrain = 8; % # no. of training cells one side4
%     numGuard = 3; % # no. of guard/gap cells one side4
%     P_far = 1e-5; %4.35 desired false alarm rate 1e-3
% 
%     avg_signal = (abs(range_fft_avg)).^2; %256x1  
%     %signal=abs(range_fft(:,1,1,1)); %256x1
%     [detection_CA, thresh_CA, k] = CFAR_1D_CA_(avg_signal,numGuard,numTrain,P_far);
%     disp(['no 0f targets idtenified:', num2str(k)]);
% 
%     fig_CA = figure('Name', sprintf('File %d: CA-CFAR Range - %s', file_idx, file_label));
%     plot(range_axis,10*log10(avg_signal),'b'); hold on;
%     plot(range_axis,10*log10(thresh_CA),'r');
%     plot(range_axis(detection_CA==1), 10*log10(avg_signal(detection_CA==1)), 'ko', 'MarkerFaceColor', 'g');
% 
%     % Add parameter text box
%     param_text = sprintf('Train Cells: %d  Guard Cells: %d  P_{fa} %0.0e',numTrain, numGuard, P_far);
%     text(0.02, 0.98, param_text, 'Units', 'normalized', ...
%         'VerticalAlignment', 'top', 'BackgroundColor', 'white', ...
%         'EdgeColor', 'black', 'FontSize', 10);
% 
%     xlabel('Range (m)');
%     ylabel('Amplitude in dB');
%     legend('Signal', 'Threshold', 'Detections');
%     title(sprintf('CA-CFAR Detection - %s', file_label));
%     % --- Save the figure ---
% filename = sprintf('CA_CFAR_%s_T%d_G%d_Pfa%.1e.png', ...
%                    file_label, numTrain, numGuard, P_far);  % dynamic name
% exportgraphics(gcf, fullfile('C:\Users\Satyam Pandey\OneDrive\Pictures\matlab saved', filename), 'Resolution', 300);      % high-quality PNG

    %% OS-CFAR Example Usage
% CFAR parameters
% numTrain = 8;   % Number of training cells on each side
% numGuard = 2;  % Number of guard cells on each side
% P_fa = 1e-5;    % Desired false alarm rate
% k_percentile=75;
% k = floor(((2*numTrain)/100)*k_percentile);    % Order statistic (k-th smallest out of 2*numTrain)
% % Your signal
% avg_signal = (abs(range_fft_avg)).^2; % 256x1
% 
% % Run OS-CFAR
% [detection_OS, thresh_OS, num_targets] = OS_CFAR_1D(avg_signal, numGuard, numTrain, P_fa, k);
% disp(['Number of targets identified: ', num2str(num_targets)]);
% 
% % Plot results
% fig_OS = figure('Name', sprintf('File %d: OS-CFAR Range - %s', file_idx, file_label));
% plot(range_axis, 10*log10(avg_signal), 'b'); hold on;
% plot(range_axis, 10*log10(thresh_OS), 'r');
% plot(range_axis(detection_OS==1), 10*log10(avg_signal(detection_OS==1)), 'ko', 'MarkerFaceColor', 'g');
% 
% % Add parameter text box
% param_text = sprintf('Train: %d  Guard: %d  P_{fa}: %.0e  k: %d', numTrain, numGuard, P_fa, k);
% text(0.02, 0.98, param_text, 'Units', 'normalized', ...
%     'VerticalAlignment', 'top', 'BackgroundColor', 'white', ...
%     'EdgeColor', 'black', 'FontSize', 10);
% 
% xlabel('Range (m)');
% ylabel('Amplitude (dB)');
% legend('Signal', 'Threshold', 'Detections');
% title(sprintf('OS-CFAR Detection - %s', file_label));
%     % --- Save the figure ---
% filename = sprintf('OS_CFAR_%s_K%d_T%d_G%d_Pfa%.1e.png', ...
%                    file_label, k, numTrain, numGuard, P_far);
% 
% exportgraphics(gcf, fullfile('C:\Users\Satyam Pandey\OneDrive\Pictures\matlab saved', filename), 'Resolution', 300);    % high-quality PNG
% 

    %% CMLD-CFAR Parameters
numTrain = 10;      % Training cells per side
numGuard = 7;     % Guard cells per side
P_far = 1e-5;       % False alarm rate
k=48 ; % top eg.(40%) are censored
k_censor = ceil(((numTrain*2)/100)*k);  % Censor 3 largest cells (out of 10 total)
                                        

% Your signal
avg_signal = (abs(range_fft_avg)).^2;

% Run CMLD-CFAR
[detection_CMLD, thresh_CMLD, num_targets] = CFAR_1D_CMLD(avg_signal, numGuard, numTrain, P_far, k_censor);
disp(['CMLD:Number of targets identified: ', num2str(num_targets)]);

% Plot results
figure('Name', 'CMLD-CFAR Detection');
plot(range_axis, 10*log10(avg_signal), 'b'); hold on;
plot(range_axis, 10*log10(thresh_CMLD), 'r');
plot(range_axis(detection_CMLD==1), 10*log10(avg_signal(detection_CMLD==1)), ...
     'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

% Add parameter text box
param_text = sprintf('Train: %d  Guard: %d  P_{fa}: %.0e  Censored: %d/%d', ...
                     numTrain, numGuard, P_far, k_censor, 2*numTrain);
text(0.02, 0.98, param_text, 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'BackgroundColor', 'white', ...
     'EdgeColor', 'black', 'FontSize', 10);

xlabel('Range (m)');
ylabel('Amplitude (dB)');
legend('Signal', 'Threshold', 'Detections', 'Location', 'best');
title(sprintf('CMLD-CFAR Detection - %s', file_label))
grid on;
% --- Save the figure ---
filename = sprintf('CMLD_CFAR_%s_K%d_T%d_G%d_Pfa%.1e.png', ...
                   file_label, k_censor, numTrain, numGuard, P_far);   % dynamic name
exportgraphics(gcf, fullfile('C:\Users\Satyam Pandey\OneDrive\Pictures\matlab saved', filename), 'Resolution', 300);      % high-quality PNG

    %% Angle of arrival estimation
    % N_az = 180;                             % Number of angle bins
    % angle_axis = linspace(-90, 90, N_az);   % Angle  in degrees
    % lambda = c / fc;                        % Wavelength
    % d = lambda / 2;                         % Element spacing (half wavelength)
    % num_antennas = numTx * numRx;           % Virtual antenna count
    % ant_pos = (0:num_antennas-1) * d;       % Antenna positions (ULA)
    % angle_map = zeros(Num_adc_samples, Num_subframes, N_az);
    % 
    % % Range
    % for r_bin = 1:Num_adc_samples
    %     % Doppler
    %     for d_bin = 1:Num_subframes
    %         x = zeros(num_antennas, 1);  
    %         ant_idx = 1;
    %         % Tx antenna
    %         for tx = 1:numTx
    %             % Rx antenna
    %             for rx = 1:numRx
    %                 x(ant_idx) = doppler_fft(r_bin, d_bin, tx, rx);
    %                 ant_idx = ant_idx + 1; 
    %             end
    %         end
    %         % Loop over angle 
    %         for i = 1:N_az
    %             theta = angle_axis(i);
    %             steering_vec = exp(1j * 2 * pi / lambda * ant_pos.' * sind(theta));
    %             angle_map(r_bin, d_bin, i) = abs(steering_vec' * x);  % Output 
    %         end
    %     end
    % end
    % 
    % % Average over range and doppler for AOA spectrum
    % angle_spectrum = squeeze(mean(mean(angle_map, 1), 2));
    % % Plot AoA spectrum
    % figure;
    % plot(angle_axis, angle_spectrum, 'LineWidth', 2);
    % xlabel('Angle (degrees)');
    % ylabel('Power');
    % title('AoA Estimation using Delay-and-Sum Beamforming');
    % grid on;
    % %% Range angle Heatmap
    % P = squeeze(mean(angle_map, 2));   % [range × angle]
    % % Polar → Cartesian
    % [AA, RR] = meshgrid(deg2rad(angle_axis), range_axis);
    % X = RR .* cos(AA);
    % Y = RR .* sin(AA);
    % % Plot
    % figure;
    % pcolor(Y, X, 10*log10(P));
    % % pcolor(Y,X,P);
    % shading interp;
    % colormap jet;
    % colorbar;
    % xlabel('X (m)');
    % ylabel('Y (m)');
    % title('Cartesian Range–Angle Heatmap (All Dopplers)');
    % axis equal tight;
end
