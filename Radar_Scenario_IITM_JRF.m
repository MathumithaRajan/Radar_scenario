%% Radar Specifications 
% Frequency of operation = 20kHz
% Max Range = 500m 
% Range Resolution = 1 m
% Max Velocity = 150 m/s 

max_range = 500;
range_resolution = 1;
max_velocity = 150;
c = 3 * 10^8;               % Speed of light

%% target Parameters
% Initial range and velocity of the targets
target1_range = 300;        % meters (stationary target)
target1_velocity = 0;       % m/s

target2_range = 500;        % meters (moving target)
target2_velocity = 100;     % m/s


%% FMCW Waveform Generation

% Calculate Bandwidth (B), Chirp Time (Tchirp), and Slope (slope)
B = c / (2 * range_resolution);
Tchirp = 5.5 * 2 * max_range / c;
slope = B / Tchirp;

% Operating carrier frequency of Radar 
fc = 20 * 10^3;             % carrier frequency

% The number of chirps in one sequence (number of Doppler cells)
Nd = 128;  % Number of chirps

% The number of samples on each chirp (number of range cells)
Nr = 1024;  % Number of samples

% Timestamp for running the displacement scenario for every sample on each chirp
t = linspace(0, Nd * Tchirp, Nr * Nd);  % total time for samples

% Creating the vectors for Tx, Rx, and Mix based on the total samples input.
Tx = zeros(1, length(t)); % transmitted signal
Rx = zeros(1, length(t)); % received signal
Mix = zeros(1, length(t)); % beat signal

% Similar vectors for range_covered and time delay.
r_t = zeros(1, length(t));
td1 = zeros(1, length(t));
td2 = zeros(1, length(t));

%% Signal generation and Moving Target simulation

% Running the radar scenario over the time.
for i = 1:length(t)
    % Update the Range of the Targets for constant velocity.
    r_t1 = target1_range + target1_velocity * t(i);
    r_t2 = target2_range + target2_velocity * t(i);
    
    td1(i) = 2 * r_t1 / c;
    td2(i) = 2 * r_t2 / c;
    
    % Update the transmitted and received signal.
    Tx(i) = cos(2 * pi * (fc * t(i) + (slope * t(i)^2) / 2));
    Rx1 = cos(2 * pi * (fc * (t(i) - td1(i)) + (slope * (t(i) - td1(i))^2) / 2));
    Rx2 = cos(2 * pi * (fc * (t(i) - td2(i)) + (slope * (t(i) - td2(i))^2) / 2));
    
    % Combine the received signals from both targets
    Rx(i) = Rx1 + Rx2;
    
    % Mix the Transmit and Receive signals to generate the beat signal
    Mix(i) = Tx(i) .* Rx(i);
end

%% RANGE MEASUREMENT
% Reshape the vector into Nr*Nd array. Nr and Nd here define the size of Range and Doppler FFT respectively.
sig = reshape(Mix, [Nr, Nd]);

% Run the FFT on the beat signal along the range bins dimension (Nr)
sig_fft1 = fft(sig, Nr);

% Normalize
sig_fft1 = sig_fft1 / Nr;

% Take the absolute value of FFT output
sig_fft1 = abs(sig_fft1);

% Output of FFT is double-sided signal, but we are interested in only one side of the spectrum.
sig_fft1 = sig_fft1(1:Nr/2, :);

% Plotting the range
figure('Name', 'Range from First FFT');
plot(linspace(0, max_range, Nr/2), max(sig_fft1, [], 2), 'LineWidth', 2);
grid on;
axis([0 500 0 max(max(sig_fft1))]);
xlabel('Range (m)');
ylabel('Amplitude');
title('RANGE FFT');

%% Estimation of Range
% Find the peak value in the range FFT output
[max_val, range_idx] = max(max(sig_fft1, [], 2));
estimated_range = range_idx * (max_range / (Nr/2));

% Display the estimated range
fprintf('Estimated Range: %.2f m\n', estimated_range);

%% DOPPLER MEASUREMENT
% Run the FFT on the beat signal along the Doppler bins dimension (Nd)
sig_fft2 = fftshift(fft(sig_fft1, Nd, 2));

% Normalize
sig_fft2 = sig_fft2 / Nd;

% Take the absolute value of FFT output
sig_fft2 = abs(sig_fft2);

% Plotting the Doppler FFT output
doppler_axis = linspace(-max_velocity, max_velocity, Nd);
figure('Name', 'Doppler FFT');
plot(doppler_axis, max(sig_fft2, [], 1), 'LineWidth', 2);
grid on;
xlabel('Velocity (m/s)');
ylabel('Amplitude');
title('Doppler FFT');

%% Estimation of Velocity
% Find the peak value in the Doppler FFT output
[~, velocity_idx] = max(max(sig_fft2, [], 1));
estimated_velocity = (velocity_idx - Nd/2) * (2 * max_velocity / Nd);

% Display the estimated velocity
fprintf('Estimated Velocity: %.2f m/s\n', estimated_velocity);
