clc;
clear;
close all;

%% Load Archival Data and Apply Sliding Window
load('archival_data_model.mat');  % Load archival data (historical patient data)

% Parameters
windowSize = 90;
t_intervals = [5 9; 6.5 10.5; 10 14; 11.5 15.5; 16 20; 17.5 21.5; 20 24];  % Time intervals

% Initialize struct to hold sliding window data
smbgValues = struct();

% Apply sliding window to the data for valid entries
for i = 1:length(archival_data_model)
    idxData = archival_data_model(i);
    for j = 1:(max(idxData.day) - windowSize + 1)
        windowData = idxData(idxData.day >= j & idxData.day < j + windowSize, :);
        if length(windowData.day) >= 90
            smbgValues.(sprintf('idx%d_win%d', i, j)) = [windowData.day, windowData.TimeOfDay, windowData.SMBG, windowData.Hba1c_Estimated];
        end
    end
end

%% Filter Data to Retain Entries with Minimum Sample Count
filtered_data = filter_data(smbgValues, t_intervals);  % Filter data to have at least 30 samples

%% Fit Log-Normal Distributions to the Filtered Data
pdf_params = struct();
field_names = fieldnames(filtered_data);

for i = 1:length(field_names)
    filtered_data_i = filtered_data.(field_names{i});
    lognormal_params = fit_lognormal_to_smbg(filtered_data_i, t_intervals);
    pdf_params.(field_names{i}) = lognormal_params;
end

%% Visualize Results with Histograms and Log-Normal Fit
plot_histogram_with_fit(pdf_params, t_intervals);

%% Save the Processed Data
save('archival_data_model.mat', 'pdf_params');

%% Function Definitions

% Filter Data Based on Time Intervals and Minimum Sample Count
function filtered_data = filter_data(mat_data, t_intervals)
    fields = fieldnames(mat_data);
    filtered_data = struct();
    
    for i = 1:numel(fields)
        data = mat_data.(fields{i});
        times = data(:, 2);  % Extract times
        counts = zeros(size(t_intervals, 1), 1);
        
        for j = 1:size(t_intervals, 1)
            counts(j) = sum(times >= t_intervals(j, 1) & times < t_intervals(j, 2));
        end
        
        if all(counts >= 30)
            filtered_data.(fields{i}) = data;
        end
    end
end

% Fit Log-Normal Distributions to SMBG Data
function lognormal_params = fit_lognormal_to_smbg(filtered_data, t_intervals)
    time = filtered_data(:, 2);
    smbg_values = filtered_data(:, 3);
    lognormal_params = struct();
    
    for i = 1:size(t_intervals, 1)
        interval_mask = time >= t_intervals(i, 1) & time < t_intervals(i, 2);
        interval_values = smbg_values(interval_mask);
        
        pd = fitdist(log(interval_values), 'Normal');
        lognormal_params.pdf(i).mu = pd.mu;
        lognormal_params.pdf(i).sigma = pd.sigma;
    end
end

% Plot Histograms with Log-Normal Fit for Each Time Interval
function plot_histogram_with_fit(pdf_params, t_intervals)
    num_intervals = size(t_intervals, 1);
    field_names = fieldnames(pdf_params);
    
    for j = 1:length(field_names)
        first_window_data = pdf_params.(field_names{j});
        figure;
        
        for i = 1:num_intervals
            subplot(ceil(num_intervals/2), 2, i);
            hold on;
            
            time = first_window_data.SMBGs(:, 2);
            smbg_values = first_window_data.SMBGs(:, 3);
            interval_mask = time >= t_intervals(i, 1) & time < t_intervals(i, 2);
            interval_data = smbg_values(interval_mask);
            
            histogram(interval_data, 'Normalization', 'pdf', 'BinWidth', 10);
            
            mu = first_window_data.pdf(i).mu;
            sigma = first_window_data.pdf(i).sigma;
            smbg_range = linspace(30, 600, 1000);
            pdf_values = lognpdf(smbg_range, mu, sigma);
            plot(smbg_range, pdf_values, 'r', 'LineWidth', 2);
            
            title(sprintf('Interval %0.1f to %0.1f', t_intervals(i, 1), t_intervals(i, 2)));
            xlabel('SMBG Values');
            ylabel('Probability Density');
            hold off;
        end
    end
end
