clc;
clear;
close all;
addpath('../Data/');

%% Load the necessary data
load('archival_data_model.mat');  % Load the archival dataset (historical patient data)
load('single_patient_visit_data.mat');  % Load the DCCT single patient's visit data

%% Extract relevant fields from the archival data
archival_field_names = fieldnames(Archival_Data);
archival_log_mu = [Archival_Data.(archival_field_names{1}).pdf.mu]';  % Log-normal distribution mean (mu)
archival_log_sigma = [Archival_Data.(archival_field_names{1}).pdf.sigma]';  % Log-normal distribution standard deviation (sigma)

%% Extract DCCT patient data for a single visit
DCCT_SMBG = single_patient_visit_data.SMBG;  % SMBG values (Self-Monitoring Blood Glucose)
DCCT_HbA1c = single_patient_visit_data.HbA1c;  % HbA1c value (hemoglobin A1c)

%% Define similarity measure for matching (based on log-normal distributions)
Beta = 0.1;  % Margin for similarity calculation
similarity_measure = sum(log(logncdf((1 + Beta) * DCCT_SMBG', archival_log_mu, archival_log_sigma) - ...
    logncdf((1 - Beta) * DCCT_SMBG', archival_log_mu, archival_log_sigma)));

%% Calculate HbA1c distance from archival data
hba1c_distance = (Archival_Data.overal_Hba1cs - DCCT_HbA1c).^2;

%% Initialize the best index search
for n = 0.3:0.1:3
    % Sort the HbA1c distances and select the closest ones based on a threshold
    [sorted_hba1c_distance, original_indices] = sort(abs(hba1c_distance), 'ascend');
    selected_indices = original_indices(sorted_hba1c_distance < n^2);

    if ~isempty(selected_indices)
        % Find the best match by maximizing the similarity measure
        [~, best_local_index] = max(similarity_measure(selected_indices));
        best_index = selected_indices(best_local_index);
        break;
    end
end

%% Generate synthetic SMBG and time data for the matched index
best_mu = archival_log_mu(:, best_index);  % Log-normal mean (mu) for the best match
best_sigma = archival_log_sigma(:, best_index);  % Log-normal standard deviation (sigma) for the best match

days = 90;  % Generate data for 90 days
samples_per_day = 7;  % 7 SMBG samples per day
total_samples = days * samples_per_day;  % Total number of samples = 90 * 7 = 630

% Generate synthetic SMBG values using the log-normal distribution
synthesized_SMBG = lognrnd(repmat(best_mu', total_samples, 1), repmat(best_sigma', total_samples, 1));

% Bound the SMBG values between 30 and 600 mg/dL
synthesized_SMBG(synthesized_SMBG > 600) = 600;
synthesized_SMBG(synthesized_SMBG < 30) = 30;

%% Generate corresponding times for each sample (7 samples per day for 90 days)
time_mu = Archival_Data.(archival_field_names{best_index}).time.mu;  % Mean time (in hours) from archival data
time_sigma = Archival_Data.(archival_field_names{best_index}).time.sigma;  % Time standard deviation

synthesized_times = zeros(total_samples, 1);  % Preallocate the time array
for day = 1:days
    % Generate 7 sample times for each day based on normal distribution
    daily_times = randn(samples_per_day, 1) .* time_sigma + time_mu;
    daily_times(daily_times < 0) = 0;  % Ensure no negative times
    daily_times(daily_times > 24) = daily_times(daily_times > 24) - 24;  % Wrap times around 24 hours

    % Assign the daily times to the total times array
    synthesized_times((day-1)*samples_per_day + 1 : day*samples_per_day) = daily_times + (day-1) * 24;
end

%% Check for any violation in SMBG rate of change (limit: 5 mg/dL per minute)
if any(abs(diff(synthesized_SMBG)) ./ (diff(synthesized_times * 60)) > 5)
    error('Violation detected: SMBG rate of change exceeds 5 mg/dL per minute.');
end

%% Save the generated SMBG and time data to a table
generated_data = table(synthesized_times, synthesized_SMBG, ...
    'VariableNames', {'Time (hours)', 'SMBG (mg/dL)'});

% Save the generated data to a CSV file
writetable(generated_data, 'generated_SMBG_time_data.csv');
