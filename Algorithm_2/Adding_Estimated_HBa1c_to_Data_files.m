clc;
clear;
close all;
addpath('../Data/');

%% Load the archival dataset (historical patient data)
load('archival_data_model.mat');  % Assume this contains a structure with various fields

%% Parameters for A1c estimation
a1c_pop_params = [4.15, 0.35];

%% Assuming 'archival_data_model' is a structure with multiple fields
fields = fieldnames(archival_data_model);  % Get all the field names of the loaded data

% Loop over each field in the dataset
for i = 1:length(fields)
    % Extract the input data from each field (assuming it is numerical and valid for the function)
    input = archival_data_model.(fields{i});
    
    % Perform the A1c estimation using the provided function
    eA1c = A1c_Estimation(input, a1c_pop_params);
    
    % Add the eA1c to the corresponding field of the structure
    archival_data_model.(fields{i}).eA1c = eA1c;  % Assuming each field is a struct and you're adding eA1c
    
end

%% Save the modified dataset back to a .mat file
save('archival_data_model.mat', 'archival_data_model');

% Function for A1c Estimation
function eA1c = A1c_Estimation(input, parameters)
    N = length(input);
    
    eA1c = nan(1, N);  % Preallocate the output
    eA1c(1) = parameters(1) + parameters(2) * input(1);  % Initial value
    
    feA1c = parameters(1) + parameters(2) * input;  % Calculate feA1c for all values
    
    % Recursive calculation of eA1c
    for i = 2:N
        eA1c(i) = 0.9512 * eA1c(i-1) + 0.0488 * feA1c(i);
    end
end
