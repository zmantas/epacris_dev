
% Main Script to run Mie related calculations

% Clear the environment
clear all;
close all;

% Define parameters
% Set wavelengths and particle sizes, adjust as necessary
lambda = linspace(0.1, 1, 100); % Wavelength range from 0.1 to 1 micron
particle_sizes = logspace(-1, 1, 100); % Particle sizes from 0.1 to 10 microns

% Complex refractive index m = m' + i*m"
m = 1.5 + 1i * 0.01; % Example values, adjust as necessary

% Compute Mie Coefficients and Efficiencies for each wavelength and size
results = cell(length(lambda), length(particle_sizes));
for i = 1:length(lambda)
    for j = 1:length(particle_sizes)
        x = 2 * pi * particle_sizes(j) / lambda(i); % Size parameter
        coefficients = Mie_abcd(m, x);
        results{i, j} = struct('coefficients', coefficients, 'efficiencies', efficiencies);
    end
end

% Additional calculations for water particles
Mie_cal_water_r1;
