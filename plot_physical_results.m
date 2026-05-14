% MATLAB script to compare Constraint Length K=3 vs K=7
% Load both datasets
phys_data = readtable('physical_results.csv');
sim_data = readtable('results.csv', 'VariableNamingRule', 'preserve');

%% 1. Process Physical Results (K=7)
% Exclude the known synchronization outlier at 1.0m
outlierIdx = (phys_data.distance == 1.0 & phys_data.trial == 2 & ...
              strcmp(phys_data.modscheme, 'qpsk') & phys_data.useCoding == 1);
phys_data(outlierIdx, :) = []; 

% Filter for QPSK Distance Sweep
phys_qpsk = phys_data(strcmp(phys_data.test, 'distance_sweep') & strcmp(phys_data.modscheme, 'qpsk'), :);

% Group by distance for K=7 Coded and Uncoded
[G_p7, d_p7] = findgroups(phys_qpsk.distance(phys_qpsk.useCoding == 1));
ber_p7 = splitapply(@mean, phys_qpsk.ber(phys_qpsk.useCoding == 1), G_p7);

[G_pu, d_pu] = findgroups(phys_qpsk.distance(phys_qpsk.useCoding == 0));
ber_pu = splitapply(@mean, phys_qpsk.ber(phys_qpsk.useCoding == 0), G_pu);

%% 2. Process Results.csv (K=3)
% Filter for QPSK Distance Sweep and remove sync failure outliers
k3_qpsk = sim_data(strcmp(sim_data.test, 'distance_sweep_qpsk') & sim_data.M == 4, :);
k3_qpsk = k3_qpsk(k3_qpsk.ber < 0.1, :); 

% Group by distance for K=3 Coded
[G_k3, d_k3] = findgroups(k3_qpsk.distance(k3_qpsk.useCoding == 1));
ber_k3 = splitapply(@mean, k3_qpsk.ber(k3_qpsk.useCoding == 1), G_k3);

%% 3. Plotting the Comparison
figure('Color', 'w', 'Name', 'Coding Comparison');
hold on; grid on;

% Plot Uncoded Baseline
plot(d_pu, max(ber_pu, 1e-6), '--k', 'LineWidth', 1.2, 'DisplayName', 'Uncoded QPSK');

% Plot K=3 Coded (from results.csv)
plot(d_k3, max(ber_k3, 1e-6), '-s', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Coded QPSK (K=3)');

% Plot K=7 Coded (from physical_results.csv)
plot(d_p7, max(ber_p7, 1e-6), '-o', 'Color', [0 0.45 0.74], 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Coded QPSK (K=7)');

% Formatting
set(gca, 'YScale', 'log');
ylim([1e-7 1]);
yticks([1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1]);
yticklabels({'0 (Error Free)', '10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '1'});

xlabel('Distance (m)');
ylabel('Bit Error Rate (BER)');
title('Impact of Constraint Length on System Robustness');
legend('Location', 'best');

% Export
exportgraphics(gcf, './plotpdfs/k_comparison_ber.pdf', 'ContentType', 'vector');