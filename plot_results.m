%% plot_results.m
clear all; clc;

data = readtable('results.csv', 'VariableNamingRule', 'preserve');

set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultTextFontSize', 13);
set(0, 'DefaultAxesGridAlpha', 0.4);
set(0, 'DefaultAxesMinorGridAlpha', 0.3);

%% ---- 1. CP SWEEP ----
figure('Name', 'CP Length Sweep');

coded   = data(strcmp(data.test, 'cp_sweep_coded'), :);
uncoded = data(strcmp(data.test, 'cp_sweep_uncoded'), :);

cp_vals    = unique(uncoded.cplen);
ber_mean_u = zeros(size(cp_vals));
ber_std_u  = zeros(size(cp_vals));
ber_mean_c = zeros(size(cp_vals));

for i = 1:length(cp_vals)
    mask_u = uncoded.cplen == cp_vals(i);
    ber_mean_u(i) = mean(uncoded.ber(mask_u));
    ber_std_u(i)  = std(uncoded.ber(mask_u));
    mask_c = coded.cplen == cp_vals(i);
    if any(mask_c)
        ber_mean_c(i) = mean(coded.ber(mask_c));
    end
end

errorbar(cp_vals, ber_mean_u, ber_std_u, '-o', 'LineWidth', 1.5, ...
    'DisplayName', 'Uncoded'); hold on;
plot(cp_vals, ber_mean_c, '-o', 'LineWidth', 1.5, ...
    'DisplayName', 'Coded (r=1/2)');
xline(256, '--', 'LineWidth', 1.2, 'Label', 'Selected: 256', ...
    'FontSize', 9, 'LabelHorizontalAlignment', 'right');

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('CP Length (samples)');
ylabel('BER');
title('BER vs CP Length');
legend('Location', 'best');
grid on;

exportgraphics(gcf, './plotpdfs/cp_sweep.pdf', 'ContentType', 'vector');

%% ---- 2. PILOT SPACING SWEEP ----
figure('Name', 'Pilot Spacing Sweep');

pilots  = data(strcmp(data.test, 'pilot_sweep') & data.numActiveCarriers == 300, :);
sp_vals = unique(pilots.pilotSpacing);
ber_mean_p = zeros(size(sp_vals));
ber_std_p  = zeros(size(sp_vals));
tput_mean_p = zeros(size(sp_vals));

for i = 1:length(sp_vals)
    mask = pilots.pilotSpacing == sp_vals(i);
    ber_mean_p(i)  = mean(pilots.ber(mask));
    ber_std_p(i)   = std(pilots.ber(mask));
    tput_mean_p(i) = mean(pilots.throughput(mask));
end

yyaxis left
errorbar(sp_vals, ber_mean_p, ber_std_p, '-o', 'LineWidth', 1.5);
ylabel('BER');

yyaxis right
plot(sp_vals, tput_mean_p/1000, '-s', 'LineWidth', 1.5);
ylabel('Throughput (kbps)');

ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
xline(5, '--', 'LineWidth', 1.2, 'Label', 'Selected: 5', 'FontSize', 9);
xlabel('Pilot Spacing');
title('BER and Throughput vs Pilot Spacing');
legend({'BER', 'Throughput (kbps)'}, 'Location', 'best');
grid on;

exportgraphics(gcf, './plotpdfs/pilot_sweep.pdf', 'ContentType', 'vector');

%% ---- 3. CARRIER COUNT SWEEP ----
figure('Name', 'Carrier Count Sweep');

carriers = data(strcmp(data.test, 'carrier_sweep'), :);
c_vals   = unique(carriers.numActiveCarriers);
ber_mean_carr  = zeros(size(c_vals));
ber_std_carr   = zeros(size(c_vals));
tput_mean_carr = zeros(size(c_vals));

for i = 1:length(c_vals)
    mask = carriers.numActiveCarriers == c_vals(i);
    ber_mean_carr(i)  = mean(carriers.ber(mask));
    ber_std_carr(i)   = std(carriers.ber(mask));
    tput_mean_carr(i) = mean(carriers.throughput(mask));
end

yyaxis left
errorbar(c_vals, ber_mean_carr, ber_std_carr, '-o', 'LineWidth', 1.5);
ylabel('BER');

yyaxis right
plot(c_vals, tput_mean_carr/1000, '-s', 'LineWidth', 1.5);
ylabel('Throughput (kbps)');

ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
xline(500, '--', 'LineWidth', 1.2, 'Label', 'Selected: 500', 'FontSize', 9);
xlabel('Number of Active Carriers');
title('BER and Throughput vs Carrier Count');
legend({'BER', 'Throughput (kbps)'}, 'Location', 'northwest');
grid on;

exportgraphics(gcf, './plotpdfs/carrier_sweep.pdf', 'ContentType', 'vector');

%% ---- 4. FC + CARRIER HEATMAP ----
figure('Name', 'FC vs Carrier Count BER Heatmap');

fc_data = data(strcmp(data.test, 'fc_carrier_sweep') & ...
               data.M == 4 & data.useCoding == 0, :);

fc_vals = unique(fc_data.fc);
c_vals2 = unique(fc_data.numActiveCarriers);
ber_matrix = nan(length(fc_vals), length(c_vals2));

for i = 1:length(fc_vals)
    for j = 1:length(c_vals2)
        mask = fc_data.fc == fc_vals(i) & fc_data.numActiveCarriers == c_vals2(j);
        if any(mask)
            ber_matrix(i,j) = mean(fc_data.ber(mask));
        end
    end
end

imagesc(c_vals2/1000*23.4, fc_vals/1000, ber_matrix);
colormap(flipud(hot));
cb = colorbar;
cb.Label.String = 'BER';

xlabel('Approx. Bandwidth (kHz)');
ylabel('Centre Frequency (kHz)');
title('BER Heatmap: Centre Frequency vs Carrier Count');
set(gca, 'XTick', c_vals2/1000*23.4, ...
    'XTickLabel', arrayfun(@num2str, c_vals2, 'UniformOutput', false));
set(gca, 'YTick', fc_vals/1000, ...
    'YTickLabel', arrayfun(@(x) sprintf('%gk', x), fc_vals/1000, 'UniformOutput', false));

hold on;
plot(500/1000*23.4, 10, 'p', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'w', ...
    'DisplayName', 'Selected Config');
legend('Location', 'best');

exportgraphics(gcf, './plotpdfs/fc_heatmap.pdf', 'ContentType', 'vector');

%% ---- 5. DISTANCE SWEEP ----
figure('Name', 'Distance Sweep');

dist = data(strcmp(data.test, 'distance_sweep_qpsk') & data.M == 4, :);
dist = dist(dist.ber < 0.1, :);  % remove sync failure outliers

d_vals    = unique(dist.distance);
coded_d   = dist(dist.useCoding == 1, :);
uncoded_d = dist(dist.useCoding == 0, :);

bm_c = zeros(size(d_vals)); bs_c = zeros(size(d_vals));
bm_u = zeros(size(d_vals)); bs_u = zeros(size(d_vals));

for i = 1:length(d_vals)
    mc = coded_d.distance == d_vals(i);
    mu = uncoded_d.distance == d_vals(i);
    if any(mc)
        bm_c(i) = mean(coded_d.ber(mc));
        bs_c(i) = std(coded_d.ber(mc));
    end
    if any(mu)
        bm_u(i) = mean(uncoded_d.ber(mu));
        bs_u(i) = std(uncoded_d.ber(mu));
    end
end

% clamp zeros so they show on log scale
bm_c = max(bm_c, 1e-6);
bm_u = max(bm_u, 1e-6);

errorbar(d_vals, bm_u, bs_u, '-o', 'LineWidth', 1.5, 'DisplayName', 'Uncoded'); hold on;
errorbar(d_vals, bm_c, bs_c, '-o', 'LineWidth', 1.5, 'DisplayName', 'Coded (r=1/2)');

xlabel('Distance (m)');
ylabel('BER');
title('BER vs Distance — QPSK');
legend('Location', 'best');
set(gca, 'YScale', 'log');
grid on;

exportgraphics(gcf, './plotpdfs/distance_ber.pdf', 'ContentType', 'vector');

%% ---- 6. SNR vs DISTANCE ----
figure('Name', 'SNR vs Distance');

dist_clean = data(strcmp(data.test, 'distance_sweep_qpsk') & ...
                  data.M == 4 & data.ber < 0.1, :);

coded_d   = dist_clean(dist_clean.useCoding == 1, :);
uncoded_d = dist_clean(dist_clean.useCoding == 0, :);

d_vals = unique(dist_clean.distance);
snr_c = zeros(size(d_vals)); snr_u = zeros(size(d_vals));
snr_c_std = zeros(size(d_vals)); snr_u_std = zeros(size(d_vals));

for i = 1:length(d_vals)
    mc = coded_d.distance == d_vals(i);
    mu = uncoded_d.distance == d_vals(i);
    if any(mc)
        snr_c(i)     = mean(coded_d.snr(mc));
        snr_c_std(i) = std(coded_d.snr(mc));
    end
    if any(mu)
        snr_u(i)     = mean(uncoded_d.snr(mu));
        snr_u_std(i) = std(uncoded_d.snr(mu));
    end
end

errorbar(d_vals, snr_u, snr_u_std, '-o', 'LineWidth', 1.5, 'DisplayName', 'Uncoded'); hold on;
errorbar(d_vals, snr_c, snr_c_std, '-o', 'LineWidth', 1.5, 'DisplayName', 'Coded');

xlabel('Distance (m)');
ylabel('Estimated SNR (dB)');
title('SNR vs Distance');
legend('Location', 'best');
grid on;

exportgraphics(gcf, './plotpdfs/distance_snr.pdf', 'ContentType', 'vector');

%% ---- 7. BER vs DISTANCE with SNR overlay ----
figure('Name', 'BER and SNR vs Distance');

yyaxis left
errorbar(d_vals, bm_u, bs_u, '-o', 'LineWidth', 1.5, 'DisplayName', 'BER Uncoded'); hold on;
errorbar(d_vals, bm_c, bs_c, '-o', 'LineWidth', 1.5, 'DisplayName', 'BER Coded');
ylabel('BER');

yyaxis right
plot(d_vals, snr_u, '--s', 'LineWidth', 1.2, 'DisplayName', 'SNR Uncoded');
plot(d_vals, snr_c, '--s', 'LineWidth', 1.2, 'DisplayName', 'SNR Coded');
ylabel('SNR (dB)');

ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
xlabel('Distance (m)');
title('BER and SNR vs Distance');
legend('Location', 'northwest');
grid on;

exportgraphics(gcf, './plotpdfs/ber_snr_distance.pdf', 'ContentType', 'vector');

%% ---- 8. CODING GAIN vs DISTANCE ----
figure('Name', 'Coding Gain vs Distance');

coding_gain     = zeros(size(d_vals));
coding_gain_snr = zeros(size(d_vals));

for i = 1:length(d_vals)
    mc = coded_d.distance == d_vals(i);
    mu = uncoded_d.distance == d_vals(i);
    if any(mc) && any(mu)
        coding_gain(i)     = mean(uncoded_d.ber(mu)) - mean(coded_d.ber(mc));
        coding_gain_snr(i) = mean(coded_d.snr(mc))   - mean(uncoded_d.snr(mu));
    end
end

yyaxis left
bar(d_vals, coding_gain, 0.4, 'EdgeColor', 'none', 'DisplayName', 'BER Reduction');
ylabel('\Delta BER (Uncoded - Coded)');

yyaxis right
plot(d_vals, coding_gain_snr, '-o', 'LineWidth', 1.5, 'DisplayName', 'SNR Gain (dB)');
ylabel('SNR Gain (dB)');

ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
xlabel('Distance (m)');
title('Coding Gain vs Distance');
legend('Location', 'northwest');
grid on;

exportgraphics(gcf, './plotpdfs/coding_gain.pdf', 'ContentType', 'vector');

%% ---- 9. BER vs CENTRE FREQUENCY ----
figure('Name', 'BER vs Centre Frequency');

fc_data = data(strcmp(data.test, 'fc_carrier_sweep') & ...
               data.M == 4 & data.useCoding == 0, :);

carrier_counts = [100, 200, 300, 400, 500];
fc_vals = unique(fc_data.fc);

for ci = 1:length(carrier_counts)
    nc = carrier_counts(ci);
    ber_fc     = zeros(size(fc_vals));
    ber_fc_std = zeros(size(fc_vals));

    for fi = 1:length(fc_vals)
        mask = fc_data.numActiveCarriers == nc & fc_data.fc == fc_vals(fi);
        if any(mask)
            ber_fc(fi)     = mean(fc_data.ber(mask));
            ber_fc_std(fi) = std(fc_data.ber(mask));
        else
            ber_fc(fi) = NaN;
        end
    end

    valid = ~isnan(ber_fc);
    errorbar(fc_vals(valid)/1000, ber_fc(valid), ber_fc_std(valid), '-o', ...
        'LineWidth', 1.5, 'DisplayName', sprintf('%d carriers', nc));
    hold on;
end

x = xline(10, '--', 'LineWidth', 1.2, 'Label', 'Selected: 10kHz', 'FontSize', 9);
x.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlabel('Centre Frequency (kHz)');
ylabel('BER');
title('BER vs Centre Frequency by Carrier Count');
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
grid on;

exportgraphics(gcf, './plotpdfs/ber_vs_fc.pdf', 'ContentType', 'vector');