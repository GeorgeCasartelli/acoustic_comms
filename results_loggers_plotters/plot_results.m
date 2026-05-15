%% plot_results.m
clear all; clc;

data = readtable('results.csv', 'VariableNamingRule', 'preserve');

% Shared style
bgCol   = [0.12 0.12 0.16];
axCol   = [0.18 0.18 0.24];
txtCol  = [0.95 0.95 0.95];
gridCol = [0.30 0.30 0.38];
accent1 = [0.26 0.63 0.95];
accent2 = [0.95 0.45 0.15];
accent3 = [0.25 0.88 0.58];
accent4 = [0.95 0.30 0.45];

styleAx = @(ax) set(ax, ...
    'Color', axCol, 'XColor', txtCol, 'YColor', txtCol, ...
    'GridColor', gridCol, 'GridAlpha', 0.5, 'FontSize', 11, ...
    'TitleFontSizeMultiplier', 1.1);

%% ---- 1. CP SWEEP ----
figure('Color', bgCol, 'Name', 'CP Length vs BER');
ax = axes; styleAx(ax);

coded   = data(strcmp(data.test, 'cp_sweep_coded'), :);
uncoded = data(strcmp(data.test, 'cp_sweep_uncoded'), :);

cp_vals    = unique(uncoded.cplen);
ber_mean_u = zeros(size(cp_vals));
ber_std_u  = zeros(size(cp_vals));
ber_mean_c = zeros(size(cp_vals));
ber_std_c  = zeros(size(cp_vals));

for i = 1:length(cp_vals)
    mu = uncoded.cplen == cp_vals(i);
    mc = coded.cplen == cp_vals(i);
    ber_mean_u(i) = mean(uncoded.ber(mu));
    ber_std_u(i)  = std(uncoded.ber(mu));
    if any(mc)
        ber_mean_c(i) = mean(coded.ber(mc));
        ber_std_c(i)  = std(coded.ber(mc));
    end
end

errorbar(cp_vals, ber_mean_u, ber_std_u, '-o', ...
    'Color', accent2, 'MarkerFaceColor', accent2, ...
    'LineWidth', 1.5, 'DisplayName', 'Uncoded'); hold on;
errorbar(cp_vals, ber_mean_c, ber_std_c, '-o', ...
    'Color', accent3, 'MarkerFaceColor', accent3, ...
    'LineWidth', 1.5, 'DisplayName', 'Coded (r=1/2)');

xline(256, '--', 'Color', accent1, 'LineWidth', 1.2, ...
    'Label', 'Selected: 256', 'FontSize', 9, ...
    'LabelHorizontalAlignment', 'right');

set(gca, 'XScale', 'log');
xlabel('CP Length (samples)', 'Color', txtCol);
ylabel('BER', 'Color', txtCol);
title('BER vs CP Length', 'Color', txtCol);
legend('TextColor', txtCol, 'Color', axCol, 'EdgeColor', gridCol);
grid on;

%% ---- 2. PILOT SPACING ----
figure('Color', bgCol, 'Name', 'Pilot Spacing vs BER');
ax = axes; styleAx(ax);

pilots = data(strcmp(data.test, 'pilot_sweep') & ...
              data.numActiveCarriers == 300, :);

sp_vals    = unique(pilots.pilotSpacing);
ber_mean_p = zeros(size(sp_vals));
ber_std_p  = zeros(size(sp_vals));

for i = 1:length(sp_vals)
    mask = pilots.pilotSpacing == sp_vals(i);
    ber_mean_p(i) = mean(pilots.ber(mask));
    ber_std_p(i)  = std(pilots.ber(mask));
end

errorbar(sp_vals, ber_mean_p, ber_std_p, '-o', ...
    'Color', accent2, 'MarkerFaceColor', accent2, 'LineWidth', 1.5);

xline(5, '--', 'Color', accent1, 'LineWidth', 1.2, ...
    'Label', 'Selected: 5', 'FontSize', 9);

xlabel('Pilot Spacing', 'Color', txtCol);
ylabel('BER', 'Color', txtCol);
title('BER vs Pilot Spacing', 'Color', txtCol);
grid on;

%% ---- 3. CARRIER COUNT ----
figure('Color', bgCol, 'Name', 'Carrier Count vs BER');
ax = axes; styleAx(ax);

carriers = data(strcmp(data.test, 'carrier_sweep'), :);

c_vals         = unique(carriers.numActiveCarriers);
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
errorbar(c_vals, ber_mean_carr, ber_std_carr, '-o', ...
    'Color', accent2, 'MarkerFaceColor', accent2, 'LineWidth', 1.5);
ylabel('BER', 'Color', txtCol);
set(gca, 'YColor', txtCol);

yyaxis right
plot(c_vals, tput_mean_carr/1000, '-s', ...
    'Color', accent3, 'MarkerFaceColor', accent3, 'LineWidth', 1.5);
ylabel('Throughput (kbps)', 'Color', txtCol);
set(gca, 'YColor', txtCol);

xline(500, '--', 'Color', accent1, 'LineWidth', 1.2, ...
    'Label', 'Selected: 500', 'FontSize', 9);

xlabel('Number of Active Carriers', 'Color', txtCol);
title('BER and Throughput vs Carrier Count', 'Color', txtCol);
grid on;

%% ---- 4. BER vs FC LINE PLOT ----
figure('Color', bgCol, 'Name', 'BER vs Centre Frequency');
ax = axes; styleAx(ax);

fc_data = data(strcmp(data.test, 'fc_carrier_sweep') & ...
               data.M == 4 & data.useCoding == 0, :);

carrier_counts = [100, 200, 300, 400, 500];
colours = {accent1, accent2, accent3, accent4, [0.8 0.6 0.9]};
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
        'Color', colours{ci}, 'MarkerFaceColor', colours{ci}, ...
        'LineWidth', 1.5, 'DisplayName', sprintf('%d carriers', nc));
    hold on;
end

xline(10, '--', 'Color', txtCol, 'Alpha', 0.5, 'LineWidth', 1.2, ...
    'Label', 'Selected: 10kHz', 'FontSize', 9);

xlabel('Centre Frequency (kHz)', 'Color', txtCol);
ylabel('BER', 'Color', txtCol);
title('BER vs Centre Frequency by Carrier Count', 'Color', txtCol);
legend('TextColor', txtCol, 'Color', axCol, 'EdgeColor', gridCol, ...
    'Location', 'northwest');
grid on;

%% ---- 5. BER vs DISTANCE ----
figure('Color', bgCol, 'Name', 'BER vs Distance');
ax = axes; styleAx(ax);

dist_clean = data(strcmp(data.test, 'distance_sweep_qpsk') & ...
                  data.M == 4 & data.ber < 0.1, :);

coded_d   = dist_clean(dist_clean.useCoding == 1, :);
uncoded_d = dist_clean(dist_clean.useCoding == 0, :);

d_vals = unique(dist_clean.distance);
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

errorbar(d_vals, bm_u, bs_u, '-o', ...
    'Color', accent2, 'MarkerFaceColor', accent2, ...
    'LineWidth', 1.5, 'DisplayName', 'Uncoded'); hold on;
errorbar(d_vals, bm_c, bs_c, '-o', ...
    'Color', accent3, 'MarkerFaceColor', accent3, ...
    'LineWidth', 1.5, 'DisplayName', 'Coded (r=1/2)');

xlabel('Distance (m)', 'Color', txtCol);
ylabel('BER', 'Color', txtCol);
title('BER vs Distance — QPSK', 'Color', txtCol);
legend('TextColor', txtCol, 'Color', axCol, 'EdgeColor', gridCol);
grid on;

%% ---- 6. SNR vs DISTANCE ----
figure('Color', bgCol, 'Name', 'SNR vs Distance');
ax = axes; styleAx(ax);

snr_c     = zeros(size(d_vals)); snr_c_std = zeros(size(d_vals));
snr_u     = zeros(size(d_vals)); snr_u_std = zeros(size(d_vals));

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

errorbar(d_vals, snr_u, snr_u_std, '-o', ...
    'Color', accent2, 'MarkerFaceColor', accent2, ...
    'LineWidth', 1.5, 'DisplayName', 'Uncoded'); hold on;
errorbar(d_vals, snr_c, snr_c_std, '-o', ...
    'Color', accent3, 'MarkerFaceColor', accent3, ...
    'LineWidth', 1.5, 'DisplayName', 'Coded');

xlabel('Distance (m)', 'Color', txtCol);
ylabel('Estimated SNR (dB)', 'Color', txtCol);
title('SNR vs Distance', 'Color', txtCol);
legend('TextColor', txtCol, 'Color', axCol, 'EdgeColor', gridCol);
grid on;

%% ---- 7. CODING GAIN vs DISTANCE ----
figure('Color', bgCol, 'Name', 'Coding Gain vs Distance');
ax = axes; styleAx(ax);

coding_gain = zeros(size(d_vals));

for i = 1:length(d_vals)
    mc = coded_d.distance == d_vals(i);
    mu = uncoded_d.distance == d_vals(i);
    if any(mc) && any(mu)
        coding_gain(i) = mean(uncoded_d.ber(mu)) - mean(coded_d.ber(mc));
    end
end

bar(d_vals, coding_gain, 0.4, 'FaceColor', accent1, 'EdgeColor', 'none');

xlabel('Distance (m)', 'Color', txtCol);
ylabel('\Delta BER (Uncoded - Coded)', 'Color', txtCol);
title('Coding Gain vs Distance', 'Color', txtCol);
grid on;

%% ---- 8. MODULATION COMPARISON ----
figure('Color', bgCol, 'Name', 'Modulation Comparison');
ax = axes; styleAx(ax); axis off;

qpsk_final = data(strcmp(data.test, 'fc_carrier_sweep') & ...
                  data.M == 4 & data.useCoding == 1 & ...
                  data.numActiveCarriers == 500 & data.fc == 10000, :);

qam_close = data(strcmp(data.test, 'fc_carrier_sweep') & ...
                 data.M == 16 & data.useCoding == 1 & ...
                 data.distance <= 0.5, :);

qam_15m   = data(strcmp(data.test, 'fc_carrier_sweep') & ...
                 data.M == 16 & data.useCoding == 1 & ...
                 data.distance == 1.5, :);

rows = {
    'QPSK 500c, 10kHz, 1.5m, coded', ...
        sprintf('%.0f bps', mean(qpsk_final.throughput)), ...
        sprintf('%.5f', mean(qpsk_final.ber));
    '16-QAM 500c, 10kHz, 0.5m, coded', ...
        sprintf('%.0f bps', mean(qam_close.throughput)), ...
        sprintf('%.5f', mean(qam_close.ber));
    '16-QAM 500c, 10kHz, 1.5m, coded', ...
        sprintf('%.0f bps', mean(qam_15m.throughput)), ...
        sprintf('%.5f', mean(qam_15m.ber));
};

nRows  = size(rows, 1);
yStart = 0.78;
yStep  = 0.18;

text(0.05, 0.93, 'Modulation Scheme Comparison', 'Color', txtCol, ...
    'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized');
text(0.05, 0.85, 'Config', 'Color', [0.65 0.65 0.75], ...
    'FontSize', 11, 'Units', 'normalized');
text(0.60, 0.85, 'Throughput', 'Color', [0.65 0.65 0.75], ...
    'FontSize', 11, 'Units', 'normalized');
text(0.80, 0.85, 'BER', 'Color', [0.65 0.65 0.75], ...
    'FontSize', 11, 'Units', 'normalized');

for i = 1:nRows
    yPos = yStart - (i-1)*yStep;
    text(0.05, yPos, rows{i,1}, 'Color', txtCol, ...
        'FontSize', 10, 'Units', 'normalized');
    text(0.60, yPos, rows{i,2}, 'Color', accent3, ...
        'FontSize', 10, 'FontWeight', 'bold', 'Units', 'normalized');
    text(0.80, yPos, rows{i,3}, 'Color', accent2, ...
        'FontSize', 10, 'FontWeight', 'bold', 'Units', 'normalized');
end;