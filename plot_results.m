data = readtable('results.csv', 'VariableNamingRule', 'preserve');

%% CP SWEEP
figure;
coded   = data(strcmp(data.test, 'cp_sweep_coded'), :);
uncoded = data(strcmp(data.test, 'cp_sweep_uncoded'), :);

cp_vals = unique(uncoded.cplen);
ber_mean = zeros(size(cp_vals));
ber_std  = zeros(size(cp_vals));

for i = 1:length(cp_vals)
    mask = uncoded.cplen == cp_vals(i);
    ber_mean(i) = mean(uncoded.ber(mask));
    ber_std(i)  = std(uncoded.ber(mask));
end

errorbar(cp_vals, ber_mean, ber_std, '-o', 'LineWidth', 1.5, 'DisplayName', 'Uncoded'); hold on;
yline(mean(coded.ber), '--', 'DisplayName', 'Coded (BER \approx 0)');

set(gca, 'XScale', 'log');
xlabel('CP Length (samples)');
ylabel('BER');
title('BER vs CP Length');
legend; grid on;