function display_ofdm_dashboard(rx, fs, start, recoveredLen, headerSize, dataIdx, k, ...
    fc, numActiveCarriers, H_interp, rxPilots, pilotSym, pilotIdx, nfft, cplen, ...
    useCoding, useInterleaving, txMode, actualBER, R_b, modScheme, snr_est, ...
    constraintLength, serPerSubcarrier)
    %% -- Initial Setup --
    symbolLen = nfft + cplen;
    bw_theory = (numActiveCarriers * fs) / nfft;

    %% -- Create Dashboard Figure --
    figure('Name', 'OFDM Analysis Dashboard Pro', 'Units', 'normalized', 'OuterPosition', [0.05 0.05 0.9 0.9]);
    
    % Theme Colors
    bgCol   = [1.00 1.00 1.00];  axCol   = [0.98 0.98 0.98]; 
    txtCol  = [0.12 0.12 0.18];  gridCol = [0.88 0.88 0.88];
    accent1 = [0.00 0.45 0.74];  accent2 = [0.85 0.33 0.10]; 
    accent3 = [0.15 0.60 0.35];  accent4 = [0.75 0.00 0.15]; 
    
    set(gcf, 'Color', bgCol);
    styleAx = @(ax) set(ax, 'Color', axCol, 'XColor', txtCol, 'YColor', txtCol, ...
        'GridColor', gridCol, 'GridAlpha', 0.8, 'FontSize', 9, ...
        'TitleFontSizeMultiplier', 1.1, 'FontName', 'Helvetica');

    %% 1. TIME DOMAIN (Signal Overview)
    ax1 = subplot(2, 3, 1); styleAx(ax1);
    t_axis = (0:length(rx)-1)/fs;
    p1 = plot(t_axis, rx, 'Color', [0.8 0.8 0.8]); hold on;
    
    headerBits_tx = headerSize * 2;
    payloadBits_tx = recoveredLen * (1 + useCoding);
    nTxTotalFrames = ceil((headerBits_tx + payloadBits_tx) / (length(dataIdx) * k)) + 1;
    ofdmEnd = min(start + (nTxTotalFrames * symbolLen) - 1, length(rx));
    
    p2 = plot(t_axis(start:ofdmEnd), rx(start:ofdmEnd), 'Color', accent1, 'LineWidth', 1);
    xline(t_axis(start), '--', 'Color', accent4, 'LineWidth', 1.5, 'Label', 'Sync');
    
    title('Time Domain — Frame Capture'); ylabel('Amp'); xlabel('Time (s)');
    legend([p1, p2], {'Ambient', 'OFDM Payload'}, 'Location', 'northeast', 'FontSize', 7);

    %% 2. PASSBAND PSD (Spectral Occupancy)
    ax2 = subplot(2, 3, 2); styleAx(ax2);
    [pxx, f_psd] = pwelch(rx, 1024, 512, 2048, fs);
    plot(f_psd/1000, 10*log10(pxx), 'Color', accent2, 'LineWidth', 1.2); hold on;
    f_start = (fc - bw_theory/2)/1000; f_end = (fc + bw_theory/2)/1000;
    patch([f_start f_end f_end f_start], [min(ylim) min(ylim) max(ylim) max(ylim)], accent1, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
    title('Passband PSD'); ylabel('dB/Hz'); xlabel('kHz'); grid on;

    %% 3. CHANNEL RESPONSE (Denoised Magnitude)
    ax3 = subplot(2, 3, 3); styleAx(ax3);
    h_time = ifft(H_interp);
    h_clean = zeros(size(h_time));
    h_clean(1:min(cplen, length(h_time))) = h_time(1:min(cplen, length(h_time))); 
    H_smooth = fft(h_clean);
    dataFreqs = (dataIdx - nfft/2) * (fs/nfft) / 1000 + fc/1000;
    
    area(dataFreqs, 20*log10(abs(H_smooth)), -100, 'FaceColor', accent3, 'FaceAlpha', 0.1, 'EdgeColor', accent3, 'LineWidth', 1.5);
    title('Denoised Channel $|H(f)|$'); ylabel('Gain (dB)'); xlabel('kHz');
    ylim([mean(20*log10(abs(H_smooth)))-15, mean(20*log10(abs(H_smooth)))+5]); grid on;

    %% 4. IMPULSE RESPONSE (Multipath Profile)
    ax4 = subplot(2, 3, 4); styleAx(ax4);
    ir_est = abs(ifft(H_interp, 256));
    stem((0:127)/fs*1e6, ir_est(1:128), 'Color', accent2, 'Marker', '.', 'LineWidth', 0.5);
    title('Channel Impulse Response'); ylabel('Mag'); xlabel('\mus'); grid on;

    %% 5. ERROR DISTRIBUTION (Where the ICI hides)
    ax5 = subplot(2, 3, 5); styleAx(ax5);
    % We use a bar plot for errors to see specific "trouble" subcarriers
    bar(dataFreqs, serPerSubcarrier, 'FaceColor', accent4, 'EdgeColor', 'none', 'BarWidth', 1);
    title('Spectral Error Density'); ylabel('SER'); xlabel('kHz');
    ylim([0 1]); grid on;
    % Add a line for average BER for comparison
    yline(actualBER, '--', 'Color', [0.3 0.3 0.3], 'Alpha', 0.5);

    %% 6. SYSTEM METRICS PANEL
    ax6 = subplot(2, 3, 6); set(ax6, 'Color', bgCol); axis off;
    metrics = {
        '⚡ Rate',       sprintf('%.2f kbps', R_b/1000);
        '📡 Bandwidth',  sprintf('%.0f Hz', bw_theory);
        '📶 Est. SNR',   sprintf('%.2f dB', snr_est);
        '🧩 Coding',     sprintf('ConstraintLength=%d, r=1/2', constraintLength); % Assuming you pass constraintLength
        '📦 Mod',        modScheme;
        '❌ Final BER',  sprintf('%.5f', actualBER)
    };
    
    text(0.1, 0.95, 'Link Diagnostics', 'Color', txtCol, 'FontSize', 12, 'FontWeight', 'bold', 'Units', 'normalized');
    for i = 1:size(metrics, 1)
        yPos = 0.82 - (i-1)*0.13;
        text(0.1, yPos, metrics{i,1}, 'Color', [0.4 0.4 0.5], 'FontSize', 10, 'Units', 'normalized');
        text(0.55, yPos, metrics{i,2}, 'Color', txtCol, 'FontSize', 10, 'FontWeight', 'bold', 'Units', 'normalized');
    end
end