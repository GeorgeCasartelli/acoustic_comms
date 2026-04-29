clear all; clc; close all;

%{

OFDM SIGNAL RECEIVER

MUST BE MATCHED WITH ofdm_transmit.m

%}

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end

%% --== DEFINE SIGNAL PARAMS ==--

M = 4;
k = log2(M);
nfft = 2048;
cplen = 1024;
fs = 48000;
fc = 10000;

%% --== DEFINE SCRIPT PARAMS ==--

useCoding = true;
txMode = 'text';
imageSize = 128;
headerSize = 32;
useInterleaving = true;

%% --== DEFINE CARRIERS ==--

numActiveCarriers = 400;
pilotSpacing = 5;

% centre around nfft/2
activeCarriers = ((nfft/2) - (numActiveCarriers/2) : (nfft/2) + (numActiveCarriers/2)).';

activeCarriers(activeCarriers == nfft/2 + 1) = []; % remove dc 0 

% define idx for pilot, data, and null carriers
pilotIdx = activeCarriers(1:pilotSpacing:end);
dataIdx = setdiff(activeCarriers, pilotIdx);
nullIdx = setdiff(1:nfft, activeCarriers).'; % define null carrier idx


%% --== MESSAGE ==--

% THIS IS USED FOR BER CALCULATION

file = fopen("./tests/textfile.txt","r")';
msg = fscanf(file, "%c");
fclose(file);

% convert text integers
binchars = dec2bin(msg, 8); 
bits = reshape(binchars.' - '0', [], 1); 
totalBits = numel(bits); 

% format
bitsPerFrame = length(dataIdx) * k;
nFrames = ceil(numel(bits) / bitsPerFrame);

paddedBits = [bits; zeros(nFrames*bitsPerFrame - totalBits, 1)]; % pad to make square. 

bitgroups = reshape(paddedBits, k, [])'; % reshape by width k
inputSymbols = bi2de(bitgroups, 'left-msb'); % conv to int

dataIn = inputSymbols;



%% --== CONTSTELLATION ==--
constSym = pskmod((0:M-1), M, pi/4); 
% define constellation diagram scope
cdScope = comm.ConstellationDiagram( ...
    'SamplesPerSymbol',1,...
    ShowReferenceConstellation=false,...
    ReferenceConstellation=constSym);

%% --== RECORDER ==--
recorder = audiorecorder(fs,16,1);
if strcmp(txMode, 'image')
    recordDuration = 30;
else
    recordDuration = 10;
end
recordblocking(recorder, recordDuration);
rx = getaudiodata(recorder);

% [rx, fs] = audioread("ofdm_signal.wav");

%% --== RX ==--
preamble = mod(0:numActiveCarriers-1, 4).';
preambleSignal = pskmod(preamble, M, pi/4);
preambleData = preambleSignal(1:length(dataIdx));


t_rx = (0:length(rx)-1)' /fs; % new time vector
rx_mixed = rx .* exp(-1j*2*pi*fc*t_rx); % mix down to baseband

% filter
lpCutoff = 10000; % may need adjusting
[b, a] = butter(6, lpCutoff/(fs/2));
rx_baseband = filter(b, a, rx_mixed);

% rx_baseband = rx_mixed;
pilotSym = repmat(pskmod(0, M, pi/4), length(pilotIdx), 1);
preamble_bb = ofdmmod(preambleData, nfft, cplen, nullIdx, pilotIdx, pilotSym); % make a preamble reference

%% --== ALIGN W/ PREAMBLE ==--

% cross correlate received with known preamble
[xc, lags] = xcorr(rx_baseband, preamble_bb);
[~, idx] = max(abs(xc)); 
start = lags(idx);

rxAligned = rx_baseband(start:end); % trim before preamble

% trim down to whole number of ofdm symbols for ofdmdemod
symbolLen = nfft + cplen;
numSymbolsReceived = floor(length(rxAligned) / symbolLen);
rxTrimmed = rxAligned(1 : numSymbolsReceived * symbolLen);

%% --== OFDM DEMOD ==--
[x1, rxPilots] = ofdmdemod(rxTrimmed, nfft, cplen, 0, nullIdx, pilotIdx);
 
rxPreamble = x1(:, 1); % get first symbol (preamble) 
rxDataSyms = x1(:, 2:end); % get payload

% disp(preamble);
% disp(pskdemod(rxPreamble, M, pi/4));

nDataSyms = size(rxDataSyms, 2);
pilots = repmat(pskmod(0,M,pi/4), length(pilotIdx), nDataSyms+1);
nSyms = size(pilots, 2);

x1_equalised = zeros(length(dataIdx), nDataSyms);



for sym = 1:nDataSyms
    H_pilots = rxPilots(:, sym+1) ./ pilots(:, sym+1); % channel estimate at the pilot
    H_interp = interp1(pilotIdx, H_pilots, dataIdx, 'linear', 'extrap'); % interpolate est. to data subc.
    x1_equalised(:, sym) = rxDataSyms(:, sym) ./ H_interp; % apply channel correction
end


rxData = pskdemod(x1_equalised, M, pi/4);



%% ---- DATA CHECK AND RECOVERY ----
% isequal(rxData(:),inputSymbols

% convert symbols back to bits
rawData = de2bi(rxData(:), k, 'left-msb');
allRxBits = reshape(rawData.', [], 1); % flatten to single bitstream

trellis = poly2trellis(3, [ 6 7 ]);
tbdepth = 34;

if useCoding
    % Header is convolutionally encoded at 1/2 rate
    headerCodedLen = headerSize * 2;
    headerPart = allRxBits(1:headerCodedLen);
    decodedHeader = vitdec(headerPart, trellis, 10, 'trunc', 'hard');
    recoveredLen = bi2de(decodedHeader(1:headerSize)', 'left-msb');
else
    % Header is raw bits, just read directly
    headerPart = allRxBits(1:headerSize);
    recoveredLen = bi2de(headerPart', 'left-msb');
end

fprintf("Header decoded! Payload length: %d bits\n", recoveredLen);

if useCoding
    payloadCodedLen = recoveredLen * 2;
    payloadPart = allRxBits(headerCodedLen + 1 : headerCodedLen + payloadCodedLen);

    if useInterleaving
        payloadPart = randdeintrlv(payloadPart, 12345);
    end

    finalBits = vitdec(payloadPart, trellis, 34, 'trunc', 'hard');
    finalBits = finalBits(1:recoveredLen);
else
    % Payload starts immediately after the raw header
    payloadPart = allRxBits(headerSize + 1 : headerSize + recoveredLen);

    if useInterleaving
        payloadPart = randdeintrlv(payloadPart, 12345);
    end

    finalBits = payloadPart;
end

% remodulate for constellation plot 
if length(allRxBits) < headerSize
    disp('Error: Could not even decode the 16-bit header.');
    return;
end

%% --== DE INTERLEAV ==--
% if useInterleaving
%     N = length(allRxBits);
% 
%     rng(42);
%     elements = randperm(N);
% 
%     allRxBits = deintrlv(allRxBits, elements);
% 
% end


%% --==DECODE ==--
% if useCoding
%     decodedBits = vitdec(allRxBits, trellis, tbdepth, 'trunc', 'hard');    
% else
%     decodedBits = allRxBits;
% end
% 
% recoveredLen = bi2de(decodedBits(1:headerSize)', 'left-msb'); % get msg length
% fprintf('Header decoded! Message length: %d bits\n', recoveredLen);
% 
% payloadBits = decodedBits(headerSize+1:end);
% 
% if length(payloadBits) < (recoveredLen)
%     fprintf("Warning: only got %d of %d bits\n", length(payloadBits), recoveredLen);
%     recoveredLen = length(payloadBits);
% end
% 
% finalBits = payloadBits(1: recoveredLen); %skip header, grab data
% 
if strcmp(txMode, 'image')
    rawBytes = bi2de(reshape(finalBits, 8, []).', 'left-msb');
    recoveredImg = uint8(reshape(rawBytes, imageSize, imageSize));
    imshow(recoveredImg);
else
    % inverse of tx, convert back to text
    outputBits_reshape = reshape(finalBits, 8, [])';
    charArray = char(outputBits_reshape + '0');
    outputText = bin2dec(charArray).';
    fprintf('Recovered: %s\n', char(outputText));
end

% print bit error
if strcmp(txMode, 'text')
    if length(bits) == length(finalBits)     
        bitErrors = sum(bits ~= finalBits);
        actualBER = bitErrors / recoveredLen;
        fprintf('Bit Error Rate (BER): %.6f\n', actualBER);
    else
        fprintf("Unmatching array sizes bits: $d , finalBits: $d\n", length(bits), length(finalBits));
    end

else 
    fprintf("Receiving unknown symbol. Cannot get BER\n");
end
%% --== PLOTS ==--
framesNeeded = ceil((recoveredLen + headerSize) / (length(dataIdx) * k));
validSyms = x1_equalised(:, 1:framesNeeded);
cdScope(validSyms(:));

[R_b] = data_rate_calc(fs, nfft, cplen, length(dataIdx), 0.5, k);


%% --== DASHBOARD ==--
figure('Name', 'OFDM Demo Day Analysis', 'Units', 'normalized', 'OuterPosition', [0.05 0.05 0.9 0.9]);

% Shared dark theme colours
bgCol   = [0.12 0.12 0.16];
axCol   = [0.18 0.18 0.24];
txtCol  = [0.95 0.95 0.95];
gridCol = [0.30 0.30 0.38];
accent1 = [0.26 0.63 0.95];  % blue
accent2 = [0.95 0.45 0.15];  % orange
accent3 = [0.25 0.88 0.58];  % green
accent4 = [0.95 0.30 0.45];  % red

set(gcf, 'Color', bgCol);

%% -- Helper to style axes --
styleAx = @(ax) set(ax, ...
    'Color', axCol, 'XColor', txtCol, 'YColor', txtCol, ...
    'GridColor', gridCol, 'GridAlpha', 0.5, 'FontSize', 10, ...
    'TitleFontSizeMultiplier', 1.1);

%% -- 1. TIME DOMAIN --
ax1 = subplot(2, 3, 1); styleAx(ax1);
t_axis = (0:length(rx)-1)/fs;

% Calculate exact number of transmitted symbols
% Header is always coded (2x), payload depends on useCoding
headerBits_tx = headerSize * 2;
if useCoding
    payloadBits_tx = recoveredLen * 2;
else
    payloadBits_tx = recoveredLen;
end
totalTxBits = headerBits_tx + payloadBits_tx;

% How many OFDM data frames does that occupy (add 1 for preamble)
bitsPerFrame   = length(dataIdx) * k;
nTxDataFrames  = ceil(totalTxBits / bitsPerFrame);
nTxTotalFrames = nTxDataFrames + 1; % +1 for preamble symbol

% Total samples in the OFDM block
symbolLen    = nfft + cplen;
ofdmBlockLen = nTxTotalFrames * symbolLen;
ofdmEnd      = min(start + ofdmBlockLen - 1, length(rx));

% Full signal in dim grey
plot(t_axis, rx, 'Color', [0.45 0.45 0.55], 'LineWidth', 0.5); hold on;

% Highlight the detected OFDM block in blue
t_ofdm = t_axis(start:ofdmEnd);
plot(t_ofdm, rx(start:ofdmEnd), 'Color', accent1, 'LineWidth', 1.2);

% Mark start with a vertical line
xline(t_axis(start), '--', 'Color', accent3, 'LineWidth', 1.2, 'Label', 'Sync');

title('Time Domain — Received Signal', 'Color', txtCol);
xlabel('Time (s)', 'Color', txtCol);
ylabel('Amplitude', 'Color', txtCol);
legend({'Raw RX', 'OFDM Block'}, 'TextColor', txtCol, 'Color', axCol, 'EdgeColor', gridCol);
grid on;

%% -- 2. PASSBAND PSD --
ax2 = subplot(2, 3, 2); styleAx(ax2);
[pxx, f_psd] = pwelch(rx, 1024, 512, 2048, fs);
plot(f_psd/1000, 10*log10(pxx), 'Color', accent2, 'LineWidth', 1.4); hold on;

bw_theory = (numActiveCarriers * fs) / nfft;
f_start   = (fc - bw_theory/2) / 1000;
f_end     = (fc + bw_theory/2) / 1000;
ylims     = ylim;
patch([f_start f_end f_end f_start], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
      accent1, 'FaceAlpha', 0.15, 'EdgeColor', accent1, 'LineWidth', 1);

% mark carrier and band edges explicitly
xline(fc/1000,           '-',  'Color', txtCol,  'Alpha', 0.5, ...
      'Label', sprintf('f_c = %g kHz', fc/1000),       'FontSize', 8);
xline(f_start,           '--', 'Color', accent1, 'Alpha', 0.7, ...
      'Label', sprintf('%.1f kHz', f_start),            'FontSize', 8);
xline(f_end,             '--', 'Color', accent1, 'Alpha', 0.7, ...
      'Label', sprintf('%.1f kHz', f_end),              'FontSize', 8);

title(sprintf('Passband PSD  |  OBW: %.0f Hz', bw_theory), 'Color', txtCol);
xlabel('Frequency (kHz)', 'Color', txtCol);
ylabel('Power (dB/Hz)',   'Color', txtCol);
xlim([(fc - fs/4)/1000,  (fc + fs/4)/1000]);
legend({'PSD', 'Occupied BW'}, 'TextColor', txtCol, 'Color', axCol, 'EdgeColor', gridCol);
grid on;

%% -- 3. CHANNEL FREQUENCY RESPONSE --
ax3 = subplot(2, 3, 3); styleAx(ax3);
H_mag_dB = 20*log10(abs(H_interp));

% convert subcarrier index to actual frequency in kHz
dataFreqs = (dataIdx - nfft/2) * (fs/nfft) / 1000 + fc/1000;
plot(dataFreqs, H_mag_dB, 'Color', accent3, 'LineWidth', 1.3);
title('Channel Freq. Response', 'Color', txtCol);
xlabel('Frequency (kHz)', 'Color', txtCol);
ylabel('|H| (dB)',        'Color', txtCol);
grid on;

%% -- 4. CHANNEL IMPULSE RESPONSE --
ax4 = subplot(2, 3, 4); styleAx(ax4);
ir_est  = abs(ifft(H_interp, 256));
ir_show = ir_est(1:48);
t_ir_us = (0:47) / fs * 1e6; % convert samples to microseconds
stem(t_ir_us, ir_show, 'Color', accent2, 'MarkerFaceColor', accent2, ...
     'MarkerSize', 4, 'LineWidth', 1.5);
title('Channel Impulse Response', 'Color', txtCol);
xlabel('Delay (\mus)',  'Color', txtCol);
ylabel('Magnitude',    'Color', txtCol);
grid on;

%% -- 5. SNR PER PILOT --
ax5 = subplot(2, 3, 5); styleAx(ax5);
nReceivedFrames  = size(rxPilots, 2);
idealPilotsFull  = repmat(pilotSym, 1, nReceivedFrames);
errors           = rxPilots - idealPilotsFull;

snr_per_pilot = 10*log10( mean(abs(rxPilots).^2, 2) ./ ...
                           max(mean(abs(errors).^2,  2), 1e-10) );

plot(pilotIdx, snr_per_pilot, 'Color', accent1, 'LineWidth', 1.3); hold on;
yline(mean(snr_per_pilot), '--', 'Color', accent4, 'LineWidth', 1.2, ...
      'Label', sprintf('Mean %.1f dB', mean(snr_per_pilot)), ...
      'LabelHorizontalAlignment', 'left', 'FontSize', 9);

title('SNR per Pilot Subcarrier', 'Color', txtCol);
xlabel('Subcarrier Index', 'Color', txtCol);
ylabel('SNR (dB)',         'Color', txtCol);
grid on;

% Overall scalar SNR (used in metrics panel)
P_signal = mean(abs(rxPilots(:)).^2);
P_noise  = mean(abs(errors(:)).^2);
if P_noise == 0; P_noise = 1e-10; end
snr_est  = 10*log10(P_signal / P_noise);

%% -- 6. METRICS PANEL --
ax6 = subplot(2, 3, 6);
set(ax6, 'Color', axCol, 'XColor', 'none', 'YColor', 'none'); axis off;

% Card-style rows
metrics = {
    '⚡  Data Rate',         sprintf('%.2f kbps',  R_b/1000);
    '📡  Occupied BW',       sprintf('%.0f Hz',    bw_theory);
    '📶  Est. SNR',          sprintf('%.2f dB',    snr_est);
    '🔢  Spectral Eff.',     sprintf('%.3f b/s/Hz',R_b/bw_theory);
    '🧩  Coding',            ternary(useCoding, 'Conv. r=1/2', 'None');
    '🔀  Interleaving',      ternary(useInterleaving, 'On', 'Off');
    '📦  Mode',              txMode;
};

if strcmp(txMode, 'text')
    metrics{end+1, 1} = '❌  BER';
    metrics{end,   2} = sprintf('%.5f', actualBER);
end

nRows  = size(metrics, 1);
yStart = 0.93;
yStep  = 0.88 / nRows;

text(0.0, 1.01, 'System Metrics', 'Color', txtCol, 'FontSize', 13, ...
     'FontWeight', 'bold', 'Units', 'normalized');

for i = 1:nRows
    yPos = yStart - (i-1)*yStep;
    % dim label
    text(0.02, yPos, metrics{i,1}, 'Color', [0.65 0.65 0.75], ...
         'FontSize', 10, 'Units', 'normalized');
    % bright value
    text(0.58, yPos, metrics{i,2}, 'Color', txtCol, ...
         'FontSize', 10, 'FontWeight', 'bold', 'Units', 'normalized');
end