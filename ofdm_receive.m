clear all; clc;

%{

OFDM SIGNAL RECEIVER

MUST BE MATCHED WITH ofdm_transmit.m

%}


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

tbdepth = 34;
trellis = poly2trellis(3, [ 6 7 ]);

% remodulate for constellation plot 
if length(allRxBits) < headerSize
    disp('Error: Could not even decode the 16-bit header.');
    return;
end


if useCoding
    decodedBits = vitdec(allRxBits, trellis, tbdepth, 'trunc', 'hard');    
    recoveredLen = bi2de(decodedBits(1:headerSize)', 'left-msb'); % get msg length
    fprintf('Header decoded! Message length: %d bits\n', recoveredLen);
    payloadBits = decodedBits(headerSize+1:end);
else
    recoveredLen = bi2de(allRxBits(1:headerSize)', 'left-msb');
    payloadBits = allRxBits(headerSize+1:end);
    fprintf('Header decoded! Message length: %d bits\n', recoveredLen);
end

if length(payloadBits) < (recoveredLen)
    fprintf("Warning: only got %d of %d bits\n", length(payloadBits), recoveredLen);
    recoveredLen = length(payloadBits);
end

finalBits = payloadBits(1: recoveredLen); %skip header, grab data
    
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
if length(bits) == length(finalBits)     
    bitErrors = sum(bits ~= finalBits);
    actualBER = bitErrors / recoveredLen;
    fprintf('Bit Error Rate (BER): %.4f\n', actualBER);
else
    fprintf("Unmatching array sizes bits: $d , finalBits: $d\n", length(bits), length(finalBits));
end

%% --== PLOTS ==--
framesNeeded = ceil((recoveredLen + headerSize) / (length(dataIdx) * k));
validSyms = x1_equalised(:, 1:framesNeeded);
cdScope(validSyms(:));

[R_b, eff] = data_rate_calc(fs, nfft, cplen, length(dataIdx), numActiveCarriers, 0.5, k);


%% --== DASHBOARD ==--

figure('Name', 'OFDM Receiver Dashboard', 'Position', [100, 100, 1400, 900]);

% 1. Cross correlation (preamble sync)
subplot(2, 3, 1);
plot(lags, abs(xc));
xlabel('Lag (samples)');
ylabel('|Correlation|');
title('Preamble Cross-Correlation');
xline(start, 'r--', 'Sync Point');
grid on;

% 2. Frequency spectrum of received baseband
% spectrum gets its own wide plot at the top
subplot(2, 3, [1 2]); % spans two columns
nfft_plot = 4096;
f = linspace(-fs/2, fs/2, nfft_plot);
Rx_fft = fftshift(fft(rx_baseband, nfft_plot));
plot(f/1000, 20*log10(abs(Rx_fft)));
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Received Signal Spectrum');
ylim([max(20*log10(abs(Rx_fft))) - 60, max(20*log10(abs(Rx_fft))) + 5]);

carrierFreqs = (dataIdx - nfft/2) * (fs/nfft) + fc;
pilotFreqs = (pilotIdx - nfft/2) * (fs/nfft) + fc;

xline(carrierFreqs/1000, 'b', 'Alpha', 0.1, 'LineWidth', 0.5);
xline(pilotFreqs/1000, 'r', 'Alpha', 0.3, 'LineWidth', 0.8);

legend('Spectrum', 'Data carriers', 'Pilot carriers', 'Location', 'northwest');
grid on;

% 3. Received signal time domain
subplot(2, 3, 3);
t_plot = (0:length(rx)-1)/fs;
plot(t_plot, rx);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Audio (Time Domain)');
xline(start/fs, 'r--', 'Sync Point');
grid on;

% 4. Channel estimate magnitude across subcarriers
subplot(2, 3, 4);
H_mean = zeros(length(pilotIdx), 1);
for sym = 1:nDataSyms
    H_mean = H_mean + abs(rxPilots(:, sym+1) ./ pilots(:, sym+1));
end
H_mean = H_mean / nDataSyms;
plot(pilotIdx, H_mean, 'o-');
xlabel('Subcarrier Index');
ylabel('|H|');
title('Mean Channel Estimate at Pilots');
grid on;

% 5. EVM per subcarrier
subplot(2, 3, 5);
idealSyms = pskmod(pskdemod(x1_equalised(:, 1:framesNeeded), M, pi/4), M, pi/4);
evm_per_carrier = mean(abs(x1_equalised(:, 1:framesNeeded) - idealSyms).^2, 2);
plot(dataIdx, 10*log10(evm_per_carrier));
xlabel('Subcarrier Index');
ylabel('EVM (dB)');
title('Error Vector Magnitude per Subcarrier');
grid on;

% 6. BER summary text box
subplot(2, 3, 6);
axis off;
if length(bits) == length(finalBits)
    berText = sprintf('%.4f', actualBER);
else
    berText = 'N/A';
end
text(0.5, 0.9, 'System Summary', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');
text(0.5, 0.82, sprintf("Data Rate: %d", R_b), ...
    "FontSize", 11, "HorizontalAlignment", "center");
text(0.5, 0.70, sprintf('Active Carriers: %d', numActiveCarriers), ...
    'FontSize', 11, 'HorizontalAlignment', 'center');
text(0.5, 0.58, sprintf('Pilot Spacing: %d', pilotSpacing), ...
    'FontSize', 11, 'HorizontalAlignment', 'center');
text(0.5, 0.46, sprintf('Data Carriers: %d', length(dataIdx)), ...
    'FontSize', 11, 'HorizontalAlignment', 'center');
text(0.5, 0.34, sprintf('Coding: %s', mat2str(useCoding)), ...
    'FontSize', 11, 'HorizontalAlignment', 'center');
text(0.5, 0.22, sprintf('BER: %s', berText), ...
    'FontSize', 11, 'HorizontalAlignment', 'center');
text(0.5, 0.10, sprintf('Recovered %d / %d bits', length(finalBits), recoveredLen), ...
    'FontSize', 11, 'HorizontalAlignment', 'center');
