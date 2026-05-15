clear all; clc; close all;

%{

OFDM SIGNAL RECEIVER

MUST BE MATCHED WITH ofdm_transmit.m

%}

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end

addpath("results_loggers_plotters");

%% --== DEFINE SCRIPT PARAMS ==--

useCoding = true;
txMode = 'image';
imageSize = 128;
headerSize = 32;
useInterleaving = false;
modScheme = "qpsk";
constraintLength = 7;

phyParams.test = 'wall-reflections';
phyParams.trial = 2;

phyParams.distance = 1.5;
phyParams.volume = 100;

phyParams.txAngle = 30;
phyParams.rxAngle = -30;

phyParams.blocked = 0;
phyParams.blockerType = "none";

phyParams.reflectionMode = "wall-only-45";

phyParams.modscheme = modScheme;
phyParams.useCoding = useCoding;

phyParams.notes = "unblocked_reflection";

%% --== DEFINE SIGNAL PARAMS ==--

if strcmp(modScheme, "qpsk") 
    M = 4;
elseif strcmp(modScheme, "16qam")
    M = 16;
end
k = log2(M);
nfft = 2048;
cplen = 256;
fs = 48000;
fc = 10000;


%% --== DEFINE CARRIERS ==--

numActiveCarriers = 500;
pilotSpacing = 5;

% centre around nfft/2
activeCarriers = ((nfft/2) - (numActiveCarriers/2) : (nfft/2) + (numActiveCarriers/2)).';

activeCarriers(activeCarriers == nfft/2 + 1) = []; % remove dc 0 

% define idx for pilot, data, and null carriers
pilotIdx = activeCarriers(1:pilotSpacing:end);
dataIdx = setdiff(activeCarriers, pilotIdx);
nullIdx = setdiff(1:nfft, activeCarriers).'; % define null carrier idx


%% --== MESSAGE (FOR BER) ==--

file = fopen("./tests/textfile.txt","r")';
msg = fscanf(file, "%c");
fclose(file);
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
if strcmp(modScheme, "qpsk")
    constSym = pskmod((0:M-1), M, pi/4); 
elseif strcmp(modScheme, "16qam")
    constSym = qammod((0:M-1), M, "UnitAveragePower", true);
end

% define constellation diagram scope
cdScope = comm.ConstellationDiagram( ...
    'SamplesPerSymbol',1,...
    ShowReferenceConstellation=false,...
    ReferenceConstellation=constSym);

%% --== RECORDER ==--
recorder = audiorecorder(fs,16,1);
recordDuration = ternary(strcmp(txMode, 'image'), 20, 15);
recordblocking(recorder, recordDuration);
rx = getaudiodata(recorder);

% [rx, fs] = audioread("ofdm_signal.wav");

%% --== RX ==--
rng(42); 

% Generate a random integer sequence for the preamble
preamble = randi([0 M-1], numActiveCarriers, 1);

% preamble = mod(0:numActiveCarriers-1, M).';
if strcmp(modScheme, "qpsk")
    preambleSignal = pskmod(preamble, M, pi/4);
elseif strcmp(modScheme, "16qam")
    preambleSignal = qammod(preamble, M, "UnitAveragePower",true);
end
preambleData = preambleSignal(1:length(dataIdx));

t_rx = (0:length(rx)-1)' /fs; % new time vector
rx_mixed = rx .* exp(-1j*2*pi*fc*t_rx); % mix down to baseband

% filter
lpCutoff = 10000; % may need adjusting
[b, a] = butter(6, lpCutoff/(fs/2));
rx_baseband = filter(b, a, rx_mixed);

% rx_baseband = rx_mixed;

if strcmp(modScheme, "qpsk")
    pilotSymbol = pskmod(0,M,pi/4);
elseif strcmp(modScheme, "16qam")
    pilotSymbol = qammod(0, M, "UnitAveragePower",true);
end


pilotSyms = repmat(pilotSymbol, length(pilotIdx), 1);


preamble_bb = ofdmmod(preambleData, nfft, cplen, nullIdx, pilotIdx, pilotSyms); % make a preamble reference

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

%% Plot Correlation Statistic (Synchronization Diagnostic)
figure('Name', 'Receiver Synchronization');
hold on; grid on;

% We only care about the magnitude of the correlation
correlationMag = abs(xc);

% Plot the full correlation result
plot(lags / fs * 1000, correlationMag, 'Color', [0 0.45 0.74], 'LineWidth', 1.5);

% Mark the detected peak with a hollow circle (consistent with your style!)
[maxVal, maxIdx] = max(correlationMag);
peakLag = lags(maxIdx) / fs * 1000; % Convert to ms
plot(peakLag, maxVal, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Detected Start');

% Formatting
title('Synchronization: Cross-Correlation Statistic');
xlabel('Time Lag (ms)');
ylabel('Correlation Magnitude');
legend('Correlation Result', 'Detected Sync Point');

% Zoom in on the area around the peak for the report
xlim([peakLag - 10, peakLag + 10]);

%% --== OFDM DEMOD ==--
[x1, rxPilots] = ofdmdemod(rxTrimmed, nfft, cplen, 0, nullIdx, pilotIdx);
rxPreamble = x1(:, 1); % get first symbol (preamble) 
rxDataSyms = x1(:, 2:end); % get payload

nDataSyms = size(rxDataSyms, 2);
pilots = repmat(pilotSymbol, length(pilotIdx), nDataSyms+1);
nSyms = size(pilots, 2);

x1_equalised = zeros(length(dataIdx), nDataSyms);
for sym = 1:nDataSyms
    H_pilots = rxPilots(:, sym+1) ./ pilots(:, sym+1); % channel estimate at the pilot
    H_interp = interp1(pilotIdx, H_pilots, dataIdx, 'linear', 'extrap'); % interpolate est. to data subc.
    x1_equalised(:, sym) = rxDataSyms(:, sym) ./ H_interp; % apply channel correction
end


if strcmp(modScheme, "qpsk")
    rxData =pskdemod(x1_equalised, M, pi/4);
else
    rxData = qamdemod(x1_equalised, M, "UnitAveragePower", true);
end
%% --== DECODE AND DE INTERLEAVE ==--

% convert symbols back to bits
rawData = de2bi(rxData(:), k, 'left-msb');
allRxBits = reshape(rawData.', [], 1); % flatten to single bitstream

if constraintLength == 3
    trellis = poly2trellis(constraintLength, [ 6 7 ]);
elseif constraintLength == 7
    trellis = poly2trellis(7, [171 133]);
end

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

% skip header decode for tests
% recoveredLen = 28672;
recoveredLen = 131072;

interleaveSeed = 12345;

fprintf("Header decoded! Payload length: %d bits\n", recoveredLen);

if useCoding
    payloadCodedLen = recoveredLen * 2;
    payloadPart = allRxBits(headerCodedLen + 1 : headerCodedLen + payloadCodedLen);

    if useInterleaving
        payloadPart = randdeintrlv(payloadPart, interleaveSeed);
    end

    finalBits = vitdec(payloadPart, trellis, 34, 'trunc', 'hard');
    finalBits = finalBits(1:recoveredLen);
else
    % Payload starts immediately after the raw header
    payloadPart = allRxBits(headerSize + 1 : headerSize + recoveredLen);

    if useInterleaving
        payloadPart = randdeintrlv(payloadPart, interleaveSeed);
    end

    finalBits = payloadPart;
end

% remodulate for constellation plot 
if length(allRxBits) < headerSize
    disp('Error: Could not even decode the 16-bit header.');
    return;
end

%% --== DATA REFORMAT ==--

recoveredImg = [];

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

%% --== BER ==--

actualBER = 0.0;

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

[R_b] = data_rate_calc(fs, nfft, cplen, length(dataIdx), 0.5 + 0.5 * ~useCoding, k);

%% --== SER ==--
if strcmp(txMode, 'text')
    txBits = [convenc(de2bi(totalBits, 32, 'left-msb').', trellis); ...
              randintrlv(convenc(bits, trellis), interleaveSeed)];
    
    if ~useCoding
        txBits = [de2bi(totalBits, 32, 'left-msb').'; randintrlv(bits, interleaveSeed)];
    end
    if ~useInterleaving
        txBits = [convenc(de2bi(totalBits, 32, 'left-msb').', trellis); convenc(bits, trellis)];
    end
    if ~useCoding && ~useInterleaving
        txBits = [de2bi(totalBits, 32, 'left-msb').'; bits];
    end

    bitsPerFrame = length(dataIdx) * k;
    txBits = [txBits; zeros(ceil(numel(txBits)/bitsPerFrame)*bitsPerFrame - numel(txBits), 1)];
    idealDataSymbols = reshape(bi2de(reshape(txBits, k, []).', 'left-msb'), length(dataIdx), []);

    compareLen = min(size(rxData, 2), size(idealDataSymbols, 2));
    rxTest = rxData(:, 1:compareLen);
    txTest = idealDataSymbols(:, 1:compareLen);

    serPerSubcarrier = sum(rxTest ~= txTest, 2) / compareLen;
    actualSER = sum(rxTest(:) ~= txTest(:)) / numel(rxTest);
    fprintf('SER: %.6f\n', actualSER);
else
    serPerSubcarrier = zeros(length(dataIdx), 1);
    actualSER = NaN;
    fprintf('SER: N/A (image mode)\n');
end

%% --== SNR ESTIMATION ==--
if strcmp(modScheme, 'qpsk')
    idealRef = pskmod((0:M-1), M, pi/4);
elseif strcmp(modScheme, '16qam')
    idealRef = qammod((0:M-1), M, 'UnitAveragePower', true);
end

syms    = validSyms(:);
[~, idx] = min(abs(syms - idealRef(:).'), [], 2);
nearest  = idealRef(idx);
snr_est  = 10 * log10(mean(abs(nearest).^2) / mean(abs(syms - nearest(:)).^2));

%% --== DASHBOARD ==--

display_ofdm_dashboard(rx, fs, start, recoveredLen, headerSize, dataIdx, k, ...
    fc, numActiveCarriers, H_interp, rxPilots, pilotSyms, pilotIdx, nfft, cplen, ...
    useCoding, useInterleaving, txMode, actualBER, R_b, modScheme, snr_est, ...
    constraintLength, serPerSubcarrier, validSyms(:), ternary(strcmp(txMode,'image'), recoveredImg, []));

%% --== LOG TEST (PARAM TUNING) ==--

params.test           = 'distance_sweep_qpsk';
params.M              = M;
params.cplen          = cplen;
params.pilotSpacing   = pilotSpacing;
params.numActiveCarriers = numActiveCarriers;
params.useCoding      = useCoding;
params.fc             = fc;
params.distance       = 3.0;
params.volume         = 150; % percent
params.constraintLength = constraintLength;
params.modscheme     = modScheme;

% log_result('results.csv', params, actualBER, snr_est, R_b);

%% --== LOG TEST (PHY EVAL) ==--

% phyParams.test = 'distance_qpsk';
% phyParams.trial = 3;
% 
% phyParams.distance = 0.5;
% phyParams.volume = 100;
% 
% phyParams.txAngle = 0;
% phyParams.rxAngle = 0;
% 
% phyParams.blocked = 0;
% phyParams.blockerType = "none";
% 
% phyParams.reflectionMode = "direct";
% 
% phyParams.modscheme = modScheme;
% phyParams.useCoding = useCoding;
% 
% phyParams.notes = "clear_los";

% log_physical_results( ...
%     './results_loggers_plotters/physical_results.csv', ...
%     phyParams, ...
%     actualBER, ...
%     actualSER, ...
%     snr_est, ...
%     R_b);