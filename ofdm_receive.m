clear all; clc;

%{

OFDM SIGNAL RECEIVER

MUST BE MATCHED WITH ofdm_transmit.m

%}

showPlots = true;

%% --== DEFINE SIGNAL PARAMS ==--

M = 4;
k = log2(M);
nfft = 2048;
cplen = 512;
fs = 48000;
fc = 10000;

%% --== DEFINE CARRIERS ==--

numActiveCarriers = 100;
pilotSpacing = 5;

% centre around nfft/2
% activeCarriers = ((nfft/2) - (numActiveCarriers/2) : (nfft/2) + (numActiveCarriers/2)).';
% activeCarriers(activeCarriers == nfft/2 + 1) = []; % remove dc 0 

halfCarriers = numActiveCarriers/2;
posCarriers = 2 : (halfCarriers + 1);
negCarriers = (nfft - halfCarriers + 1) : nfft;

activeCarriers = [ posCarriers, negCarriers].';

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

% create pilot data ( symbol 00 ) based on nFrames
% pilots = repmat(pskmod(0,M,pi/4),length(pilotIdx),nFrames+1);


%% --== CONTSTELLATION ==--
constSym = pskmod((0:M-1), M, pi/4); 
% define constellation diagram scope
cdScope = comm.ConstellationDiagram( ...
    'SamplesPerSymbol',1,...
    ShowReferenceConstellation=false,...
    ReferenceConstellation=constSym);

%% --== RECORDER ==--
recorder = audiorecorder(fs,16,1);
recordDuration = 5;
recordblocking(recorder, recordDuration);
rx = getaudiodata(recorder);

% [rx, fs] = audioread("ofdm_signal.wav");


%% --== RX ==--
preamble = mod(0:numActiveCarriers-1, 4).';
preambleSignal = pskmod(preamble, M, pi/4);
preambleData = preambleSignal(1:length(dataIdx));


t_rx = (0:length(rx)-1)' /fs; % new time vector
rx_bb = rx .* exp(-1j*2*pi*fc*t_rx); % mix down to baseband

% filterdata
lpCutoff = 6000; % may need adjusting
[b, a] = butter(6, lpCutoff/(fs/2));
rx_bb = filter(b, a, rx_bb);

% rx_baseband = rx_mixed;
pilotRef = pskmod(0,M,pi/4) * ones(length(pilotIdx), 1);
preamble_bb = ofdmmod(preambleData, nfft, cplen, nullIdx, pilotIdx, pilotRef); % make a preamble reference

%% --== SCHMIDL COX DECODE ==--
symbolLen = nfft + cplen;
P = zeros(length(rx_bb) - nfft, 1);
E = zeros(length(rx_bb) - nfft, 1);

for n = 1:length(rx_bb) - nfft - 1
    window = rx_bb(n : n + nfft/2 - 1);
    nextWindow = rx_bb(n + nfft/2 : n + nfft - 1);
    P(n) = sum(conj(window) .* nextWindow);
    E(n) = sum(abs(nextWindow).^2);
end

% stop from dividing by 0
timingMetric = zeros(size(P));
energyThreshold = 0.05 * max(E.^2); % only calc when signal exists
validIdx = (E.^2) > energyThreshold; 
timingMetric(validIdx) = abs(P(validIdx)).^2 ./ E(validIdx).^2;

smoothingWindow = floor(cplen/2);
timingMetricSmooth = movmean(timingMetric, smoothingWindow);

threshold = 0.9 * max(timingMetric);
above = timingMetric > threshold;
risingEdges = find(diff(above) == 1);
% add this to rx temporarily
fallingEdges = find(diff(above) == -1);

if isempty(risingEdges) || isempty(fallingEdges)
    disp("No valid S&C preamble");
    return;
end

safetyMargin = 0;
scIdx = risingEdges(1) + cplen - safetyMargin;
% deltaF = 0;
deltaF = angle(P(scIdx)) * fs / (pi*nfft);
fprintf('CFO estimate: %.2f Hz\n', deltaF);

t_rx = (0:length(rx_bb)-1)' / fs;
rx_bb = rx_bb .* exp(-1j * 2 * pi * deltaF * t_rx);
% 



dataStart = scIdx + nfft;

% % searchRange = scIdx : scIdx + (nfft + cplen) * 2;
% % [xc, lags] = xcorr(rxCorrected(searchRange), preamble_bb);
% % [~, localIdx] = max(abs(xc));
% dataStart = scIdx + (nfft+cplen);
% 
rxAligned = rx_bb(dataStart:end);
% 
% % trim to be whole number of ofdm symbols
% symbolLen = nfft + cplen;
numSymbolsReceived = floor(length(rxAligned) / symbolLen);
% rxTrimmed = rxAligned(1 : numSymbolsReceived * symbolLen);
rxTrimmed = rxAligned(1:numSymbolsReceived*symbolLen);

fprintf('Symbols received: %d / %d expected\n', numSymbolsReceived, nFrames);

fprintf('Rising edge: %d\n', risingEdges(1));
fprintf('Falling edge: %d\n', fallingEdges(1));
fprintf('Plateau width: %d\n', fallingEdges(1) - risingEdges(1));
fprintf('Symbols received: %d / %d\n', numSymbolsReceived, nFrames);

% reconstruct what rx_bb should look like at dataStart
% compare energy of each received symbol
figure;
symbolEnergies = zeros(numSymbolsReceived, 1);
for i = 1:numSymbolsReceived
    chunk = rxTrimmed((i-1)*symbolLen+1 : i*symbolLen);
    symbolEnergies(i) = mean(abs(chunk).^2);
end
plot(symbolEnergies);
title('Energy per received symbol - should be flat');

% also check first symbol CP vs end of symbol (should match for correct alignment)
sym1 = rxTrimmed(1:symbolLen);
cp = sym1(1:cplen);
tail = sym1(end-cplen+1:end);
fprintf('CP match error: %.6f\n', norm(cp - tail));

if showPlots
    figure;
    t_metric = (0:length(timingMetric)-1) / fs;
    subplot(2,1,1);
    plot(t_metric, timingMetric);
    title('Schmidl-Cox Timing Metric');
    xlabel('Time (s)');
    xline(scIdx/fs, 'r--', 'scIdx');

    subplot(2,1,2);
    plotRange = max(1,scIdx-700) : min(length(rx_bb), scIdx+symbolLen*3);
    plot(real(rx_bb(plotRange)));
    xline(701, 'r--', 'scIdx');
    xline(701+(nfft), 'g--', 'dataStart');
    xline(701 + (nfft/2), 'p--', 's&c middle')
    title('Baseband signal around sync point');
end

%% --== ALIGN W/ PREAMBLE ==--

% % cross correlate received with known preamble
% [xc, lags] = xcorr(rxCorrected, preamble_bb);
% [~, idx] = max(abs(xc)); 
% start = lags(idx) + length(preamble_bb);
% 
% rxAligned = rxCorrected(start:end); % trim before preamble
% 
% % trim down to whole number of ofdm symbols for ofdmdemod
% symbolLen = nfft + cplen;
% numSymbolsReceived = floor(length(rxAligned) / symbolLen);
% rxTrimmed = rxAligned(1 : numSymbolsReceived * symbolLen);
% 
% 
% % rxAligned = rx_baseband;
% 

%% --== OFDM DEMOD ==--
[x1, rxPilots] = ofdmdemod(rxTrimmed, nfft, cplen, 0, nullIdx, pilotIdx);

pilots = repmat(pskmod(0,M,pi/4), length(pilotIdx), size(rxPilots,2));

% disp(preamble);
% disp(pskdemod(rxPreamble, M, pi/4));

% nSyms = size(pilots, 2);
nDataSyms = size(x1, 2);
x1_equalised = zeros(length(dataIdx), nDataSyms);

for sym = 1:nDataSyms
    H_pilots = rxPilots(:, sym) ./ pilots(:, sym); % channel estimate at the pilot
    H_interp = interp1(pilotIdx, H_pilots, dataIdx, 'linear', 'extrap'); % interpolate est. to data subc.
    x1_equalised(:, sym) = x1(:, sym) ./ H_interp; % apply channel correction
end

rxData = pskdemod(x1_equalised, M, pi/4);
cdScope(x1_equalised(:)); % plot on scope

%% ---- DATA CHECK AND RECOVERY ----

% convert symbols back to bits
rawData = de2bi(rxData(:), k, 'left-msb');
allRxBits = reshape(rawData.', [], 1); % flatten to single bitstream

if length(allRxBits) < 16
    disp('Error: not enough bits for header');
    return;
end

recoveredLen = bi2de(allRxBits(1:16)', 'left-msb'); % get msg length
% recoveredLen = totalBits; % bypass header, use known length
fprintf('Header decoded! Message length: %d bits\n', recoveredLen);
    
payloadBits = allRxBits(17:end);

if length(payloadBits) < recoveredLen
    fprintf('Warning: only got %d of %d bits\n', length(payloadBits), recoveredLen);
    recoveredLen = length(payloadBits);
end

finalBits = payloadBits(1:recoveredLen);

% inverse of tx, convert back to text
outputBits_reshape = reshape(finalBits, 8, [])';
charArray = char(outputBits_reshape + '0');
outputText = bin2dec(charArray).'; 
fprintf('Recovered: %s\n', char(outputText));

% print bit error
compareLen = min(length(bits), length(finalBits));
bitErrors = sum(bits(1:compareLen) ~= finalBits(1:compareLen));
actualBER = bitErrors / compareLen;
fprintf('Bit Error Rate (BER): %.4f\n', actualBER);
% else
%     disp('Error: Not enough bits received to match header length.');
% end

