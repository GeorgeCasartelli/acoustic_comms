clear all; clc; close all;

%{

OFDM SIGNAL RECEIVER

MUST BE MATCHED WITH ofdm_transmit.m

%}

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end

%% --== DEFINE SCRIPT PARAMS ==--

useCoding = true;
txMode = 'text';
imageSize = 128;
headerSize = 32;
useInterleaving = true;
modScheme = "qpsk";
constraintLength = 7;

phyParams.test = 'blocking';
phyParams.trial = 1;

phyParams.distance = 3.0;
phyParams.volume = 100;

phyParams.txAngle = 0;
phyParams.rxAngle = 0;

phyParams.blocked = 0;
phyParams.blockerType = "none";

phyParams.reflectionMode = "wall";

phyParams.modscheme = modScheme;
phyParams.useCoding = useCoding;

phyParams.notes = "wall_reflection_rx_infront";

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
recordDuration = ternary(strcmp(txMode, 'image'), 30, 6);
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

nDataSyms = size(rxDataSyms, 2);
pilotSym = pskmod(0,M,pi/4);
pilots = repmat(pilotSym, length(pilotIdx), nDataSyms+1);
nSyms = size(pilots, 2);


%% --== EQUALISE FROM PILOTS ==-- 
x1_equalised = zeros(length(dataIdx), nDataSyms);
for sym = 1:nDataSyms
    H_pilots = rxPilots(:, sym+1) ./ pilots(:, sym+1); % channel estimate at the pilot
    H_interp = interp1(pilotIdx, H_pilots, dataIdx, 'linear', 'extrap'); % interpolate est. to data subc.
    x1_equalised(:, sym) = rxDataSyms(:, sym) ./ H_interp; % apply channel correction
end


rxData = ternary(strcmp(modScheme, "qpsk"), pskdemod(x1_equalised, M, pi/4), qamdemod(x1_equalised, M, "UnitAveragePower", true));
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
recoveredLen = 28672;

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

%% --== SNR ESTIMATION ==--
% calc SNR based on the deviation of received pilots from the ideal pilot symbol
errors  = rxPilots - pilotSym;
P_sig   = mean(abs(rxPilots(:)).^2);
P_noise = mean(abs(errors(:)).^2);

snr_est = 10 * log10(P_sig / P_noise);
if P_noise == 0; P_noise = 1e-10; end


%% --== RECONSTRUCT IDEAL SYMBOLS  ==--

% rebuild tx side symbol indices for SER comparison
% trellis = poly2trellis(7, [171 133]);
% interleaveSeed = 12345;

% encode known symbols
if useCoding
    tx_header_re = convenc(de2bi(totalBits, 32, 'left-msb').', trellis);
    payload_coded_re = convenc(bits, trellis); % 'bits' are your final decoded bits
else
    tx_header_re = de2bi(totalBits, 32, 'left-msb').';
    payload_coded_re = bits;
end

% interleave known symbols
if useInterleaving
    payload_int_re = randintrlv(payload_coded_re, interleaveSeed);
else
    payload_int_re = payload_coded_re;
end

tx_bits_re = [tx_header_re; payload_int_re];
bitsPerFrame = length(dataIdx) * k;
nFramesData = ceil(numel(tx_bits_re) / bitsPerFrame);
paddedBits_re = [tx_bits_re; zeros(nFramesData * bitsPerFrame - numel(tx_bits_re), 1)];

bitgroups_re = reshape(paddedBits_re, k, [])';
inputSyms_re = bi2de(bitgroups_re, 'left-msb');
idealDataSymbols = reshape(inputSyms_re, length(dataIdx), nFramesData);

%% --== ALIGN AND COMPARE ==--

compareLen = min(size(rxData,2), nFramesData);

rxTest = rxData(:, 1:compareLen);
txTest = idealDataSymbols(:, 1:compareLen);

% symbol error rate per subcarrier
serPerSubcarrier = sum(rxTest ~= txTest, 2) / compareLen;

% overall SER
totalSymbolErrors = sum(rxTest(:) ~= txTest(:));
totalSymbols = numel(rxTest);

actualSER = totalSymbolErrors / totalSymbols;

fprintf('Symbol Error Rate (SER): %.6f\n', actualSER);

%% --== DASHBOARD ==--

display_ofdm_dashboard(rx, fs, start, recoveredLen, headerSize, dataIdx, k, ...
    fc, numActiveCarriers, H_interp, rxPilots, pilotSym, pilotIdx, nfft, cplen, ...
    useCoding, useInterleaving, txMode, actualBER, R_b, modScheme, snr_est, ...
    constraintLength, serPerSubcarrier);


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

log_physical_results( ...
    'physical_results.csv', ...
    phyParams, ...
    actualBER, ...
    actualSER, ...
    snr_est, ...
    R_b);