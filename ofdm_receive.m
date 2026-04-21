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

%% --== DEFINE CARRIERS ==--

numActiveCarriers = 100;
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

% create pilot data ( symbol 00 ) based on nFrames
pilots = repmat(pskmod(0,M,pi/4),length(pilotIdx),nFrames+1);


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
rx_mixed = rx .* exp(-1j*2*pi*fc*t_rx); % mix down to baseband

% filter
lpCutoff = 2000; % may need adjusting
[b, a] = butter(6, lpCutoff/(fs/2));
rx_baseband = filter(b, a, rx_mixed);

% rx_baseband = rx_mixed;

preamble_bb = ofdmmod(preambleData, nfft, cplen, nullIdx, pilotIdx, pilots(:,1)); % make a preamble reference

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

disp(preamble);
disp(pskdemod(rxPreamble, M, pi/4));

nSyms = size(pilots, 2);
nDataSyms = min(size(rxDataSyms, 2), nFrames);
x1_equalised = zeros(length(dataIdx), nDataSyms);

for sym = 1:nDataSyms
    H_pilots = rxPilots(:, sym+1) ./ pilots(:, sym+1); % channel estimate at the pilot
    H_interp = interp1(pilotIdx, H_pilots, dataIdx, 'linear', 'extrap'); % interpolate est. to data subc.
    x1_equalised(:, sym) = rxDataSyms(:, sym) ./ H_interp; % apply channel correction
end

rxData = pskdemod(x1_equalised, M, pi/4);

%% ---- DATA CHECK AND RECOVERY ----
% isequal(rxData(:),inputSymbols)

cdScope(x1_equalised(:)); % plot on scope

% convert symbols back to bits
rawData = de2bi(rxData(:), k, 'left-msb');
allRxBits = reshape(rawData.', [], 1); % flatten to single bitstream

if length(allRxBits) >= 16

    recoveredLen = bi2de(allRxBits(1:16)', 'left-msb'); % get msg length
    fprintf('Header decoded! Message length: %d bits\n', recoveredLen);
    
    if length(allRxBits) >= (16 + recoveredLen)
        finalBits = allRxBits(17 : 17 + recoveredLen - 1); %skip header, grab data
        
        % inverse of tx, convert back to text
        outputBits_reshape = reshape(finalBits, 8, [])';
        charArray = char(outputBits_reshape + '0');
        outputText = bin2dec(charArray).';
        fprintf('Recovered: %s\n', char(outputText));
        
        % print bit error
        bitErrors = sum(bits ~= finalBits);
        actualBER = bitErrors / recoveredLen;
        fprintf('Bit Error Rate (BER): %.4f\n', actualBER);
    else
        disp('Error: Not enough bits received to match header length.');
    end
else
    disp('Error: Could not even decode the 16-bit header.');
end
