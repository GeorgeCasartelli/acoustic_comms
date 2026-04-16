clear all; clc;

M = 4;
k = log2(M);
nfft = 1024;
cplen = 32;
fs = 48000;
fc = 6000;

%% ---- CARRIER DEFINITIONS ----

numActiveCarriers = 42;
pilotSpacing = 10;

% indices for first and second half of points
activeCarriers = ((nfft/2) - (numActiveCarriers/2) : (nfft/2) + (numActiveCarriers/2)).';

activeCarriers(activeCarriers == nfft/2 + 1) = []; % remove dc 0 
pilotIdx = activeCarriers(1:pilotSpacing:end);
dataIdx = setdiff(activeCarriers, pilotIdx);

nullIdx = setdiff(1:nfft, activeCarriers).'; % define null carrier idx



%% ---- MESSAGE ----

file = fopen("textfile.txt","r")';
msg = fscanf(file, "%c");
fclose(file);

binchars = dec2bin(msg, 8); 
bits = reshape(binchars.' - '0', [], 1); 
totalBits = numel(bits); % will be included in packet header eventually

bitsPerFrame = length(dataIdx) * k;
nFrames = ceil(numel(bits) / bitsPerFrame);
paddedBits = [bits; zeros(nFrames*bitsPerFrame - totalBits, 1)]; % pad to make square. 
bitgroups = reshape(paddedBits, k, [])'; % reshape by width k
inputSymbols = bi2de(bitgroups, 'left-msb'); % conv to int


pilots = repmat(pskmod(0,M,pi/4),length(pilotIdx),nFrames+1);

%% ---- CONTSTELLATION -----
constSym = pskmod((0:M-1), M, pi/4); 
% define constellation diagram scope
cdScope = comm.ConstellationDiagram( ...
    'SamplesPerSymbol',1,...
    ShowReferenceConstellation=false,...
    ReferenceConstellation=constSym);

%% ---- RECORDER -----
recorder = audiorecorder(fs,16,1);
recordDuration = 5;
recordblocking(recorder, recordDuration);
rx = getaudiodata(recorder);

% [rx, fs] = audioread("ofdm_signal.wav");


%% ---- RX -----
preamble = mod(0:numActiveCarriers-1, 4).';
preambleSignal = pskmod(preamble, M, pi/4);
preambleData = preambleSignal(1:length(dataIdx));


t_rx = (0:length(rx)-1)' /fs; % new time vector
rx_mixed = rx .* exp(-1j*2*pi*fc*t_rx); % mix down to baseband

lpCutoff = 2000; % Adjust based on your numActiveCarriers width
[b, a] = butter(6, lpCutoff/(fs/2));
rx_baseband = filter(b, a, rx_mixed);

% rx_baseband = rx_mixed;

preamble_bb = ofdmmod(preambleData, nfft, cplen, nullIdx, pilotIdx, pilots(:,1)); % make a preamble reference

% rx_baseband = rx_mixed;

%% ---- ALIGN ----
[xc, lags] = xcorr(rx_baseband, preamble_bb);
[~, idx] = max(abs(xc));
start = lags(idx);

rxAligned = rx_baseband(start:end); % align based on preamble reference and rx

% trim down signal to be multiple of symbollen for ofdmdemod
symbolLen = nfft + cplen;
numSymbolsReceived = floor(length(rxAligned) / symbolLen);
rxTrimmed = rxAligned(1 : numSymbolsReceived * symbolLen);
% 
% x1 = ofdmdemod(rxTrimmed, nfft, cplen, 0, nullIdx);
[x1, rxPilots] = ofdmdemod(rxTrimmed, nfft, cplen, 0, nullIdx, pilotIdx);
 
rxPreamble = x1(:, 1); % get first symbol (preamble) 
rxDataSyms = x1(:, 2:end); % get payload

disp(preamble);
disp(pskdemod(rxPreamble, M,pi/4));

% un rotate the signal using known preamble and received
H = rxPreamble ./ preambleData; 
x1_equalised = rxDataSyms ./ H;

rxData = pskdemod(x1_equalised, M, pi/4);

%% ---- DATA CHECK AND RECOVERY ----
isequal(rxData(:),inputSymbols)

cdScope(x1_equalised(:));
% cdScope(x1(:));
% scatterplot(x1);

drawnow;

output = de2bi(rxData(:), k, 'left-msb');
outputBits = reshape(output.', [], 1);
finalBits = outputBits(1:totalBits);

outputBits_reshape = reshape(finalBits, 8, [])';
charArray = char(outputBits_reshape + '0');
outputText = bin2dec(charArray).';

fprintf('Recovered: %s\n', char(outputText));