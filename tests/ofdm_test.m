% clear all; clc;

M = 4;
k = log2(M);
nfft = 1024;
cplen = 32;
fs = 48000;
fc = 9000;

% define active carriers around the midpoint of nfft val
numDataCarriers = 20;
firstHalf = (nfft/2 - 10) : (nfft/2 - 1);
secondHalf = (nfft/2 + 1) : (nfft/2 + 10);
activeCarriers = [firstHalf, secondHalf];
nullIdx = setdiff(1:nfft, activeCarriers).';


% msg = 'bonjourno my name is giorgio. ';
file = fopen("textfile.txt","r")';
msg = fscanf(file, "%c");
fclose(file);


constSym = pskmod((0:M-1), M, pi/4); 
% define constellation diagram scope
cdScope = comm.ConstellationDiagram( ...
    'SamplesPerSymbol',1,...
    ShowReferenceConstellation=false,...
    ReferenceConstellation=constSym);

%% ---- TX ----
binchars = dec2bin(msg, 8); % get kx8 matrix
bits = reshape(binchars.' - '0', [], 1); % stack vertical
totalBits = numel(bits); % will be included in packet header eventually

bitsPerFrame = numDataCarriers * k;
nFrames = ceil(numel(bits) / bitsPerFrame);
paddedBits = [bits; zeros(nFrames*bitsPerFrame - totalBits, 1)]; % pad to make square. 
bitgroups = reshape(paddedBits, k, [])'; % reshape by width k
inputSymbols = bi2de(bitgroups, 'left-msb'); % conv to int

% organise data into carriers x symbols matrix
qpskSig = reshape(pskmod(inputSymbols, M, pi/4), numDataCarriers, nFrames);  

tx_bb = ofdmmod(qpskSig, nfft, cplen, nullIdx);

% up convert to passband

% preamble = ones(samplesPerSymbol, 1);
t = (0:length(tx_bb)-1)'/fs;
txPassband = real(tx_bb .* exp(1j*2*pi*fc*t));
txPassband = txPassband / max(abs(txPassband)) * 0.9; % normalise

signal = awgn(txPassband, 100);
audiowrite("ofdm_signal.wav", signal, fs);
disp("Tx done");

%% ---- RX ----
[rx, fs] = audioread("ofdm_signal.wav");

rx_mixed = 2 * rx .* exp(-1j*2*pi*fc*t);

% filter out image freqs
bw = (numDataCarriers/nfft) * fs; 
rx_baseband = lowpass(rx_mixed, bw/2, fs);

x1 = ofdmdemod(rx_baseband, nfft, cplen, 0, nullIdx);

x1_normalised = x1 / max(abs(x1(:))); % was at 40, not at 1
rxData = pskdemod(x1, M, pi/4);

%% ---- DATA CHECK AND RECOVERY ----
isequal(rxData(:),inputSymbols)

cdScope(x1_normalised(:));

drawnow;

output = de2bi(rxData(:), k, 'left-msb');
outputBits = reshape(output.', [], 1);
finalBits = outputBits(1:totalBits);

outputBits_reshape = reshape(finalBits, 8, [])';
charArray = char(outputBits_reshape + '0');
outputText = bin2dec(charArray).';

fprintf('Recovered: %s\n', char(outputText));

