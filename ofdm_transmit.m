clear all; clc;
%{

OFDM SIGNAL TRANSMITTER

MUST BE MATCHED WITH ofdm_receive.m

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

%centre around nfft/2 
activeCarriers = ((nfft/2) - (numActiveCarriers/2) : (nfft/2) + (numActiveCarriers/2)).';

activeCarriers(activeCarriers == nfft/2 + 1) = []; % remove dc 0 

% define idx for pilot, data and null carriers
pilotIdx = activeCarriers(1:pilotSpacing:end);
dataIdx = setdiff(activeCarriers, pilotIdx);
nullIdx = setdiff(1:nfft, activeCarriers).'; 

%% --== MESSAGE ==--

% msg = 'bass is really cool yo i swear yo';
file = fopen("./tests/textfile.txt","r")';
msg = fscanf(file, "%c");
fclose(file);

%% --== CONSTELLATION ==--

constSym = pskmod((0:M-1), M, pi/4); 
% define constellation diagram scope
cdScope = comm.ConstellationDiagram( ...
    'SamplesPerSymbol',1,...
    ShowReferenceConstellation=false,...
    ReferenceConstellation=constSym);


%% --== TX - QPSK ==--

% translate message into integers 0->M-1
binchars = dec2bin(msg, 8); 
bits = reshape(binchars.' - '0', [], 1); % convert from ascii nums to int
totalBits = numel(bits); 

% HEADER
headerBits = de2bi(totalBits, 16, 'left-msb').';

allBits = [ headerBits; bits ]; % stack header and payload

trellis = poly2trellis(3, [6 7]);
codedBits = convenc(allBits, trellis);


% FORMAT
bitsPerFrame = length(dataIdx) * k;
nFrames = ceil(numel(codedBits) / bitsPerFrame);
requiredTotalBits = nFrames * bitsPerFrame;

paddedBits = [codedBits; zeros(requiredTotalBits - numel(codedBits), 1)]; % pad to make square. 

bitgroups = reshape(paddedBits, k, [])'; % reshape by width k
inputSymbols = bi2de(bitgroups, 'left-msb'); % conv to int


% PREAMBLE
preamble = mod(0:numActiveCarriers-1, 4).'; % length of carriers, makes one symbol
preambleSignal = pskmod(preamble, M, pi/4);

% MODULATE
qpskSig = reshape(pskmod(inputSymbols, M, pi/4), length(dataIdx), nFrames);  
preambleData = preambleSignal(1:length(dataIdx)); % make the same length as data

% STACK
qpskSigFull = [preambleData, qpskSig];



%% --== TX - OFDM ==--

% DEFINE PILOTS
pilots = repmat(pskmod(0,M,pi/4),length(pilotIdx),nFrames+1);

tx_bb = ofdmmod(qpskSigFull, nfft, cplen, nullIdx, pilotIdx, pilots);


symbolLen = nfft + cplen;

% WINDOW 
roll = 16; % try 8–32 samples

win = ones(symbolLen,1);
rc = (0:roll-1)'/roll;
taper = 0.5*(1-cos(pi*rc)); % raised cosine

win(1:roll) = taper;
win(end-roll+1:end) = flipud(taper);

% apply per symbol
tx_bb_win = tx_bb;
numSyms = floor(length(tx_bb)/symbolLen);

for i = 1:numSyms
    idx = (i-1)*symbolLen + (1:symbolLen);
    tx_bb_win(idx) = tx_bb(idx).*win;
end

tx_bb = tx_bb_win;

%% --== TX - AUDIO PREP ==--

% add silence
silenceDuration = 2;
silenceSamples = round(silenceDuration * fs);
silence = zeros(silenceSamples, 1);

% upconvert to carrier
t = (0:length(tx_bb)-1)'/fs;

txPassband = real(tx_bb .* exp(1j*2*pi*fc*t));
txPassband = txPassband / max(abs(txPassband)) * 0.9;
txPassband = [ txPassband]; % not necessarily needed

% signal = awgn(txPassband, 50); % awgn if wanted
signal = txPassband;
audiowrite("ofdm_signal.wav", signal, fs);
disp("Tx done");

player = audioplayer(signal, fs);

% play(player);