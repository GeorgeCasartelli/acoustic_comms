clear all 
clc

%{

Script for receiving and decoding from recorded audio file

Current state takes external audio source recorded elsewhere (e.g.
Audacity)

then uses a correlation (xcorr) value from a known preamble to align the
received signal

Signal then has header read for data payload length, and data payload is
extracted from 1rxSymbols(headerSymbols+1:headerSymbols+rxDataLength);

Phase offset is corrected by: 

phaseError = angle(sum(preambleSymbols .* conj(rxPreamble)));
rxSymbols = rxSymbols .* exp(1j * phaseError); ,


which rotates signal back to expected 4 quadrants


%}


fs = 48000;
fc = 1000;
Rs = 500; % symbol rate in baud

% symbol vars
M = 4;
k = log2(M);
phaseshift = pi/M;

% rcc filter params
rolloff = 0.25;
span = 12; % symbols
sps = fs / Rs; % = 480 samples per symbol
delay = span * sps / 2; % delay of rcc filter 

% define preamble and rccfilter
% preamble = randi([0 M-1], 32, 1);
preamble = [0;1;2;3;0;1;2;3;0;1;2;3;0;1;2;3;3;2;1;0;3;2;1;0;3;2;1;0;3;2;1;0];
rrcfilter = rcosdesign(rolloff, span, sps, 'sqrt');

preambleSymbols = pskmod(preamble,M,phaseshift); % modulate preamble into qpsk

preambleRef = upfirdn(preambleSymbols, rrcfilter, sps, 1); 
preambleRef = preambleRef(delay+1:end-delay); % 

%audio recorder

recorder = audiorecorder(fs,16,1);
recordDuration = 8;
recordblocking(recorder, recordDuration);
rx = getaudiodata(recorder);

% read audio file
% [rx, fs] = audioread("untitled.wav");
% rx = txNoise;

% figure(2);
% clf
% pwelch(rx, [], [], [], fs);
% title('rx psd');

t = (0:length(rx)-1).' / fs;
rx_bb = rx .* exp(-1j*2*pi*fc*t);

rxFiltered = filter(rrcfilter, 1, rx_bb); % matched filter
rxFiltered = rxFiltered(delay+1:end); % delay compensate to peak of rrc

[xc, lags] = xcorr(rxFiltered, preambleRef);
[~, idx] = max(abs(xc));

start = lags(idx);
rxAligned = rxFiltered(start:end);
% rxFiltered = rxFiltered(lags(idx)+length(preambleRef):end);

% figure(3); clf;
% plot(abs(xc));
% title('cross-correlation — preamble search');
% xlabel('lag'); ylabel('|xcorr|');

rxSymbols = rxAligned(1:sps:end);
rxPreamble = rxSymbols(1:length(preambleSymbols));

phaseError = angle(sum(preambleSymbols .* conj(rxPreamble)));
rxSymbols = rxSymbols .* exp(1j * phaseError);

rxSymbols = rxSymbols(length(preambleSymbols)+1:end);

headerSymbols = 8;
rxHeader = rxSymbols(1:headerSymbols);
% scatterplot(symbols);

% scatterplot(rxHeader);

headerData = pskdemod(rxHeader, M, phaseshift);
headerBitsRx = de2bi(headerData, k, 'left-msb');
headerBitsRx = reshape(headerBitsRx.', 1, []);
rxDataLength = bi2de(headerBitsRx, 'left-msb');
disp(rxDataLength);


rxPayload = rxSymbols(headerSymbols+1:headerSymbols+rxDataLength);
% scatterplot(rxPayload);
% expectedSymbols = length(dataIn);
dataOut = pskdemod(rxPayload, M, phaseshift);

if length(dataOut) < rxDataLength
    dataOut = [dataOut; zeros(rxDataLength - length(dataOut), 1)]; % pad output if too short
end

% len = min(length(dataOut), length(dataIn));

% BER = sum(dataOut(1:len) ~= dataIn(1:len)) / len;
% disp(BER)

output = de2bi(dataOut, 'left-msb');
% numChars = floor(size(output, 1) / 8);
% outputTrunc = output(1:numChars*8, :);
outputBits = reshape(output.', 8, [])';
charArray = char(outputBits + '0');
outputText = bin2dec(charArray).';
disp(char(outputText));
