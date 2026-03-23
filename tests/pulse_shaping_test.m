clear all 
clc

fs = 48000;
fc = 800;
N = 1000;
Rs = 100; % symbol rate in baud

% symbol vars
M = 4;
k = log2(M);

% rcc filter params
rolloff = 0.25;
span = 12; % symbols
sps = fs / Rs; % = 480 samples per symbol
delay = span * sps / 2; % delay of rcc filter 

% define preamble and rccfilter
preamble = randi([0 M-1], 32, 1);
rrcfilter = rcosdesign(rolloff, span, sps, 'sqrt');

% generate message and encode to binary then to int [0 M-1]
msg = 'I love my girlfriend! She is so sexy!';
binchars = dec2bin(msg, 8); % get kx8 matrix
bits = reshape(binchars.' - '0', [], 1); % stack vertical
bitgroups = reshape(bits, k, [])'; % reshape by width k

inputSymbols = bi2de(bitgroups, 'left-msb'); % conv to int

% generate header
headerLength = 16;
dataLength = length(bitgroups);
dataLengthBits = dec2bin(dataLength);
if length(dataLengthBits) < headerLength
    dataLengthBits = [repmat('0', 1, headerLength - length(dataLengthBits)), dataLengthBits];

end
disp(dataLengthBits);
% headerBitChars = dec2bin(dataLengthBits);
headerBits = reshape(dataLengthBits.' - '0', [], 1);
headerBitGroups = reshape(headerBits, k, [])';
header = bi2de(headerBitGroups, 'left-msb');

% plot rrc
[H, f] = freqz(rrcfilter, 1, 4096, fs);

% figure;
% plot(f, 20*log10(abs(H)));
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% title('RRC Filter Frequency Response');
% grid on;
% xlim([0 Rs]);  % zoom to one symbol rate worth


phaseshift = pi/M;
% dataIn = randi([0 M-1], N, 1);
dataIn=inputSymbols;

symbols = pskmod(inputSymbols, M, phaseshift); % modulate data into qpsk symbs
preambleSymbols = pskmod(preamble,M,phaseshift); % modulate preamble into qpsk
headerSymbols = pskmod(header, M, phaseshift);

txSymbols = [preambleSymbols; headerSymbols; symbols]; % package preamble w/ symbols
% txSymbols = [preambleSymbols; symbols];
% txSymbols = symbols;

% upsample preamble and delay compensate for rrc filter
preambleRef = upfirdn(preambleSymbols, rrcfilter, sps, 1); 
preambleRef = preambleRef(delay+1:end-delay); % 

txShaped = upfirdn(txSymbols, rrcfilter, sps, 1); % upsample & shaping

t = (0:length(txShaped)-1).' / fs;

% carrier modulation
txPassband = real(txShaped .* exp(1j*2*pi*fc*t));
txPassband = txPassband / max(abs(txPassband));

silenceDuration = 2;
silenceSamples = round(silenceDuration * fs);
silence = zeros(silenceSamples, 1);

txPassband = [silence; txPassband];

snr_range = [-10:5:20];
ber_measured = zeros(size(snr_range));


% disp(snr_range(i))

% txNoise = awgn(txPassband, 5);
txNoise = txPassband;

audiowrite('bpsk_tx.wav', txNoise, fs);
disp("Tx done");

% scatterplot(symbols);
% 
% figure(2); clf%plot waveform
% plot(t(1:sps*20), txPassband(1:sps*20));
% xlabel("time");
% ylabel("amplitide")
% title("First 10 symbols")

% figure(2); % spectro
% clf
% pwelch(txPassband, [], [], [], fs);
% title('psd')


%receive

[rx, fs] = audioread("bpsk_tx.wav");
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

timingOffset = 0;
rxSymbols = rxAligned(1:sps:end);
rxSymbols = rxSymbols(length(preambleSymbols)+1:end);
rxHeader = rxSymbols(1:8);
rxPayload = rxSymbols(9:end);
% scatterplot(symbols);
% scatterplot(rxSymbols);
% scatterplot(rxHeader);

headerData = pskdemod(rxHeader, M, phaseshift);
headerBitsRx = de2bi(headerData, k, 'left-msb');
headerBitsRx = reshape(headerBitsRx.', 1, []);
rxDataLength = bi2de(headerBitsRx, 'left-msb');
disp(rxDataLength);


expectedSymbols = length(dataIn);
dataOut = pskdemod(rxPayload, M, phaseshift);

if length(dataOut) < rxDataLength
    dataOut = [dataOut; zeros(rxDataLength - length(dataOut), 1)]; % pad output if too short
end

len = min(length(dataOut), length(dataIn));

BER = sum(dataOut(1:len) ~= dataIn(1:len)) / len;
disp(BER)

output = de2bi(dataOut, 'left-msb');
% numChars = floor(size(output, 1) / 8);
% outputTrunc = output(1:numChars*8, :);
outputBits = reshape(output.', 8, [])';
charArray = char(outputBits + '0');
outputText = bin2dec(charArray).';
disp(char(outputText));
% expectedStart = silenceSamples; % roughly where we expect it
% fprintf('Expected start: ~%d samples\n', silenceSamples + delay);
% fprintf('Detected start: %d samples\n', start);
% fprintf('Error: %d samples\n', abs(start - (silenceSamples - delay)));
