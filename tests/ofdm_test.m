clear all; clc;

M = 4;
k = log2(M);
nfft = 2048;
cplen = 1024;
fs = 48000;
fc = 10000;

% define active carriers around the midpoint of nfft val
numActiveCarriers = 100;
pilotSpacing = 5;

% indices for first and second half of points
activeCarriers = ((nfft/2) - (numActiveCarriers/2) : (nfft/2) + (numActiveCarriers/2)).';

activeCarriers(activeCarriers == nfft/2 + 1) = []; % remove dc 0 
pilotIdx = activeCarriers(1:pilotSpacing:end);
dataIdx = setdiff(activeCarriers, pilotIdx);

nullIdx = setdiff(1:nfft, activeCarriers).'; % define null carrier idx


% msg = 'bass is really cool yo i swear yo';
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

% translate message into symbols
binchars = dec2bin(msg, 8); 
bits = reshape(binchars.' - '0', [], 1); 
totalBits = numel(bits); % will be included in packet header eventually

% header
headerBits = de2bi(totalBits, 16, 'left-msb').';

% allBits = [ headerBits; headerBits; bits ];
allBits = [ headerBits; bits ];

bitsPerFrame = length(dataIdx) * k;
nFrames = ceil(numel(allBits) / bitsPerFrame);
requiredTotalBits = nFrames * bitsPerFrame;

paddedBits = [allBits; zeros(requiredTotalBits - numel(allBits), 1)]; % pad to make square. 

bitgroups = reshape(paddedBits, k, [])'; % reshape by width k
inputSymbols = bi2de(bitgroups, 'left-msb'); % conv to int

% preamble = [0;1;2;3;0;1;2;3;0;1;2;3;0;1;2;3;3;2;1;0;3;2;1;0;3;2;1;0;3;2;1;0];
preamble = mod(0:numActiveCarriers-1, 4).';
preambleSignal = pskmod(preamble, M, pi/4);




% organise data into carriers x symbols matrix
qpskSig = reshape(pskmod(inputSymbols, M, pi/4), length(dataIdx), nFrames);  
preambleData = preambleSignal(1:length(dataIdx));

qpskSigFull = [preambleData, qpskSig];


pilots = repmat(pskmod(0,M,pi/4),length(pilotIdx),nFrames+1);

tx_bb = ofdmmod(qpskSigFull, nfft, cplen, nullIdx, pilotIdx, pilots);

symbolLen = nfft + cplen;
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

% add silence to start
silenceDuration = 2;
silenceSamples = round(silenceDuration * fs);
silence = zeros(silenceSamples, 1);

t = (0:length(tx_bb)-1)'/fs;

% upconvert to carrier
txPassband = real(tx_bb .* exp(1j*2*pi*fc*t));
txPassband = txPassband / max(abs(txPassband)) * 0.9;
txPassband = [silence; txPassband];

% sprinkle some wgn
% signal = awgn(txPassband, 50);
signal = txPassband;
audiowrite("ofdm_signal.wav", signal, fs);
disp("Tx done");

%% ---- RX ----
[rx, fs] = audioread("ofdm_signal.wav");

t_rx = (0:length(rx)-1)' /fs; % new time vector
rx_mixed = rx .* exp(-1j*2*pi*fc*t_rx); % mix down to baseband
rx_baseband = rx_mixed;

preamble_bb = ofdmmod(preambleData, nfft, cplen, nullIdx, pilotIdx, pilots(:,1)); % make a preamble reference

% rx_baseband = rx_mixed;

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

nSyms = size(pilots, 2);
nDataSyms = min(size(rxDataSyms, 2), nFrames);

disp(preamble);
disp(pskdemod(rxPreamble, M,pi/4));

% un rotate the signal using known preamble and received
H = rxPreamble ./ preambleData; 
x1_equalised = rxDataSyms ./ H;

rxData = pskdemod(x1_equalised, M, pi/4);

rxPilotSym = pskdemod(rxPilots,M,pi/4);
originalPilotSym = pskdemod(pilots, M, pi/4);

if isequal(rxPilotSym,originalPilotSym)
    disp("perfect pilots!");
else
    disp("uh oh bad pilots");
end

%% ---- DATA CHECK AND RECOVERY ----
isequal(rxData(:),inputSymbols)

cdScope(x1_equalised(:));
% cdScope(x1(:));
% scatterplot(x1);

drawnow;

%% ---- DATA RECOVERY WITH HEADER (RX) ----

% symbols convert bak to bits
rawBits = de2bi(rxData(:), k, 'left-msb');
allRxBits = reshape(rawBits.', [], 1);


if length(allRxBits) >= 16
    h1 = allRxBits(1:16);
    % h2 = allRxBits(17:32);
    % h3 = allRxBits(33:48);
    % headerVote = round((h1 + h2 + h3) / 3);
    headerVote = h1;
    recoveredLen = bi2de(headerVote', 'left-msb'); 
    fprintf('Header given message length: %d bits\n', recoveredLen);
    
    % grab message bits
    if length(allRxBits) >= (16 + recoveredLen)
        finalBits = allRxBits(17 : 17 + recoveredLen - 1);
        
        % conv back to text
        outputBits_reshape = reshape(finalBits, 8, [])';
        charArray = char(outputBits_reshape + '0');
        outputText = bin2dec(charArray).';
        fprintf('Recovered: %s\n', char(outputText));
        
        % calc ber
        bitErrors = sum(bits ~= finalBits);
        actualBER = bitErrors / totalBits;
        fprintf('Bit Error Rate (BER): %.4f\n', actualBER);
    else
        disp('Error: Length mismatch');
    end
else
    disp('Error: Couldnt decode header');
end


% output = de2bi(rxData(:), k, 'left-msb');
% outputBits = reshape(output.', [], 1);
% finalBits = outputBits(1:totalBits);
% 
% outputBits_reshape = reshape(finalBits, 8, [])';
% charArray = char(outputBits_reshape + '0');
% outputText = bin2dec(charArray).';
% 
% fprintf('Recovered: %s\n', char(outputText));


% figure('Name', 'Power Spectral Density');
% [pxx, f] = pwelch(txPassband, hanning(1024), 512, 1024, fs); 
% 
% plot(f/1000, 10*log10(pxx), 'LineWidth', 1.5, 'Color', [0 0.447 0.741]);
% grid on;
% xlabel('Frequency (kHz)');
% ylabel('Power/Frequency (dB/Hz)');
% title('PSD of the Transmitted OFDM Signal');
% 
% % Zoom in on the area of interest (fc +/- a few kHz)
% xlim([(fc-5000)/1000, (fc+5000)/1000]);
% 
% %% ---- VISUALIZING THE 2*fc COMPONENT ----
% % Calculate the spectrum of the mixed (but unfiltered) signal
% [pxx_mixed, f_mixed] = pwelch(rx_mixed, hanning(1024), 512, 1024, fs, 'centered');
% 
% figure('Name', 'Justification for Lowpass Filter');
% plot(f_mixed/1000, 10*log10(pxx_mixed), 'LineWidth', 1.5);
% grid on;
% hold on;
% 
% % Label the parts
% xlabel('Frequency (kHz)');
% ylabel('Magnitude (dB)');
% title('Spectrum of rx\_mixed (Before Filtering)');
% 
% % Draw a box or line around the parts to justify the filter
% xline(0, '--r', 'Baseband (Data)', 'LabelVerticalAlignment', 'bottom');
% xline(-2*fc/1000, '--k', 'Image at -2fc'); % If using complex mixer
% xline(2*fc/1000, '--k', 'Image at 2fc');
% 
% legend('Signal Spectrum', 'Target Data (0 Hz)', 'Unwanted Image');