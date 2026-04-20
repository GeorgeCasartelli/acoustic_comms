fs = 48000;
duration = 10;
f1 = 20; f2 = 20000;
t = (0:1/fs:duration)';

sweep = sin(2*pi* duration/log(f2/f1) * exp(t/duration * log(f2/f1)) - 1);
sweep = sweep / max(abs(sweep)) * 0.9;
audiowrite('sweep.wav', sweep, fs);

% player = audioplayer(sweep, fs);
recorder = audiorecorder(fs, 16, 1);

record(recorder);
% play(player);
pause(duration + 2);
stop(recorder);

recorded = getaudiodata(recorder);

N = min(length(sweep), length(recorded));
H = fft(recorded(1:N)) ./ fft(sweep(1:N));
freqs = (0:N-1) * fs/N;

H_mag = 20*log10(abs(H(1:N/2)));
freqs_plot = freqs(1:N/2)/1000;

smoothed = smoothdata(H_mag, 'gaussian', 200);

plot(freqs_plot, smoothed)
xlabel('Frequency (kHz)')
ylabel('Magnitude (dB)')
title('Room Response (smoothed)')
xlim([0 20])
ylim([-40 40])
grid on
xline(6, 'r--', '6kHz');
xline(9, 'r--', '9kHz');
xline(12, 'r--', '12kHz');
xline(14, 'r--', '14kHz');