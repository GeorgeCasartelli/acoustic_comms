function [R_b, efficiency] = data_rate_calc(fs, nfft, cplen, dataCarriers, activeCarriers, R, bpsc)
%DATA_RATE_CALC Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    fs
    nfft
    cplen
    dataCarriers
    activeCarriers
    R
    bpsc
end

arguments (Output)
    R_b
    efficiency
end

symbolSize = nfft + cplen;
symbolTime = symbolSize / fs;
rawSymbolRate = (dataCarriers * bpsc) / symbolTime;
R_b = rawSymbolRate * R;

cpEfficiency = nfft / (nfft + cplen);
pilotEfficiency = dataCarriers / activeCarriers;
codingEfficiency = R;
efficiency = cpEfficiency * pilotEfficiency* codingEfficiency;

fprintf("---------DATA RATE ANALYSIS------------\n");
fprintf("Symbol duration:   %.2f ms\n", symbolTime * 1000);
fprintf("Raw data rate:     %.2f bps\n", R_b);
fprintf("CP Efficiency:     %.2f%%\n", cpEfficiency * 100);
fprintf("Pilot Efficiency:  %.2f%%\n", pilotEfficiency * 100);
fprintf("Coding Efficiency: %.2f%%\n", codingEfficiency * 100);
fprintf("Efficiency:        %.2f%%\n", efficiency * 100);
fprintf("Throughput:        %.2f bps\n", R_b * pilotEfficiency * cpEfficiency);
fprintf("---------------------------------------\n");
end