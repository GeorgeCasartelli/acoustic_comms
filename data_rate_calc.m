function [R_b] = data_rate_calc(fs, nfft, cplen, dataCarriers, R, bpsc)
%DATA_RATE_CALC Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    fs
    nfft
    cplen
    dataCarriers
    R
    bpsc
end

arguments (Output)
    R_b
end

symbolSize = nfft + cplen;
symbolTime = symbolSize / fs;
rawSymbolRate = (dataCarriers * bpsc) / symbolTime;
R_b = rawSymbolRate * R;

end