function log_result(filename, params, ber, snr, throughput)
    if ~isfile(filename)
        fid = fopen(filename, 'w');
        fprintf(fid, 'test,ber,snr,throughput,M,cplen,pilotSpacing,numActiveCarriers,useCoding,fc,distance,volume\n');
        fclose(fid);
    end
    fid = fopen(filename, 'a');
    fprintf(fid, '%s,%.6f,%.4f,%.4f,%d,%d,%d,%d,%d,%d,%.1f,%d\n', ...
        params.test, ber, snr, throughput, ...
        params.M, params.cplen, params.pilotSpacing, ...
        params.numActiveCarriers, params.useCoding, ...
        params.fc, params.distance, params.volume);
    fclose(fid);
end