clear; clc; close all

%% Parameters
M = 16;
nSym = 20e5;
sps = 4;
timingErr_ = [0 2];
SNRdB = 0:2:10;
flag1 = 0;

%% Sync schemes
TED_name = "Gardner (non-data-aided)";

%% RRC
txfilter = comm.RaisedCosineTransmitFilter(...
    "OutputSamplesPerSymbol", sps, "RolloffFactor", 0.35);
rxfilter = comm.RaisedCosineReceiveFilter(...
    "InputSamplesPerSymbol", sps, "DecimationFactor", 1.0, "RolloffFactor", 0.35);

BER = zeros(1, length(SNRdB));

for k = 1:length(timingErr_)
    timingErr = timingErr_(k);
    fprintf("\n------ %s -------\n", TED_name);

    for i = 1: length(SNRdB)
        snr = SNRdB(i);
        while true
            symbolSync = comm.SymbolSynchronizer(...
                "SamplesPerSymbol", sps,...
                "NormalizedLoopBandwidth", 0.01, ...
                "DampingFactor", 1.0, ...
                "TimingErrorDetector", TED_name);
            fprintf("Processing %ddB...\n", snr);
            data = randi([0 M-1], nSym, 1);
            modSig = dvbsapskmod(data, M, 's2');

            fixedDelay = dsp.Delay(timingErr);
            fixedDelaySym = ceil(fixedDelay.Length/sps);

            txSig = txfilter(modSig);
            delayedSig = fixedDelay(txSig);

            rxSig = awgn(delayedSig, snr, 'measured');

            rxSample = rxfilter(rxSig);

            rxSync = symbolSync(rxSample);

            recData = dvbsapskdemod(rxSync, M, 's2');

            sysDelay = dsp.Delay(fixedDelaySym + txfilter.FilterSpanInSymbols/2 + ...
                rxfilter.FilterSpanInSymbols/2);
            S = sysDelay(data);

            if (size(S, 1) == size(recData, 1))
                [~, ber] = biterr(S, recData);
                break;
                %         else
                %             flag1 = flag1+1;
                %             if (size(recData, 1) > size(S, 1))
                %                 rm_symbol = size(recData,1) - size(S, 1);
                %                 [~, ber] = biterr(S, recData(1:end-rm_symbol));
                %             else
                %                 add_symbol = size(S, 1) - size(recData, 1);
                %                 recData = [zeros(add_symbol, 1); recData];
                %                 [~, ber] = biterr(S, recData(1:end-rm_symbol));
                %             end
            end
        end
        BER(k, i) = ber;
    end
end

%% Plot
figure(1);
semilogy(SNRdB, BER(1,:), '-'); hold on
semilogy(SNRdB, BER(2,:), 'o');