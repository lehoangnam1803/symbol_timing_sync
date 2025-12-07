clear;clc;close all

M = 2;         % Modulation order for BPSK
nSym = 5e6;  % Number of symbols in a packet
sps = 4;       % Samples per symbol
timingErr = 2; % Samples of timing error
SNRdB = 0:2:10;      % Signal-to-noise ratio (dB)

TED_name = [{'Early-Late (non-data-aided)'}, {'Gardner (non-data-aided)'},...
    {'Mueller-Muller (decision-directed)'}, {'Zero-Crossing (decision-directed)'}];

BER = zeros(4, length(SNRdB));

for ii = 1:4



    for i = 1:length(SNRdB)
        snr = SNRdB(i);
        flag1 = 0;
        flag2 = 0;
        while true
            txfilter = comm.RaisedCosineTransmitFilter(...
                'OutputSamplesPerSymbol',sps);
            rxfilter = comm.RaisedCosineReceiveFilter(...
                'InputSamplesPerSymbol',sps,'DecimationFactor',1);
            symbolSync = comm.SymbolSynchronizer(...
                'SamplesPerSymbol',sps, ...
                'NormalizedLoopBandwidth',0.01, ...
                'DampingFactor',1.0, ...
                'TimingErrorDetector',TED_name{1});

            data = randi([0 M-1],nSym,1);
            modSig = pskmod(data,M);
            fixedDelay = dsp.Delay(timingErr);
            fixedDelaySym = ceil(fixedDelay.Length/sps); % Round fixed delay to nearest integer in symbols
            txSig = txfilter(modSig);
            delayedSig = fixedDelay(txSig);
            rxSig = awgn(delayedSig,snr,'measured');
            rxSample = rxfilter(rxSig);
            %     scatterplot(rxSample(10000:end),2)
            rxSync = symbolSync(rxSample);
            %     scatterplot(rxSync(10000:end),2)
            recData = pskdemod(rxSync,M);
            sysDelay = dsp.Delay(fixedDelaySym + txfilter.FilterSpanInSymbols/2 + ...
                rxfilter.FilterSpanInSymbols/2);
            if (size(sysDelay(data),1) == size(recData,1))
                [numErr1,ber1] = biterr(sysDelay(data),recData);
                if (ber1 < 0.4)
                    break;
                else
                    flag2 = 1;
                end
            else
                flag1 = 1;
            end
        end
        BER(ii,i) = ber1;
    end
end

figure;
for i = 1:4
    semilogy(SNRdB, BER(i,:), '-o'); hold on
    legend(TED_name{i})
end
legend(TED_name{1}, TED_name{2}, TED_name{3}, TED_name{4})
%semilogy(SNRdB, BER(1,:), '-o');
grid on

