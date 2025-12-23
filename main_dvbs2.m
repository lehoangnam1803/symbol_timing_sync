clear, clc, close all

%% Parameters
L        = 32;         % Oversampling factor
M        = 16;          % Constellation order

nSymbols = 1e5;        % Number of transmit symbols
Bn_Ts    = 0.01;       % Loop noise bandwidth (Bn) times symbol period (Ts)
eta      = 1;          % Loop damping factor

rollOff  = 0.25;        % Pulse shaping roll-off factor
timeOffset = 20;       % Simulated channel delay in samples
fsOffsetPpm = 0;       % Sampling clock frequency offset in ppm
rcDelay  = 10;         % Raised cosine (combined Tx/Rx) delay

EsN0     = 10:2:16;         % Target Es/N0
% snr      = EsN0 - 10*log10(L);
Ex       = 1;          % Average symbol energy

TED      = 'MLTED';    % TED (MLTED, ELTED, ZCTED, GTED, or MMTED)
intpl    = 3;          % 0) Polyphase; 1) Linear; 2) Quadratic; 3) Cubic
forceZc  = 0;          % Use to force zero-crossings and debug self-noise

falsecnt = 0;
Nmc = 1e4;
FSP = zeros(1, length(EsN0));

%% System Objects
% Tx Filter
TXFILT = comm.RaisedCosineTransmitFilter( ...
    'OutputSamplesPerSymbol', L, ...
    'RolloffFactor', rollOff, ...
    'FilterSpanInSymbols', rcDelay);

% Rx Matched Filter (MF)
RXFILT = comm.RaisedCosineReceiveFilter( ...
    'InputSamplesPerSymbol', L, ...
    'DecimationFactor', 1, ...
    'RolloffFactor', rollOff, ...
    'FilterSpanInSymbols', rcDelay);
mf = RXFILT.coeffs.Numerator; % same as "rcosdesign(rollOff, rcDelay, L)"

% Digital Delay
DELAY = dsp.Delay(timeOffset);

%% Loop Constants
% Time-error Detector Gain (TED Gain)
Kp = calcTedKp(TED, rollOff);

% Scale Kp based on the average symbol energy (at the receiver)
K  = 1; % Assume channel gain is unitary (or that an AGC is used)
Kp = K * Ex * Kp;
% NOTE: if using the GTED when K is not 1, note the scaling is K^2 not K.

% Counter Gain
K0 = -1;
% Note: this is analogous to VCO or DDS gain, but in the context of timing
% recovery loop.

% PI Controller Gains:
[ K1, K2 ] = piLoopConstants(Kp, K0, eta, Bn_Ts, L);

% fprintf("Loop constants:\n");
% fprintf("K1 = %g; K2 = %g; Kp = %g\n", K1, K2, Kp);

% MATLAB's symbol synchronizer for comparison
tedMap = containers.Map({'ELTED', 'ZCTED', 'GTED', 'MMTED'}, ...
    {'Early-Late (non-data-aided)', ...
    'Zero-Crossing (decision-directed)', ...
    'Gardner (non-data-aided)', ...
    'Mueller-Muller (decision-directed)'
    });


for ii = 1:length(EsN0)
    snr = EsN0(ii) - 10*log10(L);
    for i = 1:Nmc
        % fprintf("Processing i = %d ...\n", i);
        %% Simulation: Tx -> Channel -> Rx Matched Filtering -> Symbol Synchronizer
        % Tx Filter
        % txSig = step(TXFILT, modSig);
        [txSig, data] = dvbs2_waveform(L, rollOff);

        % Sampling clock frequency offset
        fsRatio = 1 + (fsOffsetPpm * 1e-6); % Rx/Tx clock frequency ratio
        tol = 1e-9;
        [P, Q] = rat(fsRatio, tol); % express the ratio as a fraction P/Q
        txResamp = resample(txSig, P, Q);

        % Channel
        delaySig = step(DELAY, txResamp);
        % delaySig = [zeros(timeOffset, 1); txResamp];
        txSigPower = 1 / sqrt(L);
        [rxSeq, noiseVar] = awgn(delaySig, snr, 'measured');

        % Rx matched filter (MF)
        mfOut = step(RXFILT, rxSeq);
        % mfOut = mfOut(rcDelay + 1:end);

        %% Symbol Timing Recovery
        % Downsampled symbols without symbol timing recovery
        rxNoSync = downsample(mfOut, L);

        % Downsampled symbols with perfect symbol timing recovery
        rxPerfectSync = downsample(mfOut, L, timeOffset);

        % Our symbol timing recovery implementation
        [ rxSync1, mu ] = symbolTimingSync(TED, intpl, L, rxSeq, mfOut, K1, K2, ...
            rollOff, rcDelay);

        mu_ref = mod(timeOffset / L, 1);
        mu_ss = mean(mu(end-1001:end));
        isFalse = abs(mu_ss - mu_ref);
        if isFalse > 0.25
            falsecnt = falsecnt + 1;
        end
    end
    FSP(1, ii) = falsecnt / Nmc;
    % disp(FSP)
end
%% ===================== BER / SER COMPUTATION =====================

% grpDelaySym = rcDelay;
% rxSync1 = rxSync1(grpDelaySym+1:end);
% rxPerfectSync = rxPerfectSync(grpDelaySym+1:end);
% [bits,FramesLost,pktCRCStat] = dvbs2BitRecover(rxSync1,10^(-EsN0/10));

% [~, ber] = biterr(bits{1,1}, data);
%
% fprintf("\nFrames Lost: %d", FramesLost);
% fprintf("\nBER = %d\n", ber);
fprintf("False synchronization prob.: %.4f\n", FSP);
fprintf("SNR (dB) = %.4f\n", snr);