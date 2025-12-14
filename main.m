clear, clc, close all

%% Debug Configuration

debug_tl_static  = 0; % Show static debug plots after sync processing
debug_tl_runtime = 0; % Open scope for debugging of sync loop iterations

%% Parameters
L        = 32;         % Oversampling factor
M        = 16;          % Constellation order
N        = 2;          % Dimensions per symbol (1 for PAM, 2 for QAM)
nSymbols = 1e5;        % Number of transmit symbols
Bn_Ts    = 0.01;       % Loop noise bandwidth (Bn) times symbol period (Ts)
eta      = 1;          % Loop damping factor
rollOff  = 0.2;        % Pulse shaping roll-off factor
timeOffset = 25;       % Simulated channel delay in samples
fsOffsetPpm = 0;       % Sampling clock frequency offset in ppm
rcDelay  = 10;         % Raised cosine (combined Tx/Rx) delay
EsN0     = 10;         % Target Es/N0
Ex       = 1;          % Average symbol energy
TED      = 'GTED';    % TED (MLTED, ELTED, ZCTED, GTED, or MMTED)
intpl    = 3;          % 0) Polyphase; 1) Linear; 2) Quadratic; 3) Cubic
forceZc  = 0;          % Use to force zero-crossings and debug self-noise

%% System Objects
% Tx Filter
TXFILT = comm.RaisedCosineTransmitFilter( ...
    'OutputSamplesPerSymbol', L, ...
    'RolloffFactor', rollOff, ...
    'FilterSpanInSymbols', rcDelay);

% Rx Matched Filter (MF)
%
% NOTE: in most simulations, the decimation factor would be L below.
% However, here we want to process the filtered sequence as-is, without the
% downsampling stage. The symbol timing recovery loop processes the
% fractionally-spaced sequence and handles the downsampling process.
RXFILT = comm.RaisedCosineReceiveFilter( ...
    'InputSamplesPerSymbol', L, ...
    'DecimationFactor', 1, ...
    'RolloffFactor', rollOff, ...
    'FilterSpanInSymbols', rcDelay);
mf = RXFILT.coeffs.Numerator; % same as "rcosdesign(rollOff, rcDelay, L)"

% Digital Delay
DELAY = dsp.Delay(timeOffset);

% Reference constellation for MER measurement
if (N==2)
    const = qammod(0:M-1,M);
else
    const = pammod(0:M-1,M);
end
Ksym = modnorm(const, 'avpow', Ex);
const = Ksym * const;

% MER Meter
mer = comm.MER;
mer.ReferenceSignalSource = 'Estimated from reference constellation';
mer.ReferenceConstellation = const;

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

fprintf("Loop constants:\n");
fprintf("K1 = %g; K2 = %g; Kp = %g\n", K1, K2, Kp);

% MATLAB's symbol synchronizer for comparison
tedMap = containers.Map({'ELTED', 'ZCTED', 'GTED', 'MMTED'}, ...
    {'Early-Late (non-data-aided)', ...
     'Zero-Crossing (decision-directed)', ...
     'Gardner (non-data-aided)', ...
     'Mueller-Muller (decision-directed)'
     });
if strcmp(TED, "MLTED")
    warning("MLTED not supported by MATLAB's synchronizer - using ZCTED");
    matlabTed = "ZCTED";
else
    matlabTed = TED;
end
SYMSYNC = comm.SymbolSynchronizer(...
    'TimingErrorDetector', tedMap(matlabTed), ...
    'SamplesPerSymbol', L, ...
    'NormalizedLoopBandwidth', Bn_Ts, ...
    'DampingFactor', eta, ...
    'DetectorGain', Kp);

%% Random Transmit Symbols
if (forceZc)
    % Force zero-crossings within the transmit symbols. Use to eliminate
    % the problem of self-noise and debug the operation of the loop
    data = zeros(nSymbols, 1);
    data(1:2:end) = M-1;
else
    data = randi([0 M-1], nSymbols, 1);
end

if (N==2)
    modSig = Ksym * qammod(data, M);
else
    modSig = real(Ksym * pammod(data, M));
end
% Ensure the average symbol energy is unitary, otherwise the loop constants
% must be altered (because Kp, the TED gain, must scale accordingly).

%% Simulation: Tx -> Channel -> Rx Matched Filtering -> Symbol Synchronizer
% Tx Filter
txSig = step(TXFILT, modSig);

% Sampling clock frequency offset
%
% The frequencies produced by the Tx and Rx sampling clocks are often
% significantly distinct, unless both sides adopt high-accuracy oscillators
% (e.g., atomic clocks) or clock disciplining mechanisms, such as with
% GPSDOs. Simulate this relative frequency offset by resampling the signal.
fsRatio = 1 + (fsOffsetPpm * 1e-6); % Rx/Tx clock frequency ratio
tol = 1e-9;
[P, Q] = rat(fsRatio, tol); % express the ratio as a fraction P/Q
txResamp = resample(txSig, P, Q);

% Channel
delaySig = step(DELAY, txResamp);
txSigPower = 1 / sqrt(L);
rxSeq = awgn(delaySig, EsN0, txSigPower);

% Rx matched filter (MF)
mfOut = step(RXFILT, rxSeq);

%% Symbol Timing Recovery
% Downsampled symbols without symbol timing recovery
rxNoSync = downsample(mfOut, L);

% Downsampled symbols with perfect symbol timing recovery
rxPerfectSync = downsample(mfOut, L, timeOffset);

% Our symbol timing recovery implementation
[ rxSync1 ] = symbolTimingSync(TED, intpl, L, rxSeq, mfOut, K1, K2, ...
    const, Ksym, rollOff, rcDelay, debug_tl_static, debug_tl_runtime);
% MATLAB's implementation
rxSync2 = step(SYMSYNC, mfOut);

%% MER Measurement and Constellation Plots
skip = 0.2 * nSymbols; % skip the initial transitory when plotting

fprintf("\nMeasured MER:\n")
fprintf("No Timing Correction: %.2f dB\n", mer(rxNoSync(skip:end)))
fprintf("Ideal Timing Correction: %.2f dB\n", mer(rxPerfectSync(skip:end)))
fprintf("Our %s Timing Recovery: %.2f dB\n", TED, mer(rxSync1(skip:end)))
fprintf("MATLAB's %s Timing Recovery: %.2f dB\n", ...
    matlabTed, mer(rxSync2(skip:end)))

if (debug_tl_static)
    figure(1)
    scatterplot(rxNoSync(skip:end))
    title('No Timing Correction');
    figure(2)
    scatterplot(rxPerfectSync(skip:end))
    title('Ideal Timing Correction');
    figure(3)
    scatterplot(rxSync1(skip:end))
    title(sprintf('Our %s Timing Recovery', TED));
    figure(4)
    scatterplot(rxSync2(skip:end))
    title(sprintf('MATLAB''s %s', matlabTed));
end
%% ===================== BER / SER COMPUTATION =====================

% 1. Bỏ transient do RRC (rất quan trọng)
% Group delay của cặp RRC Tx + Rx (tính theo SYMBOL)
grpDelaySym = rcDelay;  

rxSync1 = rxSync1(grpDelaySym+1:end);
rxSync2 = rxSync2(grpDelaySym+1:end);
rxPerfectSync = rxPerfectSync(grpDelaySym+1:end);

% 2. Căn chỉnh độ dài
Nmin = min([ ...
    length(data), ...
    length(rxSync1), ...
    length(rxSync2), ...
    length(rxPerfectSync) ...
]);

data_ref = data(1:Nmin);

rxSync1 = rxSync1(1:Nmin);
rxSync2 = rxSync2(1:Nmin);
rxPerfectSync = rxPerfectSync(1:Nmin);

% 3. Demodulation (QUAN TRỌNG: chia lại Ksym)
if (N == 2)
    rxSym1 = qamdemod(rxSync1 / Ksym, M, 'UnitAveragePower', false);
    rxSym2 = qamdemod(rxSync2 / Ksym, M, 'UnitAveragePower', false);
    rxSymP = qamdemod(rxPerfectSync / Ksym, M, 'UnitAveragePower', false);
else
    rxSym1 = pamdemod(real(rxSync1 / Ksym), M);
    rxSym2 = pamdemod(real(rxSync2 / Ksym), M);
    rxSymP = pamdemod(real(rxPerfectSync / Ksym), M);
end

% 4. SER (symbol error rate) – nên xem trước
[~, ser1] = symerr(data_ref, rxSym1);
[~, ser2] = symerr(data_ref, rxSym2);
[~, serP] = symerr(data_ref, rxSymP);

% 5. BER (bit error rate)
k = log2(M);

txBits  = de2bi(data_ref, k, 'left-msb');
rxBits1 = de2bi(rxSym1,   k, 'left-msb');
rxBits2 = de2bi(rxSym2,   k, 'left-msb');
rxBitsP = de2bi(rxSymP,   k, 'left-msb');

[~, ber1] = biterr(txBits(:), rxBits1(:));
[~, ber2] = biterr(txBits(:), rxBits2(:));
[~, berP] = biterr(txBits(:), rxBitsP(:));

% 6. In kết quả
fprintf("\n===== SER =====\n");
fprintf("Custom Timing Sync : %.3e\n", ser1);
fprintf("MATLAB Timing Sync : %.3e\n", ser2);
fprintf("Perfect Timing     : %.3e\n", serP);

fprintf("\n===== BER =====\n");
fprintf("Custom Timing Sync : %.3e\n", ber1);
fprintf("MATLAB Timing Sync : %.3e\n", ber2);
fprintf("Perfect Timing     : %.3e\n", berP);
