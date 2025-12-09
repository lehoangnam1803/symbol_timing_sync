clear; clc; close all

P = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
     25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
     25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
      9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
     24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
      2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];
blockSize = 27;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);
cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);
cfgLDPCDec = ldpcDecoderConfig(pcmatrix);

M = 64;
maxnumiter = 10;
snr = [0.5 1 1.5 2];
numframes = 20000;

for ii = 1:length(snr)
    data = randi([0 1],cfgLDPCEnc.NumInformationBits,numframes,'int8');

    % Transmit and receive with LDPC coding
    encodedData = ldpcEncode(data,cfgLDPCEnc);
    modSignal = qammod(encodedData,M,InputType='bit');
    [rxsig] = awgn(modSignal,snr(ii));
    llrOut = qamdemod(rxsig,M, ...
        OutputType='approxllr');%, ... % Or use 'llr' for exact but slower LLR calculation
        %NoiseVariance=noisevar);
    rxbits = ldpcDecode(llrOut,cfgLDPCDec,maxnumiter);
    fprintf(['SNR = %2.1f\n   Coded: Error rate = %1.6f, ' ...
        'Number of errors = %d\n'], ...
        snr(ii),nnz(data~=rxbits)/numel(data),nnz(data~=rxbits));

    % Transmit and receive with no LDPC coding
    noCoding = qammod(data,M,InputType='bit');
    rxNoCoding = awgn(noCoding,snr(ii));
    rxBitsNoCoding = qamdemod(rxNoCoding,M,OutputType='bit');

    fprintf(['Noncoded: Error rate = %1.6f, ' ...
        'Number of errors = %d\n\n'], ...
        nnz(data~=rxBitsNoCoding)/numel(data),nnz(data~=rxBitsNoCoding))
end