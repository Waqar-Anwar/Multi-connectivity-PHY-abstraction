function [PER]= LDPC_multiconnectivity_simulator(Sim_Para)
%% This function simulates a transmission of LDPC encoded data.
% BER and PER are calculated.


% Create Rayleigh fading channel
rayChan = Sim_Para.RayleighChan;


r = Sim_Para.Rate;
m = Sim_Para.M;
numDC = Sim_Para.NumDataCarriers;
numData = Sim_Para.PSDULength * 8;

% Rate adaption to preserve integer as output length
if round(numData/r/log2(m)) ~= numData/r/log2(m)
    r = numData/log2(m)/round(numData/r/log2(m));
end

% Create AWGN channel
awgn = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SignalPower', 1);
awgn.SNR = Sim_Para.SNR;

% Calculation of the number of padding bits

numPadbit = ceil(numData/numDC/log2(m)) * numDC*log2(m) - numData;

lengthModData = (numData+numPadbit)/log2(m)/r;

hError = comm.ErrorRate;


if m == 2
    modulation = 'pi/2-BPSK';
elseif m == 4
    modulation = 'QPSK';
else
    modulation = strcat(num2str(m),'QAM');
end

rv = 0; % Redundancy version for rate match, {0, 1, 2, 3}
nlayers = 1; % Number of layers

% Initialization of LDPC
cbsInfo = nrDLSCHInfo(numData,r);


Packet_error = 0;

%% Main processing loop
for jj=1:Sim_Para.MaxPacketTx
    
    
    dataIn = randi([0 1], numData, 1);
    
    
    % CRC attachement and Code Block Segmentation
    tbIn = nrCRCEncode(dataIn,cbsInfo.CRC);
    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
    
    % LDPC encoding
    dataEncode = nrLDPCEncode(cbsIn,cbsInfo.BGN);
    
    % Perform rate match
    if floor(cbsInfo.K/r) ~= ceil(cbsInfo.K/r)
        r = cbsInfo.K / round(cbsInfo.K/r);
    end
    dataEncode = nrRateMatchLDPC(dataEncode,cbsInfo.K/r,rv,modulation,nlayers);
    
    dataEncode = dataEncode(:); % Serial data stream for modulation
    
    
    % Modulation
    modData = nrSymbolModulate(dataEncode,modulation);
    
    
    padBitDC = ceil(length(modData)/numDC)*numDC - length(modData);
    modData(end+1:end+padBitDC) = 0;
    modData=reshape(modData,numDC,[]);
    
    % Rayleigh fading channel + AWGN
    for l=1:Sim_Para.n_link
        
        release(rayChan);
        rayChan.Seed=l*jj;
        
        [~,H_est(:,:,l)] = channel(rayChan,modData,Sim_Para);
    end
    
    if(strcmp(Sim_Para.diversity,'SC'))
        H_sc = max(H_est,[],3);
        rxNoiselessWaveform = modData.*H_sc;
        rxWaveform = awgn(rxNoiselessWaveform);
        
        T_MMSE = conj(H_sc)./(conj(H_sc).*H_sc + Sim_Para.N0(it));
        EqRxData = T_MMSE .* rxWaveform;
        csi = conj(H_sc).*H_sc + Sim_Para.N0(it);
        csi = csi(1:end-padBitDC).';
        
    elseif(strcmp(Sim_Para.diversity,'MRC'))
        rxNoiselessWaveform = repmat(modData,1,1,Sim_Para.n_link).*H_est;
        rxNoiselessWaveform = reshape(rxNoiselessWaveform,64,[]);
        rxWaveform = awgn(rxNoiselessWaveform);
        
        rxWaveform = reshape(rxWaveform,64,[],Sim_Para.n_link);
        norm_f = sum(conj(H_est).*H_est,3);
        EqRxData = sum(conj(H_est).*rxWaveform,3)./norm_f;
        csi = norm_f;
        csi = csi(1:end-padBitDC).';
    else
        rxNoiselessWaveform = repmat(modData,1,1,Sim_Para.n_link).*abs(H_est);
        rxNoiselessWaveform = reshape(rxNoiselessWaveform,64,[]);
        rxWaveform = awgn(rxNoiselessWaveform);
        rxWaveform = reshape(rxWaveform,64,[],Sim_Para.n_link);
        EqRxData = sum(rxWaveform,3)./Sim_Para.n_link;
        csi = mean(abs(H_est),3).*Sim_Para.n_link;
        csi = csi(1:end-padBitDC).';
        
    end
    
    
    % Demodulation
    RxData = EqRxData(1:end-padBitDC).';
    
    demodData = nrSymbolDemodulate(RxData,modulation,Sim_Para.N0);
    
    csi = nrLayerDemap(csi);
    Qm = length(demodData) / length(RxData);
    csi = reshape(repmat(csi{1}.',Qm,1),[],1);
    ulschLLRs = demodData .* csi;
    
    
    % Rate recovery
    rxcodedcbs = nrRateRecoverLDPC(ulschLLRs,numData,r,rv,...
        modulation,nlayers);
    
    % LDPC decoding and code block desegmentation
    rxcbs = nrLDPCDecode(rxcodedcbs,cbsInfo.BGN,100);
    dataOut = nrCodeBlockDesegmentLDPC(rxcbs,cbsInfo.BGN,numData+cbsInfo.L);
    
    Error = isequal(dataIn(1:numData),dataOut(1:numData));
    
    %     bit_Error = hError(int8(dataIn(1:numData)),dataOut(1:numData));
    if(~Error)
        Packet_error = Packet_error+1;
    end
    
    
    if(Packet_error >= Sim_Para.MaxNumError)
        break;
    end
    
end
PER = Packet_error/jj;

% BER = bit_Error(1);

end

