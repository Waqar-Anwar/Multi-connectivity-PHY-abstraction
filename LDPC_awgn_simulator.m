function [PER]= LDPC_awgn_simulator(Sim_Para)
%% This function simulates a transmission of LDPC encoded data.
% BER and PER are calculated.

r = Sim_Para.Rate;
m = Sim_Para.M;
numDC = Sim_Para.NumDataCarriers;
numData = Sim_Para.PSDULength * 8;

% Rate adaption to preserve integer as output length
if round(numData/r) ~= numData/r
    r = numData/round(numData/r);
end

% Create AWGN channel
awgn = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SignalPower', 1);
awgn.SNR = Sim_Para.SNR;

numPadbit = 16; 
hError = comm.ErrorRate;


if m == 2
    modulation = 'pi/2-BPSK';
    gfdm.mu = 2;
elseif m == 4
    modulation = 'QPSK';
else
    modulation = strcat(num2str(m),'QAM');
end

rv = 0; % Redundancy version for rate match, {0, 1, 2, 3}
nlayers = 1; % Number of layers

% Initialization of LDPC
cbsInfo = nrDLSCHInfo(Sim_Para.PSDULength*8,r);

Packet_error = 0;

%% Main processing loop
for jj=1:Sim_Para.MaxPacketTx
    
    dataIn = randi([0 1], Sim_Para.PSDULength*8, 1);
    
    % CRC attachement and Code Block Segmentation
    tbIn = nrCRCEncode(dataIn,cbsInfo.CRC);
    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
    
    % LDPC encoding
    dataEncode = nrLDPCEncode(cbsIn,cbsInfo.BGN);
    
    % Perform rate match
    dataEncode = nrRateMatchLDPC(dataEncode,Sim_Para.PSDULength*8/r,rv,modulation,nlayers);
    dataEncode = dataEncode(:); % Serial data stream for modulation
    
    % Modulation
    modData = nrSymbolModulate(dataEncode,modulation);
    
    modData = complex(modData); % for correct calculation with BPSK
    
    % Channel passing
    
    % Just AWGN channel
    EqRxData = awgn(modData);
    RxData = EqRxData(:);
    
    
    % Demodulation
    demodData = nrSymbolDemodulate(RxData,modulation,Sim_Para.N0);
    
    ulschLLRs=demodData;
    
    % Rate recovery
    rxcodedcbs = nrRateRecoverLDPC(ulschLLRs, ...
        Sim_Para.PSDULength*8,r,rv,modulation,nlayers);
    
    
    % LDPC decoding and code block desegmentation
    rxcbs = nrLDPCDecode(rxcodedcbs,cbsInfo.BGN,25);
    dataOut = nrCodeBlockDesegmentLDPC(rxcbs,cbsInfo.BGN,Sim_Para.PSDULength*8+cbsInfo.L);
    
    Error = isequal(dataIn,dataOut(1:Sim_Para.PSDULength*8));
    
    bit_Error = hError(int8(dataIn),dataOut(1:Sim_Para.PSDULength*8));
    if(~Error)
        Packet_error = Packet_error+1;
    end
    
    if(Packet_error >= Sim_Para.MaxNumError)
        break;
    end
end
PER = Packet_error/jj;

end



