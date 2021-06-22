
snr = 1;
PER_OFDM = 0;

    
PSDULength = 300; % Number of bytes
NumCarriers     = 64;
NumDataCarriers = 64;
CarrierStep = 4; % Just for SC-FDM, every n-th subcarrier is occupied

N0 = 10.^(-snr/10);
r = 1/2; % Code rate
m = 4; % Modulation order {2, 4, 16, 64}

% Different path gains with their coherence bandwidth
% 20 MHz
% pathdelay = [0 4 7]*1e-7;
% pathgains = [0 -35 -40];
% 10 MHz
% pathdelay = [0 4 7 10]*1e-7;
% pathgains = [0 -30 -36 -40];
% 5 MHz
% pathdelay = [0 4 7 10]*1e-7;
% pathgains = [0 -25 -30 -32];
% 2 MHz
% pathdelay = [0 4 7 10]*1e-7;
% pathgains = [0 -15 -20 -30];
% 1 MHz
pathdelay = [0 4 7 10]*1e-7;
pathgains = [0 -10 -14 -18];

% Maximum Doppler Shift (standard:10)
maxDoppler = 2000;

% Create Rayleigh fading channel
rayChan = comm.RayleighChannel('SampleRate',1e7,'MaximumDopplerShift',maxDoppler,...
    'PathDelays',pathdelay,'AveragePathGains',pathgains,...
    'PathGainsOutputPort',true,'RandomStream','mt19937ar with seed','Seed',5);


MaxNumError = 1e2; % 1e3
MaxPacketTx = 1e4; % 1e5

DataInd = 1:NumCarriers;

% Simulation parameters
Sim_Para.PSDULength         = PSDULength;
Sim_Para.NumCarriers        = NumCarriers;
Sim_Para.NumDataCarriers    = NumDataCarriers;
Sim_Para.N0                 = N0;
Sim_Para.SNR                = snr;
Sim_Para.MaxNumError        = MaxNumError;
Sim_Para.MaxPacketTx        = MaxPacketTx;
Sim_Para.DataInd            = DataInd;
Sim_Para.AWGN               = false; % Defines if just AWGN channel is used or Rayleigh
Sim_Para.Rate               = r;
Sim_Para.M                  = m;    % Modulation order
Sim_Para.NoPilots           = true;
Sim_Para.RayleighChan       = rayChan;
Sim_Para.ImperfectEq        = false;
Sim_Para.diversity          = 'MRC'; % SC, EGC, MRC
Sim_Para.n_link             =3;

  

if Sim_Para.AWGN 
[PER_OFDM]        = LDPC_awgn_simulator(Sim_Para) 
else
[PER_OFDM]        = LDPC_multiconnectivity_simulator(Sim_Para) 
end


