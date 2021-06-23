chan = load('channel_k0_2k_fr64_L4.mat');


sim.snr = -10:1:15; % number of SNR points
sim.num_ch = 1e2; % % Number of random channel realization
sim.L = 4; % multi-connectivity order, choose 1, 2 or 4
sim.combining = 3; % 1 for slection combining, 2 EGC and 3 for MRC

sim.mcs = 1;   % 1 for 1/2 QPSK, 2 for 1/2 16-QAM, 3 for 3/4 16-QAM, 4 for 3/4 64-QAM

fading = 4 ; % 1 for DFF (doubly flat fading), 2 for FSF (frequency selective fading), 3 for TSF (time selective fadding), 4 for DSF (doubly flat fadding)

if fading == 1;
    sim.channel= chan.H_est(1,1,1:L,1:num_ch);
    sim.ICI = 1.2882e-06; % ICI at 100 Hz Doppler with 156.25 KHz carrier spacing
    sim.sigma_sym = abs(sinc(100*(1/156.25e3)))^2;
elseif fading == 2;
    sim.channel= chan.H_est(:,1,1:L,1:num_ch);
    sim.ICI = 1.2882e-06 ; % ICI at 100 Hz Doppler with 156.25 KHz carrier spacing
    sim.sigma_sym = abs(sinc(100*(1/156.25e3)))^2;
elseif fading == 3;
    sim.channel= chan.H_est(1,:,1:L,1:num_ch);
    sim.ICI = 5.1286e-04; % ICI at 2 KHz Doppler with 156.25 KHz carrier spacing
    sim.sigma_sym = abs(sinc(2000*(1/156.25e3)))^2;
elseif fading == 4;
    sim.channel= chan.H_est(:,:,1:L,1:num_ch);
    sim.ICI = 5.1286e-04; % ICI at 2 KHz Doppler with 156.25 KHz carrier spacing
    sim.sigma_sym = abs(sinc(2000*(1/156.25e3)))^2;
else
    error('Fading type is not a valid number please choose from 1-4')
end

per = run_fading_simulator (sim)