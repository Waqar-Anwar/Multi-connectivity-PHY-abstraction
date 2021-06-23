b = load('Results_OFDM_awgn_mcs1.mat');
chan = load('channel_k0_2k_fr64_L4.mat');

beta_eEESM = 2; % 
snr = -10:1:15; % Number of SNR points
num_ch = 1e4;  % Number of random channel realization

snr_awgn = b.snr_awgn_mcs1;
per_awgn = b.per_awgn_mcs1;

ICI = 3.185e-5; % ICI at 2 KHz Doppler with 156.25 KHz carrier spacing
sigma_sym = abs(sinc(500*(1/156.25e3)))^2;

L = 2; % multi-connectivity order
Com = 1; % 1 for slection combining, 2 EGC and 3 for MRC

alpha_ofdm = 1;
channel= chan.H_est(:,:,1:L,1:num_ch);


for i=1:length(snr)
    
    pn = 10.^(-snr(i)/10);
    %         noise = (ICI + pn)/sigma_sym;
    
    if (L==1)
        snr_real = sum((1/pn).*((abs(channel)).^2),3);
        var_mcs1_l1 = var(10*log10(snr_real),1,[1 2 3]);
        alpha_ofdm = var_mcs1_l1.*-0.00 + 1.23;
        esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eEESM/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eEESM)./sqrt(((2.*snr_real)./beta_eEESM) + 1),[1 2]).^(-2))) - 1));
        loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
        per_pla_mcs1_l1_dsf(i)=mean(per_awgn(loc_d));
    elseif (L==2)
        if(Com ==1)
            channel_sc = max(((abs(channel)).^2),[],3);
            snr_real = (1/pn).*channel_sc;
            var_mcs1_l2 = var(10*log10(snr_real),1,[1 2 3]);
            alpha_ofdm = var_mcs1_l2*-0.000 + 1.21;
            esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eEESM/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eEESM)./sqrt(((2.*snr_real)./beta_eEESM) + 1),[1 2]).^(-2))) - 1));
            loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
            per_pla_mcs1_sc_l2_dsf(i)=mean(per_awgn(loc_d));
        elseif(Com ==2)
            channel_egc = ((sum(abs(channel),3)).^2)/L;
            snr_real = (1/pn).*channel_egc;
            var_mcs1_l2 = var(10*log10(snr_real),1,[1 2 3]);
            alpha_ofdm = var_mcs1_l2*-0.000 + 1.21;
            esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eEESM/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eEESM)./sqrt(((2.*snr_real)./beta_eEESM) + 1),[1 2]).^(-2))) - 1));
            loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
            per_pla_mcs1_egc_l2_dsf(i)=mean(per_awgn(loc_d));
        else
            snr_real = sum((1/pn).*((abs(channel)).^2),3);
            var_mcs1_l2 = var(10*log10(snr_real),1,[1 2 3]);
            alpha_ofdm = var_mcs1_l2*-0.000 + 1.21;
            esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eEESM/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eEESM)./sqrt(((2.*snr_real)./beta_eEESM) + 1),[1 2]).^(-2))) - 1));
            loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
            per_pla_mcs1_mrc_l2_dsf(i)=mean(per_awgn(loc_d));
        end
    elseif (L==4)
        if(Com == 1)
            channel_sc = max(((abs(channel)).^2),[],3);
            snr_real = (1/pn).*channel_sc;
            var_mcs1_l4 = var(10*log10(snr_real),1,[1 2 3]);
            alpha_ofdm = var_mcs1_l4*-0.0007 + 1.195;
            esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eEESM/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eEESM)./sqrt(((2.*snr_real)./beta_eEESM) + 1),[1 2]).^(-2))) - 1));
            loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
            per_pla_mcs1_sc_l4_dsf(i)=mean(per_awgn(loc_d));
        elseif (Com == 2)
            channel_egc = ((sum(abs(channel),3)).^2)/L;
            snr_real = (1/pn).*channel_egc;
            var_mcs1_l4 = var(10*log10(snr_real),1,[1 2 3]);
            alpha_ofdm = var_mcs1_l4*-0.0007 + 1.195;
            esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eEESM/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eEESM)./sqrt(((2.*snr_real)./beta_eEESM) + 1),[1 2]).^(-2))) - 1));
            loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
            per_pla_mcs1_egc_l4_dsf(i)=mean(per_awgn(loc_d));
        else
            snr_real = sum((1/pn).*((abs(channel)).^2),3);
            var_mcs1_l4 = var(10*log10(snr_real),1,[1 2 3]);
            alpha_ofdm = var_mcs1_l4*-0.0007 + 1.195;
            esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eEESM/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eEESM)./sqrt(((2.*snr_real)./beta_eEESM) + 1),[1 2]).^(-2))) - 1));
            loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
            per_pla_mcs1_mrc_l4_dsf(i)=mean(per_awgn(loc_d));
        end
        
    end
    
end