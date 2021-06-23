function per = run_fading_simulator(sim)

awgn_result = load('AWGN_look_up_300_bytes.mat');

beta = [2 10 10 42]; % default value of beta: 2 for QPSK, 10 for 16-QAM, 42 for 16-QAM and 170 for 64-QAM

beta_eesm = beta(sim.mcs);

snr_awgn = awgn_result.snr_awgn(sim.mcs,:);
per_awgn = awgn_result.per_awgn(sim.mcs,:);
 
    
    
for i=1:length(sim.snr)
        
        noise = 10.^(-sim.snr(i)/10);
        pn = (sim.ICI + noise)/sim.sigma_sym;
        if (sim.L==1)
            snr_real = sum((1/pn).*((abs(sim.channel)).^2),3);
            var_mcs1_l1 = var(10*log10(snr_real),1,[1 2 3]);
            alpha_ofdm = var_mcs1_l1.*-0.00 + 1.23;
            esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eesm/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eesm)./sqrt(((2.*snr_real)./beta_eesm) + 1),[1 2]).^(-2))) - 1));
            loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
            per(i)=mean(per_awgn(loc_d));
        elseif (sim.L==2)
            if(sim.combining ==1)
                channel_sc = max(((abs(sim.channel)).^2),[],3);
                snr_real = (1/pn).*channel_sc;
                var_mcs1_l2 = var(10*log10(snr_real),1,[1 2 3]);
                alpha_ofdm = var_mcs1_l2*-0.000 + 1.21;
                esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eesm/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eesm)./sqrt(((2.*snr_real)./beta_eesm) + 1),[1 2]).^(-2))) - 1));
                loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
                per(i)=mean(per_awgn(loc_d));
            elseif(sim.combining ==2)
                channel_egc = ((sum(abs(sim.channel),3)).^2)/sim.L;
                snr_real = (1/pn).*channel_egc;
                var_mcs1_l2 = var(10*log10(snr_real),1,[1 2 3]);
                alpha_ofdm = var_mcs1_l2*-0.000 + 1.21;
                esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eesm/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eesm)./sqrt(((2.*snr_real)./beta_eesm) + 1),[1 2]).^(-2))) - 1));
                loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
                per(i)=mean(per_awgn(loc_d));
            elseif (sim.combining == 3)
                snr_real = sum((1/pn).*((abs(sim.channel)).^2),3);
                var_mcs1_l2 = var(10*log10(snr_real),1,[1 2 3]);
                alpha_ofdm = var_mcs1_l2*-0.000 + 1.21;
                esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eesm/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eesm)./sqrt(((2.*snr_real)./beta_eesm) + 1),[1 2]).^(-2))) - 1));
                loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
                per(i)=mean(per_awgn(loc_d));
            else
                error('Combining scheme is not a valid scheme, please choose 1 for SC, 2 for EGC and 3 for MRC')
            end
        elseif (sim.L==4)
            if(sim.combining == 1)
                channel_sc = max(((abs(sim.channel)).^2),[],3);
                snr_real = (1/pn).*channel_sc;
                var_mcs1_l4 = var(10*log10(snr_real),1,[1 2 3]);
                alpha_ofdm = var_mcs1_l4*-0.0007 + 1.195;
                esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eesm/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eesm)./sqrt(((2.*snr_real)./beta_eesm) + 1),[1 2]).^(-2))) - 1));
                loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
                per(i)=mean(per_awgn(loc_d));
            elseif (sim.combining == 2)
                channel_egc = ((sum(abs(sim.channel),3)).^2)/sim.L;
                snr_real = (1/pn).*channel_egc;
                var_mcs1_l4 = var(10*log10(snr_real),1,[1 2 3]);
                alpha_ofdm = var_mcs1_l4*-0.0007 + 1.195;
                esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eesm/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eesm)./sqrt(((2.*snr_real)./beta_eesm) + 1),[1 2]).^(-2))) - 1));
                loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
                per(i)=mean(per_awgn(loc_d));
            elseif(sim.combining)
                snr_real = sum((1/pn).*((abs(sim.channel)).^2),3);
                var_mcs1_l4 = var(10*log10(snr_real),1,[1 2 3]);
                alpha_ofdm = var_mcs1_l4*-0.0007 + 1.195;
                esnr_mrc_mcs1 = 10.*log10((alpha_ofdm.*beta_eesm/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta_eesm)./sqrt(((2.*snr_real)./beta_eesm) + 1),[1 2]).^(-2))) - 1));
                loc_d = knnsearch(snr_awgn',esnr_mrc_mcs1(:));
                per(i)=mean(per_awgn(loc_d));   
            else
                error('Combining scheme is not a valid scheme, please choose 1 for SC, 2 for EGC and 3 for MRC')
            end
        else
            error('Please choose multi-connectivity order 1, 2 or 4')    
        end
     end
 