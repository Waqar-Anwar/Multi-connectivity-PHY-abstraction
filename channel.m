function [txdata, H_out] = channel(Channel_in, data, Sim_Para)

N = length(data(:,1)); % Number of fft points

% if nargin > 2
%     if strcmp(Sim_Para.Waveform,'GFDM')
%         N = length(data(:,1)); % GFDM needs higher frequency solution
%     else
%         s = size(data);
%         N = s(1);
%     end
% end


Sig = ones(length(data(1,:))*length(data(:,1)),1);

[~,pathgains] = Channel_in(Sig);
dim = size(pathgains);
numTaps = dim(2);

Path_re = reshape(pathgains,[length(data(:,1)),length(data(1,:)),numTaps]);
b(:,1)= Path_re(floor(N/2),:,1);
b(:,5)=Path_re(floor(N/2),:,2);
b(:,8)=Path_re(floor(N/2),:,3);
if numTaps == 4
    b(:,11)=Path_re(floor(N/2),:,4);
end

for i = 1:length(b(:,1))
    H_est(:,i) = fft(b(i,:),N);
end

H_out = H_est(1:length(data(:,1)),:);
txdata = data.*H_out;
end