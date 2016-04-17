function [out_tc] = shift_tc(in_tc,shift_time,sample_rate)

%% Set to double (if not already)
in_tc = double(in_tc);
mean_tc = mean(in_tc);
for i = 1:size(in_tc,2);
    in_tc(:,i) = in_tc(:,i) - mean_tc(i);
end
%% set the sample rate
if ~exist('sample_rate','var')
    sample_rate=0.05;
end
%% Upsample the timecourse
% Linear
X = 1:1:size(in_tc,1);
Xq = 1:sample_rate:size(in_tc,1);
% Note the transpose for sinc_interp
tmp_tc = sinc_interp(in_tc',X,Xq);
tmp_tc = tmp_tc';
%% Shift the timecourse
shift_block = round(abs(shift_time)/sample_rate);
%shift_section = repmat(mean(tmp_tc,1),shift_block,1); % fill with the mean
if shift_time > 0; % shift earlier in time
    shift_section = tmp_tc(1:shift_block,:);
    new_tc = [tmp_tc(shift_block+1:end,:);shift_section];
elseif shift_time < 0; % shift later in time
    shift_section = tmp_tc(end-(shift_block-1):end,:);
    new_tc = [shift_section;tmp_tc(1:end-shift_block,:)];
end
%% Downsample the upsampled timecourse
if shift_time ~=0
    out_tc = downsample(new_tc,1/sample_rate);
    for i = 1:size(new_tc,2);
        out_tc(:,i) = out_tc(:,i) + mean_tc(i);
    end
else
    out_tc = in_tc;
end