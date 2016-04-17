function [medianGauss,x] = plot_Gauss(x,sig,symmetric,show_all_Gauss,show_sig_values)

% Funtion that takes in cortical distances and sigma values, and plot the
% resulting Gaussian
%
%   Usage: plot_Gauss(x,sig,show_all_Gauss,show_sig_values,symmetric)
%
%   Assumes DoG model:
%   x = column vector of distances
%   sig = [sig1 sig2 sig3 sig4]
%
%   Written by Andrew S Bock May 2015

%% set defaults
if ~exist('symmetric','var')
    symmetric = 1;
end
if ~exist('show_all_Gauss','var')
    show_all_Gauss = 0;
end
if ~exist('show_sig_values','var')
    show_sig_values = 0;
end
%% Plot Gauss
tmp1 = -(x.^2);
tmp2 = (2*sig(:,1).^2)';
tmp3 = bsxfun(@rdivide,tmp1,tmp2);
tmpGauss1 = exp(tmp3);
tmpGauss1 = bsxfun(@mtimes,tmpGauss1,sig(:,3)');
tmpGauss1(isnan(tmpGauss1))=0;
% create negative Gaussian (we will subtract, hence 'negative')
tmp1 = -(x.^2);
tmp2 = (2*sig(:,2).^2)';
tmp3 = bsxfun(@rdivide,tmp1,tmp2);
tmpGauss2 = exp(tmp3);
tmpGauss2 = bsxfun(@mtimes,tmpGauss2,sig(:,4)');
tmpGauss2(isnan(tmpGauss2))=0;
% If sigma(2) is not zero, subtract tmpGauss2
tmpGauss = tmpGauss1-tmpGauss2;
% Calc mean Gauss
medianGauss = median(tmpGauss,2);
if show_all_Gauss
    figure;plot(x,tmpGauss1,'.',x,tmpGauss2,'.',x,tmpGauss,'.');
    for s = 1:size(sig,1);
        if s == 1
            legendstr = ['''' num2str(sig(s,[1 3])) '''' ',' ...
                '''' num2str(sig(s,[2 4])) '''' ',' ...
                '''' num2str(sig(s,1:4)) ''''];
        else
            legendstr = [legendstr ',' ...
                '''' num2str(sig(s,[1 3])) '''' ',' ...
                '''' num2str(sig(s,[2 4])) '''' ',' ...
                '''' num2str(sig(s,1:4)) ''''];
        end
    end
end
if symmetric
    x = [-flipud(x);x];
    medianGauss = [flipud(medianGauss);medianGauss];
end
figure;plot(x,medianGauss,'.');
for s = 1:size(sig,1);
    if s == 1
        legendstr = ['''' num2str(sig(s,1:4)) ''''];
    else
        legendstr = [legendstr ',' ...
            '''' num2str(sig(s,1:4)) ''''];
    end
end
if show_sig_values
    eval(['legend(' legendstr ')']);
end
