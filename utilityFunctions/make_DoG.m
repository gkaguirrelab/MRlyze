function [DoG_Gauss] = make_DoG(x,sig)
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
% create DoG
DoG_Gauss = tmpGauss1-tmpGauss2;