function [out] = gammapdf(in,g,beta)
%
% gamma pdf forumla
%        ((x-u)/beta)^(g-1) * exp (-(x-u)/beta)
% g(x)=  -------------------------------------
%                   beta*gammafunction(g)
% or
% log(g(x)) = (g-1)*log((x-u)) + (-((x-u)/beta)) - loggamma(g) -
% g*log(beta)

if length(g) == 1
  g = ones(size(in))*g;
end
if length(beta) == 1
  beta = ones(size(in))*beta;
end

out = zeros(size(in));
ii = find(in > 0);

out(ii) = (g(ii)-1) .* ...
          log(in(ii))-(in(ii)./beta(ii))-gammaln(g(ii))-g(ii).*log(beta(ii));
out = exp(out(ii));

ii = find(in == 0 & g < 1);
out(ii) = Inf;
ii = find(in == 0 & g == 1);
out(ii) = 1./beta(ii);