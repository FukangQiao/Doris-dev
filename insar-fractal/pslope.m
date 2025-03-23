function [C0,beta] = pslope(freq,power);
% [C0,beta] = pslope(freq,power);
% 
% function to derive the slope beta and C0 of an exponential function
% in loglog  scale
%
% S(k)  = C0^2 .* k^(-beta)
% power = C0^2 .* freq^(-beta)
%
% Ramon Hanssen, April 2000

if nargin==0, help pslope;break;end

% Don't worry about orientation
  freq   = freq(:);
  power  = power(:);

% Check if we've got a zero frequency. If yes, remove it
  if find(freq==0), 
    notzero = find(freq~=0);
    freq    = freq (notzero);
    power   = power(notzero);
  end

  logf       = log10(freq);
  logp       = log10(power);

  [p, dummy] = polyfit(logf,logp,1);
  beta       = -p(1); 

% find the C0 value (for freq = 1, or logf = 0 )
  position   = interp1(logf,[1:length(logf)],0) ;
  logC0      = interp1([1:length(logp)],logp,position) ;
  C0squared  = 10^logC0;
  C0         = sqrt(C0squared);
