function [fline] = fracline(N,beta,doplot,P0);
% [fline] = fracline(N,beta,[doplot,P0]);
%
%  Script to create a 1 dimensional fractal line 
%  with a power law behavior
%
%  beta : power law exponent for a 1D  data
%         1D:  S(k) = C0 k^(-beta)
%     [   2D:  S(k) = C1 k^(-(beta+1))  ]
%
%  N    : length of the line
%  doplot : 'y'/'n' show figure (default 'n')
%  [P0]     : multiply AMPLITUDE (!) spectrum by sqrt(P0) (default P0 = 1)
%
%  fline: real-valued output line
% 
% e.q., [fline] = fracline(128, 1.3,'y');
%
% Ramon Hanssen, May 2000
% RH 25-May-2000 13:11 : Added scaling with P0

if nargin==0,help fracline;break;end
if nargin<4, P0 = sqrt(1) ;end % NOTE: for power we use P(k) = P0 x k^(-beta)
                               % Here we use amplitude only, so sqrt(P0) in stead of P0
if nargin<3, doplot = 'n' ;end

% Check if N is even
if rem(N,2), fprintf(1,'Error: use an even value for N\n');break;end

% Simulate a uniform random signal

  h      = rand(N,1) ;
  % demean
  h = h - mean(h);
  %H      = fftshift(fft(h));
  H      = fft(h);

% scale the spectrum with the power law
  X      = [-N/2: (N/2 -1)]';
%  [X,Y]  = meshgrid(x,x);
  k      = sqrt(X.^2);
% beta   = beta+1;             % beta+1 is used as beta, since, the power exponent
                               % is defined for a 1D slice of the 2D spectrum:
                               % austin94: "Adler, 1981, shows that the surface profile 
                               %   created by the intersection of a plane and a
                               %   2-D fractal surface is itself fractal with 
                               %   a fractal dimension  equal to that of the 2D 
                               %   surface decreased by one."
  noemer = k .^(beta/2);       % The power beta/2 is used because the power spectral
                               % density is proportional to the amplitude squared 
                               % Here we work with the amplitude, instead of the power
                               % so we should take sqrt( k.^beta) = k.^(beta/2)  RH

% prevent dividing by zero;
  noemer(find(noemer==0))=1;   % This means that the zero frequency is NOT weighted

  Hnew   = sqrt(P0) .* H ./ ifftshift(noemer);
% TEST OM Hnew(2) [dat is k=1] gelijk te krijgen aan P0
  Hnew = (Hnew ./ Hnew(2)) .* sqrt(P0);

% Create spectral line by ifft
  fline = abs(ifft(Hnew));
% Remove mean to get zero-mean data
  fline = fline - mean(fline(:));

fprintf(1,'TEST\n');

      %[abs(H), abs(Hnew), ifftshift( noemer)]
      P0
      % Variance of the line
      var(fline)
      % This value has to be equal to P0
      Hnew(2).*conj(Hnew(2))
      
      (1/(2*N*(N-1)))*sum(Hnew.*conj(Hnew))
fprintf(1,'EINDE TEST\n');

  if strcmp(doplot,'y'), 
    figure(1);plot(fline);
  end

