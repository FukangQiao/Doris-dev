function [freq,powerv,powerh] = anafrac2D(fsurf,dx,figstart,rotational)
%  [freq,powerv,powerh] = anafrac2D(fsurf,dx,[figstart],[rotational])
%
%  Function to analyze the fractal characteristics of a 2D surface
%
% Draws the 1D/2D rotationally averaged power spectra of 
%    fsurf (2Dsurface) with
%    dx    sampling interval
%    [figstart] Starting number for figures  fig = 0 means no figures!!
%    [rotational] ('y','n') plot rotationally averaged spectrum on top
%                 (default = 'y')
%
% For example:  Create fractal surface with 
%     N=128;beta=8/3;dx=0.02;% 20 m sampling
%     [fsurf] = fracsurf(N,beta);
%   or 
%     N=128;beta=[5/3,8/3,2/3];regime=[60,90,100];dx=0.02;% 20 m sampling
%        [fsurf] = fracsurfatmo(N,beta,regime);
%   Then,  Analyze it using
%     anafrac2D(fsurf,dx)
%
% RH 21-May-2000 17:02 : Created

if nargin==0;help anafrac2D;break;end
if nargin==2;
  if  isempty(get(0,'Children')),
    figstart = 1;
  else
    figstart = max(get(0,'Children'))+1;
  end
end
if nargin<4, 
  rotational = 'y';
end

  NR = size(fsurf,1);
  NC = size(fsurf,2);
  if NR == NC,
     N = NR;
     L = [0 : dx : (N-1)*dx];
  else
     fprintf(1,'At the moment, only use square input size\n');
     break
  end
 
if figstart ~=0,
  figure(figstart);imagesc(L,L,fsurf);colorbar
  xlabel(['x [km]']);ylabel(['y [km]']);
end

%%%% 1D spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1D spectra ( horizontal )
  fsh = fft(fsurf')';
  % Periodogram spectral estimator
  powerh = (1/length(fsh)) * mean(( fsh .* conj(fsh) ));

% 1D spectra (vertical)
  fsv = (fft(fsurf))';
  % Periodogram spectral estimator
  powerv = (1/length(fsv)) * mean( fsv .* conj(fsv) );

% The frequency coordinate
  Nyquist = 1/(2*dx);
  df      = 1/(N*dx);
  freq  = [-Nyquist : df : ( Nyquist - df )];
  freq  = fftshift(freq);

% Use only positive part of spectrum
  powerh = powerh(1: N/2);
  powerv = powerv(1: N/2);
  freq   = freq  (1: N/2);

if figstart ~=0,
  figure(figstart+1); loglog(freq,powerh);
  hold on;
    loglog(freq,powerv,'r');
  hold off
  xlabel(['Cycles per km']);ylabel(['Power']);
end

% Calculate slopes from spectra
  [C0h,betah] = pslope(freq,powerh);
  D2h         = (7-betah)/2;
    fprintf(1,'Horizontal: C0 = %6.6f, beta = %6.3f, Fractal dim = %3.1f\n',C0h,betah,D2h);
  [C0v,betav] = pslope(freq,powerv);
  D2v         = (7-betav)/2;
    fprintf(1,'Vertical  : C0 = %6.6f, beta = %6.3f, Fractal dim = %3.1f\n',C0v,betav,D2v);


%%%% 2D spectra (rotationally averaged) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
if strcmp(rotational,'y'),
  fprintf(1,'Calculating rotionally averaged spectrum\n');

  fs2D    = fftshift(fft2(fsurf));
  freq    = [-Nyquist : df : ( Nyquist - df )];
  [FX,FY] = meshgrid(freq,freq);
  k       = sqrt( FX.^2 + FY.^2 );
  pow2D   = ( 1/N^2 ) .* (fs2D .* conj(fs2D) );
  pow2D1D = zeros(N/2-1,1);
  for i = 1 : N/2-1
     sel        = find( k >= i * df & k < (i+1) * df ); 
     pow2D1D(i) = mean(pow2D(sel));
  end
  freqk = [1:N/2-1] *df;

  [C0r,betar] = pslope(freqk,pow2D1D);
  D2r         = (7-betar+1)/2;
  fprintf(1,'Rotational: C0 = %6.6f, beta = %6.3f, Fractal dim = %3.1f\n',C0r,betar,D2r);

  if strcmp(rotational,'y'),
    if figstart ~=0,
    figure(figstart+1);hold on;loglog(freqk, pow2D1D,'m');hold off
    legend(['Hor, (',num2str(betah),')'],...
           ['Ver, (',num2str(betav),')'],...
           ['Rot, (',num2str(betar),')']);
    end
  else
    legend(['Hor, (',num2str(betah),')'],...
           ['Ver, (',num2str(betav),')']);
  end
end

% DRAW THE DIAGONAL POWER LAW LINES
if figstart ~=0,
  figure(figstart+1); 
   % fix the axis limits
     XLIM = get(gca,'xlim');
     YLIM = get(gca,'ylim');
     set(gca,'xlim',XLIM)
     set(gca,'ylim',YLIM)
   % plot power lines
     hold on
     for i = -20:25
       % turb = 10^i*XLIM.^(-5/3); % voor regime III
         turb = 10^i*XLIM.^(-8/3); % voor regime II
       h=plot(XLIM,turb);
       set(h,'linestyle',':');
     end
     hold off
end



% Calculate variace fractal surface
  fprintf(1,'Variance fractal surface               : %f\n',cov(fsurf(:)) );
  fprintf(1,'Variance cross-section fractal surface : %f\n',cov(fsurf(34,:)) );
%figure(3);hist(s(:),20)

