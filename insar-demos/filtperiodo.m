% FILTPERIODO 
%    compare periodogram, fft over 1 line, mean fft for some master/slave
%    conclude: periodogram is no improvement.
%    (mean periodogram not tested.)
%    use periodogram to estimate peak
%    Please set variables in this script
%
%    See also INSARDEMOS, RANGEDEMO, FILTERBLOCK, MYHAMMING,
%    MYRECT
%

% $Revision: 1.1 $  $Date: 2000/03/27 14:05:40 $


%%% Initialize
if (~exist('m','var') & ~exist('s','var')) clear; end; close all force; more off;
showperiodo=1;
filtadaptive=0;


L = 128;%	number of lines
P = 512;
master = 'master2_128_512';
slave  = 'slave2_128_512';

%%% Read images from disk if not present
if ( ~ exist('m','var') )
  fid   = fopen(master);
  Mtmp  = fread(fid,[2*P L],'short');	% orig. SLC
  fclose(fid);
  Mtmp  = Mtmp.';
  Mreal = Mtmp(:,1:2:2*P);
  Mimag = Mtmp(:,2:2:2*P);
  m     = Mreal + i*Mimag;
  clear Mimag Mreal Mtmp;
end

if ( ~ exist('s','var') )
  fid   = fopen(slave);
  Stmp  = fread(fid,[2*P L],'float');	% resampled
  fclose(fid);
  Stmp  = Stmp.';
  Sreal = Stmp(:,1:2:2*P);
  Simag = Stmp(:,2:2:2*P);
  s     = Sreal + i*Simag;
  clear Simag Sreal Stmp;
end


% plot variables
sub = 10;%			do a plot per sub
subplots = floor(P/sub);%	number of subplots
xaxis = 0:P-1;
% ffthifted xaxis = xaxis - P/2;


%%% Test interferogram gen. by pointwize mult.
disp ('press key to continue with complex interferogram');
pause;
close all force;
cint = m .* conj(s);%		complex interferogram

% Show phase map
figure (1);
imagesc(angle(cint)); colorbar
title ('phase of complex interferogram.');

if (showperiodo==1)
  % Show range spectrum for some (total) lines
  % use periodogram matlab signal toolbox pp3-5 (?)
  plotnumber = 1;
  nfft       = 128;
  %noverlap   = 0;
  noverlap   = nfft/2;
  noverlap   = 96;
  window     = hamming(nfft);%			weight window
  %window     = hanning(nfft);%			weight window
  %window     = boxcar(nfft);%			weight window
  fs         = size(cint,2);%				scaling spectrum
  dflag      = 'none';%				do not detrend
  xaxispsd   = (0:nfft-1) / nfft*fs;%		fftshifted??
  for i=1:L
    if (rem(i,sub)==0)
      x = cint(i,:);%				total width
      figure(plotnumber);
      hold on

      % Periodogram
      Irange = psd(x,nfft,fs,window,noverlap,dflag);
      [maxvalue,ii]=max(Irange);
      SNR = maxvalue ./ ((sum(Irange)-maxvalue)/nfft);
      shifts1(plotnumber) = (ii - nfft) * P/nfft;%	?? OJA
      snrs1  (plotnumber) = SNR;
      %subplot(3,1,1), plot(xaxispsd,10*log10(Irange),'r');
      subplot(3,1,1), plot(xaxispsd,Irange,'r');
      eval(['title (''range spectrum for interferogram line ',num2str(i),''')']);

      % Simple FFT over 1 line
      %Irange   = abs(fftshift(fft(x)));
      Irange   = fft(x);
      Irange   = Irange .* conj(Irange);%		power
      [maxvalue,ii]=max(Irange);
      SNR = maxvalue ./ ((sum(Irange)-maxvalue)/P);
      shifts2(plotnumber) = ii - P;
      snrs2  (plotnumber) = SNR;
      subplot(3,1,2), plot(xaxis,Irange,'b');
      eval(['title (''range spectrum by simple fft line ',num2str(i),''')']);

      % FFT mean over 8 lines (kind of periodogram)
      Irange   = fft(cint(i-7:i,:),[],2);%		FFT over range
      Irange   = Irange .* conj(Irange);%		power
      Irange   = ((mean(Irange,1)));%			mean over blocklines,
      [maxvalue,ii]=max(Irange);
      SNR = maxvalue ./ ((sum(Irange)-maxvalue)/P);
      shifts3(plotnumber) = ii - P;
      snrs3  (plotnumber) = SNR;
      subplot(3,1,3), plot(xaxis,Irange,'b');
      eval(['title (''range spectrum by average fft 8 line ',num2str(i),''')']);

      hold off
      plotnumber = plotnumber + 1;
    end
  end
  
  disp('maxima complex interferogram spectrum');
  shifts1= shifts1
  %snrs1  = snrs1
  shifts2= shifts2
  %snrs2  = snrs2
  shifts3= shifts3
  %snrs3  = snrs3
end


%%% Do addaptive filtering based on shifts, warn if snrs low...
% blocks 8*256, compute shift, fft(m),fft(s), trunc2zero correct part, cint2
if (filtadaptive==1)
  cint2=complex(zeros(size(cint)));
  blocklines=8;
  blockrange=256;%		4000 m
  %blockrange=128;
  numblocksazi = L/blocklines;
  numblocksrange = P/blockrange;
  %starti=floor(blocklines/2);
  starti=1;
  rangesamplerate = 18.96;
  rangebandwidth = 15.5;
  % number of zeros on one side of the dataspectrum
  numzerosdata = floor(((rangesamplerate-rangebandwidth)/rangesamplerate)*(blockrange/2));
  for i=1:numblocksazi
    startj=1;
    for i=1:numblocksrange
      partcint = cint(starti:starti+blocklines-1,startj:startj+blockrange-1);
      m2       = m(starti:starti+blocklines-1,startj:startj+blockrange-1);
      s2       = s(starti:starti+blocklines-1,startj:startj+blockrange-1);
      Irange   = fft(partcint,[],2);%			FFT over range
      Irange   = Irange .* conj(Irange);%		power
      Irange   = (fftshift(mean(Irange,1)));%		mean over blocklines,
      							%%NOTE: fftshift vector!
      [maxvalue,ii]=max(Irange);
      SNR      = maxvalue ./ ((sum(Irange)-maxvalue)/blockrange);
      shift    = ii - blockrange/2;
      
      if (SNR>2)
      % Filter master and slave with shift
      % do fftshift for convenience over all rows
      M2 = (fft(m2,[],2));%			FFT over range
      S2 = (fft(s2,[],2));%			FFT over range
      for i=1:blocklines
        M2(i,:) = fftshift(M2(i,:));%		1d shifts
        S2(i,:) = fftshift(S2(i,:));%		1d shifts
      end

      % decide which side 2b filtered
      masterzeros = 1:numzerosdata+abs(shift);
      slavezeros  = blockrange-numzerosdata-abs(shift):blockrange;
      %if (shift > 0)%					swap?
      if (shift < 0)%					swap?
        tmpvar = masterzeros;
	masterzeros = slavezeros;
        slavezeros  = tmpvar;
      end
      M2(:,masterzeros)= 0;%				spectral filter
      S2(:,slavezeros) = 0;%				spectral filter
      for i=1:blocklines
        M2(i,:) = ifftshift(M2(i,:));%		1d shifts
        S2(i,:) = ifftshift(S2(i,:));%		1d shifts
      end
      m2 = ifft(M2,[],2);%                    FFT over range
      s2 = ifft(S2,[],2);%                    FFT over range
      else
	disp('SNR < 3; no filtering');
        SNR=SNR
      end

      % Filtered cint
      cint2(starti:starti+blocklines-1,startj:startj+blockrange-1)=m2.*conj(s2);
      startj   = startj+blockrange;
    end
    starti = starti+blocklines;
  end

  % Show phase map
  figure;
  imagesc(angle(cint2)); colorbar
  title ('phase of rangefiltered complex interferogram.');
  
  % Show phase map
  figure;
  imagesc(angle(cint .* conj(cint2)))
  title ('diferecne complex interferogram.');

  var1 = var(angle(cint(:)))
  var2 = var(angle(cint2(:)))
end


%%% EOF
more on
