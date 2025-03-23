% FILTRANGE
% test parameters to use for adaptive filtering in range
% how nfft, mean over how many lines, snr < x do what?, etc.
%

% $Revision: 1.6 $  $Date: 2000/04/10 10:45:52 $
% Bert Kampes, 24-Mar-2000

%--------1---------2---------3---------4---------5---------6---------7---------8---------9
%%% Initialize
clear; close all force; more off;

%%% Initial main parameters of images.
imagedir ='/home/fmr/d1/scripts/matlab/insar/';
master   = 'master2_128_512';%          SHORT
slave    = 'slave2_128_512';%           FLOAT
%master   = 'master3_128_512';%          FLOAT
%slave    = 'slave3_128_512';%           FLOAT
L        = 128;%                        number of lines
P        = 512;
L        = 50;%                        number of lines

showspectrummaster=0;
showspectrumslave=0;
showspectrumcint=0;
filtadaptive = 1;
simulatedata = 1;
addnoise     = 0.;%			percentage of signal
addnoise     = 1.0;%			percentage of signal 1 means 50%

%%% parameters for filtblock.m
% these may be played with
nlmean       = 15;%		use mean per 9 lines
snrthreshold = 3;%		threshold on SNR for peak estimation
rsr          = 18.96;%		rangesamplerate ERS1 (leader) 
rbw          = 15.55;%		rangebandwidth ERS1 (leader) 
dooversample = 1;%            	oversample in range before interferogram gen.
dohamming    = 0;%            	apply hamming 2 times
%fftlength   = 64;%		size of partmaster
fftlength    = P;%		size of partmaster
%fftlength   = P/2;%		size of partmaster
%fftlength   = P/4;%		size of partmaster
debug = 3;
%debug = 0;


if (simulatedata)
  % properties of SAR image in range, see paper hanssen interpolation e.
  disp('_');
  disp('simulating data...');
  P4 = 4*P;
  % force odd number of nonzeros elements
  numnonzeros = ceil(P*(rbw/rsr));
  if (rem(numnonzeros,2)==0) numnonzeros = numnonzeros-1; end;
  %numzeros    = P-numnonzeros;
  offset = (numnonzeros-1)/2;

  % simulate data rule of paper RH/RB
  phase = 2*pi*rand(1,P4);
  ampli = sqrt(-log(rand(1,P4)));
  data  = ampli .* exp(i*phase);
  DATA  = fftshift(fft(data));

  % shift spectrum of slave wrt. master spectrum
  %shift = 0;%			does not look ok, but caused by simul.
  shift = round(P/2+P/6);%	dooversample is obli.
  %shift = round(P/2-P/6);%	dooversample is not obli. but yields higher snr
  shift = 20;
  %shift = 100;
  shift = 200;
  %shift = 300;
  %shift = 360;
  disp(['shift spectrum master w.r.t. slave: ' num2str(shift)]);
  M = DATA(P4/2-P/2:P4/2+P/2-1);
  S = DATA(P4/2-P/2+shift:P4/2+P/2-1+shift);

  disp(['adding gaussian noise (in percent of mean signal power): ', ...
	 num2str(100.*(1+addnoise)/2)]);
  M = M + addnoise*mean(abs(M))*complex(randn(size(M)),randn(size(M))); 
  S = S + addnoise*mean(abs(S))*complex(randn(size(S)),randn(size(S))); 

%  % hamming filtering
% BK 30-Mar-2000
% one might better perhaps dehamming before peak estimation, because
% for large shifts the power goes down of correlation.
% another approach might be to correct the correlation for the number of points
% used in the computation, so divide 0th freq. by N (2N for oversampling),
% 1st freq. by N-1, next by N-2, etc.
% or both things have to be done for best estimate of peak.
%
  deltaf   = rsr/length(M);%                              interval
  freqaxis = -rsr/2:deltaf:rsr/2-deltaf;
  if (dohamming==0)
    myfilter = myrect(freqaxis./rbw);
    disp('rect filter...');
  else
    disp('hamming filter...');
    myfilter = myhamming(freqaxis,rbw,rsr,0.75);
  end
  M = M.*myfilter;
  S = S.*myfilter;

  % signals in space domain
  m = ifft(ifftshift(M));
  s = ifft(ifftshift(S));

  % copy downwards to simulate more lines for filterblock function
  m = ones(L,1)*m;
  s = ones(L,1)*s;
%  clear data phase ampli M S offset rsr rbw P4 numnonzeros;

  %%% Plot this
  if (debug>=1)
    disp('Plotting spectra...');
    figure(10);
    subplot(3,1,1), plot(abs(fftshift(fft(m(1,:)))));
    title('master range spectrum (fftshifted)');
    subplot(3,1,2), plot(abs(fftshift(fft(s(1,:)))));
    title('slave range spectrum (fftshifted)');
    subplot(3,1,3), plot(angle(m(1,:).*conj(s(1,:))));
    title('phase of complex interferogram (not oversampled).');
  end;

else
%%% Read images from disk if not present.
  if (~exist('m','var')) m=freadbk([imagedir, master],L,'cmplxint16'); end
  %if (~exist('m','var')) m=freadbk([imagedir, master],L,'cmplxfloat32'); end
  if (~exist('s','var')) s=freadbk([imagedir, slave], L,'cmplxfloat32'); end
end


% plot variables
sub = 10;%			do a plot per sub
subplots = floor(P/sub);%	number of subplots
xaxis = 0:P-1;
xaxis = xaxis - P/2;

if (showspectrummaster==1)
  %%% Show spectrum in range for master
  %figure(1);
  total = zeros(1,P);
  %xaxis = fftshift(xaxis);
  plotnumber = 0;
  for i=1:L
    Mrange = fft(m(i,:));
    Mrange = abs(fftshift(Mrange));
    if (rem(i,sub)==0)
      plotnumber = plotnumber + 1;
      figure(plotnumber);
      %subplot(subplots,1,plotnumber),
      plot(xaxis,Mrange);
      eval(['title (''range spectrum for master line ',num2str(i),''')']);
    end
    total = total + Mrange;
  end
  %%% Averaged spectrum
  figure(plotnumber+1);
  total = total / L;
  plot(xaxis,total);
  eval(['title (''average range spectrum for all ',num2str(L),' lines (master)'')']);

  %%% Plot for all 
  figure(plotnumber+2);
  imagesc(abs(fftshift(fft(m,[],2)))); colorbar
  eval(['title (''range spectrum for all ',num2str(L),' lines (master)'')']);
end


if (showspectrumslave==1)
  %%% Show spectrum in range for slave
  disp('_');
  disp ('press key to continue with slave');
  pause;
  close all force;
  total = zeros(1,P);
  plotnumber = 0;
  for i=1:L
    Srange = fft(s(i,:));
    Srange = abs(fftshift(Srange));
    if (rem(i,sub)==0)
      plotnumber = plotnumber + 1;
      figure(plotnumber);
      plot(xaxis,Srange);
      eval(['title (''range spectrum for slave line ',num2str(i),''')']);
    end
    total = total + Srange;
  end
  %%% Averaged spectrum
  figure(plotnumber+1);
  total = total / L;
  plot(xaxis,total);
  eval(['title (''average range spectrum for all ',num2str(L),' lines (slave)'')']);
  
  %%% Plot for all 
  figure(plotnumber+2);
  imagesc(abs(fftshift(fft(m,[],2)))); colorbar
  eval(['title (''range spectrum for all ',num2str(L),' lines (slave)'')']);
end 


%%% Generate interferogram by pointwize multiplication.
cint = m .* conj(s);%		complex interferogram


if (showspectrumcint==1)
  disp('_');
  disp ('press key to continue with complex interferogram');
  pause;
  close all force;
  % Show phase map
  figure (1);
  imagesc(angle(cint)); colorbar
  title ('phase of complex interferogram.');
  % Show range spectrum for some (total) lines
  % use mean over small block in azimuth (8 lines?)
  plotnumber = 1;
  for i=1:L
    if (rem(i,sub)==0)
      Irange = fft(cint(i-7:i,:),[],2);
      %Irange = fft(cint(i-15:i,:),[],2);
      %Irange = abs(fftshift(mean(Irange,1)));
      Irange = (fftshift(mean(Irange,1)));
      Irange = Irange .* conj(Irange);
      [maxvalue,ii]=max(Irange);
      SNR = maxvalue ./ ((sum(Irange)-maxvalue)/P);
      shifts(plotnumber) = ii - P/2;
      snrs(plotnumber) = SNR;
      plotnumber = plotnumber + 1;
      figure(plotnumber);
      subplot(2,1,1), plot(xaxis,Irange,'r');
      eval(['title (''range spectrum for interferogram line ',num2str(i),''')']);
      %Irange = fft(cint(i,:));
      %Irange = abs(fftshift(Irange));
      Irange = fft(cint(i,:));
      Irange = abs(fftshift(Irange));
      hold on
      subplot(2,1,2), plot(xaxis,Irange,'b');
      hold off
    end
  end
  
  disp('maxima complex interferogram spectrum');
  shifts=shifts
  snrs=snrs
end


%%% Do addaptive filtering based on shifts, warn if snrs low...
% blocks 8*256, compute shift, fft(m),fft(s), trunc2zero correct part, cint2

if (filtadaptive==1)
  disp('_');
  disp ('press key to continue with adaptive filter');
  pause;
  close all force;

  morig = m;
  sorig = s;


  % derived pm
  blocklines      = 128;%	lines per block in filtblock routine, as large as possible
  blocklines      = L;%		lines per block in filtblock routine, as large as possible
%  returnedlines   = blocklines - nlmean + 1;

  % do the filtering
  startlinethisblock = 1;%				start at the start...
  endlinethisblock   = startlinethisblock + blocklines - 1;

  %for aziblocks = 1:blocklines:L+blocklines-1
  allazimuthdone = 0;%						
  numrangeblocks = floor(P/fftlength);%			might be last smaller one
  lastrangeblock = 0;
  %if (numrangeblocks~=ceil(P/fftlength)) lastrangeblock = 1; end;
  if (numrangeblocks*fftlength~=P) lastrangeblock = 1; end;

  for aziblocks=1:100000%					for ever
    % update index pm per new block of blocklines azimuth lines
    startpixel = 1;
    endpixel   = startpixel+fftlength-1;
    % (read from file per azimuth blocks in doris)
    if (endlinethisblock > L - floor(blocklines/4))%		last block
      allazimuthdone = 1;
      endlinethisblock = L;%					little larger
    end;
    firstreturnedline  = startlinethisblock + (nlmean-1)/2;
    lastreturnedline   = endlinethisblock   - (nlmean-1)/2;
    partmaster2 = m(startlinethisblock:endlinethisblock,:);
    partslave2  = s(startlinethisblock:endlinethisblock,:);

    %for rangeblocks = 1:fftlength:P+fftlength-1
    for rangeblocks = 1:numrangeblocks+lastrangeblock%		last=0/1
      % offer this in rangeblocks to filtblock routine, check end
      %if (rangeblocks>=P)
      if (rangeblocks==numrangeblocks+1)%			lastblock
	endpixel   = P;
	startpixel = endpixel - fftlength + 1;%			partly filtered already
      end;
      partmaster = partmaster2(:,startpixel:endpixel);
      partslave  = partslave2 (:,startpixel:endpixel);

      % call filter, returnes smaller part
      [mf,sf] = filterblock(partmaster,partslave,nlmean, ...
		    snrthreshold,rsr,rbw,dooversample,dohamming,debug);
      %[mf,sf] = filterblock(partslave,partmaster,nlmean, ...
      %		    snrthreshold,rsr,rbw,dooversample,dohamming,debug);
      m(firstreturnedline:lastreturnedline, startpixel:endpixel) = mf;
      s(firstreturnedline:lastreturnedline, startpixel:endpixel) = sf;

      % update index pm.
      startpixel = startpixel + fftlength;%		next range block
      endpixel	 = startpixel + fftlength-1;
    end;
  % update index pm.
  startlinethisblock = lastreturnedline + 1 - (nlmean-1)/2;%	centered round new frstretrnd
  endlinethisblock   = startlinethisblock + blocklines - 1;
  if (allazimuthdone==1) break; end;
  end;

  %Figures
  cintorig = morig.*conj(sorig);
  cint= m.*conj(s);
  figure;
    imagesc(angle(cintorig)); colorbar; title('orig');
  figure;
    imagesc(angle(cint)); colorbar; title('rfiltered');
  % show spectra?
  % compute snr?
  %var(angle(cintorig(:)))
  %var(angle(cint(:)))
  %residues(cintorig)
  %residues(cint)

  %coh1 =  abs((morig(100,:)*sorig(100,:)')) / ...
%	 sqrt((morig(100,:)*morig(100,:)') * (sorig(100,:)*sorig(100,:)'))
%  coh2 =  abs((m(100,:)*s(100,:)')) / ...
%	 sqrt((m(100,:)*m(100,:)') * (s(100,:)*s(100,:)'))
%
%  coh2 =  abs(sum(cint(100,:))) / ...
%	 sqrt(sum(m(100,:).^2) * sum(s(100,:).^2))

  line = floor(L/2);
  figure;
  subplot(3,1,1), plot(angle(cintorig(line,:)));
  hold on; plot(angle(cintorig(line,:)),'r+');
  subplot(3,1,2), plot(angle(cint(line,:)));
  hold on; plot(angle(cint(line,:)),'r+');
  subplot(3,1,3), plot(angle(cintorig(line,:)) - angle(cint(line,:)));
end




%%% EOF
more on
