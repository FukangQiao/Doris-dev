% Script fracdemo.m
%
% Script to demonstrate some characteristic fractal/power-law  lines and
% surfaces 
%
% Ramon Hanssen, May 2000

doplot = 'n';
Nline  = 128;
Nsurf  = 128;
betas1  = [0,  1, 2];
%ylab1   = strvcat('White noise','1/f noise','fBm');
ylab1   = strvcat('Witte ruis','1/f ruis','fBm');
titl1   = strvcat('\beta=0','\beta=1','\beta=2');
betas2  = [2/3,   5/3,   8/3];
ylab2   = strvcat('2/3','5/3','8/3');
titl2   = strvcat('\beta=2/3','\beta=5/3','\beta=8/3');

colormap gray
linecol = 'k'

for j = 1:2,
   beta = eval(['betas',num2str(j)]);
   for i=1:3
     fracdim1 = (5 - beta(i))/2;
     fracdim2 = (7 - beta(i))/2;
     [flineb0] = fracline(Nline,beta(i),doplot);
     [fsurfb0] = fracsurf(Nsurf,beta(i),doplot);   
 
     k=[1:100] ; fpowb0 = exp( -beta(i) .* log(k) ); 


     figure(j);
     subplot(3,3,((i-1)*3)+1);plot(flineb0,linecol);
        set(gca,'xlim',[1 Nline],'yticklabel',[],'xticklabel',[]);
        ylabel(deblank( eval(['ylab',num2str(j),' (i,:)'])));
        title(['Frac. Dim. ',num2str(fracdim1)])
     subplot(3,3,((i-1)*3)+2);imagesc(fsurfb0);
        set(gca,'yticklabel',[],'xticklabel',[]);
        title(['Frac. Dim. ',num2str(fracdim2)])
     subplot(3,3,((i-1)*3)+3);loglog(k,fpowb0,linecol);
        set(gca,'yticklabel',[],'xticklabel',[]);
        text(3,0.2,deblank(eval(['titl',num2str(j),' (i,:)'])));
        title(['1D energie spectrum'])
   end
   % adjust yaxis slopes
   subplot(3,3,3);yl1 = get(gca,'ylim');
   subplot(3,3,6);yl2 = get(gca,'ylim');
   subplot(3,3,9);yl3 = get(gca,'ylim');
    ylmi = min([min(yl1),min(yl2),min(yl3)]);
    ylma = max([max(yl1),max(yl2),max(yl3)]);
   subplot(3,3,3);set(gca,'ylim',[ylmi ylma]);
   subplot(3,3,6);set(gca,'ylim',[ylmi ylma]);
   subplot(3,3,9);set(gca,'ylim',[ylmi ylma]);
end


figure(1);print fracdemo1.eps -depsc
figure(2);print fracdemo2.eps -depsc

