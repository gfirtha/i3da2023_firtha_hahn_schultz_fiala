clear
close all

q = 10;
R0 = 2;
Nssd = 2400;
[x0,y0] = pol2cart(linspace(0,360*(1-1/Nssd), Nssd)'*pi/180,repmat(R0,Nssd,1));
x0 = [x0,y0];
n0 = -x0/R0;

fs = 20e3;
N = 1024;
xr = [0,0];
xs = [3,0];
omega = 2*pi*(0:N-1)'/N*fs;
H_pre = fft(get_wfs_prefilter(N,fs),N);

W0 = 1;
k0a = 2;
[amp,delay] = get_delay_and_gain(x0,n0,xs,xr);
[AA0,k0] = get_antialiasing_filters(x0(1:q:end,:),n0(1:q:end,:),xs,xr,fs,N,W0);
AA = zeros(N,Nssd);
AA(:,1:q:end) = AA0;
%D_wfs = ifft((H_pre*amp.').*exp(-1i*omega*delay'),length(omega),1,'symmetric');
D_wfs = ifft(AA.*(H_pre*amp.').*exp(-1i*omega*delay'),length(omega),1,'symmetric');
comb = zeros(size(D_wfs));
comb(:,1:q:end) = 1;
D_wfs = D_wfs.*comb;


ks = ((0:Nssd-1)-Nssd/2)/Nssd*q;
Dks = fft2(D_wfs,size(D_wfs,1),size(D_wfs,2));

ftsize = 12;
f = figure('units','normalized','outerposition',[0.2 0.3 0.35 0.4]);
set(f,'defaulttextinterpreter','latex')

pos =  [0.11,0.19,0.85,0.75];
p1 = axes('Units','normalized','Position',pos(1,:));

x = ks;
y =  omega/2/pi/fs;
pcolor( x,y, q*abs(fftshift(Dks,2)));
yl0 = [0,0.5];
ylim(yl0)
xlim([-1,1])
caxis([0,50])
shading interp

hold on
% plot([0,-k0],yl0,'-w','LineWidth',1)
% plot([-0.5/W,-0.5/W-k0],yl0,'--w','LineWidth',1)
% plot([0+0.5/W,0+0.5/W-k0],yl0,'--w','LineWidth',1)
xlabel( '$k_s/k_{s,\mathrm{s}}$ []', 'FontSize', ftsize );
ylabel( '$f / f_s$ []', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
set(gca,'XTickLabelMode','auto')
b = get(gca,'YTickLabel');
set(allAxesInFigure,'YTickLabel',b,'FontSize',ftsize);
set(gca,'YTickLabelMode','auto')


set(gcf,'PaperPositionMode','auto');
%print( '-r300', 'ideal_synthesis_spectrum' ,'-dpng')
%%


function h = get_wfs_prefilter(N,fs)
Nfilt = 128;
om = (0 : N/2 - 1)'/N*2*pi*fs;
H = sqrt(1i*om);
h = fftshift(ifft(H,N,'symmetric'));
h = h(round(end/2)-Nfilt/2+1:round(end/2)+Nfilt/2).*hann(Nfilt);
end


function [amp,delay] = get_delay_and_gain(x0,n0,xs,xr)
c = 343;
dl = sqrt(sum((x0-circshift(x0,1)).^2,2));
rho_G = sqrt(sum((x0-xr).^2,2));
A0 = 1;
kP = bsxfun( @minus, [0,0], xs );
kh = kP/norm(kP);
khn = sum(kh.*n0,2);
win = double(khn>=0);
dref = sqrt(rho_G);
delay =  (x0*kh') / c;
amp = A0*sqrt(8*pi/c).*win.*dref.*khn.*dl;
delay = delay - min(delay);
end

function [AA,k0t] = get_antialiasing_filters(x0,n0,xs,xr,fs,N,W)
Nbut = 4;
dl = sqrt(sum((x0-circshift(x0,1)).^2,2));
kP = repmat(bsxfun( @minus, [0,0], xs ),[size(x0,1),1]);
kP = kP./sqrt(sum(kP.^2,2));
v0 = [n0(:,2), -n0(:,1)];

d = sqrt(sum((xr-x0).^2,2));
[~,ix_stat] = min(abs(sum((x0 + d.*kP - xr).^2,2)));
k0t = sum(kP(ix_stat,:).*v0(ix_stat,:),2);
kht = sum( kP.*v0,2) - k0t;
[~,ixm] = min(kht);
kht(ixm) = 0;
wc = pi./dl.*343./abs(kht)/W;
w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*fs;
[Wc,W] = meshgrid(wc,w);
AA = 1./sqrt(1+(W./Wc).^(2*Nbut));

% kht = sum( kP.*v0,2);
% [~,ixm] = min(kht);
% kht(ixm) = 0;
% wc = pi./dl.*343./abs(kht);
% w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*fs;
% [Kt,W] = meshgrid(kht,w);
% AA = 1./sqrt(1+( W/343.*Kt/(pi/mean(dl))  - k0).^(2*Nbut));
end
