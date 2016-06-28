clc
clear
close all

data = load('./DomainIntg_rhoRMS.dat');
ref = load('./rhoPrime-LES-32Cube-Spyropoulos.dat');

fac = 8 * pi^3;
plot(data(:,1)/1.809, (data(:,3)./fac).^0.5,'-+');
hold on
plot(ref(:,1), ref(:,2), 'o');
ylim([0 0.20])

%
figure

%x = 0.5*(data(1:end-1,1) + data(2:end,1))/1.81;
%tke_dt = diff(data(:,4));
plot(data(:,1)/1.81, data(:,4)./data(1,4),'-+');
%ref = load('./TKE-LES-32Cube_Rizzetta.dat');
%hold on

kim = load('./turbQuan_NoSGS_wFilterPade.dat');
kim1 = load('./turbQuan_NoSGS_wFilterSFo11p.dat');

%kim_x = 0.5*(kim(1:end-1,1) + kim(2:end,1)) / 1.81;
%kim_tke_dt = diff(kim(:,2));
hold on
plot(kim(:,1)/1.81, kim(:,2)./kim(1,2), 'o');
plot(kim1(:,1)/1.81, kim1(:,2)./kim1(1,2), 'ro');

%{
figure
%kinetic energy spectrum 
spec = load('./spectrum/out.dat');
nmode = 1300;
Y = fft(spec(:,2), nmode);
Pyy = Y.*conj(Y) / nmode;
f = 1000/nmode*(0:127);
xx = 1:1:nmode;
loglog(xx, Pyy,'-o')
%xlim([0, 15])

%semilogy(abs(fft(spec(:,2))));
%xlim([0, 20])
%}


