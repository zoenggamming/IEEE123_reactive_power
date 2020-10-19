%clc; clear all; close all;
% V=4.16e3;
% Ts=50e-6;
% load('Datos_Net_123.mat');

%% Chirp data to impulse data
%Datos modulando la entrada uno
load('zon3.mat');


%plot(Zon1_PL1.Data(:,3:4))
tt=Zon3_PL5.Time;
Vy=(Zon3_PL5.Data(:,4:5)-0);
Qu=Zon3_PL5.Data(:,1:3);

%subplot(2,1,1);plot(Qu);subplot(2,1,2);plot(Vy);
ts=tt; dys=Vy; dus=Qu;
its=1;eve=20000;
nyus=length(dys(1,:)); %%Longitud del vector de medida en D-PMU, cantidad de datos del vector
nus=length(dus(1,:)); %%Longitud del vector de señal de exitacion Chirp, cantidad de datos del vector
dt=ts(104)-ts(103);   %% Periodo de muestreo, donde ts es el vector de tiempo
fs=1/dt;              %% frecuencia de muestreo de muestreo
nt=length(ts(eve:end));%%Longitud del vector de tiempo desde el momento donde inicia la señal Chirp (eve) hasta el final(end)
ts=ts(eve:end,1);  %% se ajusta el vector de tiempo para la identificacion solo desde que se inicia la señal chirp
ux=dus; %d1(ius,:); % 
stu=ux(eve:end,1); %Se ajusta el vector desde donde inicia la señal Chirp hasta el final
sigu=fft(stu); %%Se obtiene la fft de la señal Chirp usada en la identificación
yx=dys; %d1(iys,:); 
for yyy=1:2
sty=yx(eve:end,yyy);%Se ajusta el vector de medidas desde donde inicia la señal Chirp hasta el final 
sty=sty-ones(size(sty(:,1)))*mean(sty)*1;
omega=fft(sty);%%Se obtiene la fft de la señal de medidas con los D-PMU
% subplot(2,1,1);plot(stu);subplot(2,1,2);plot(sty);
sigy=ifft(omega./sigu)';
zfl=designfilt('lowpassiir','FilterOrder',15, ...
    'HalfPowerFrequency',0.0015,'DesignMethod','butter'); %0.01 0.005
st_u1(:,yyy)=filtfilt(zfl,sigy(:,:)');
end
%% 
%Datos modulando la entrada dos
tt=Zon3_PL6.Time;
Vy=(Zon3_PL6.Data(:,4:5)-0);
Qu=Zon3_PL6.Data(:,1:3);
% subplot(2,1,1);plot(Qu);subplot(2,1,2);plot(Vy);
ts=tt; dys=Vy; dus=Qu;
its=1;eve=20000;
nyus=length(dys(1,:)); 
nus=length(dus(1,:)); %% Cambio
dt=ts(104)-ts(103);
fs=1/dt;nt=length(ts(eve:end));
ts=ts(eve:end,1);
ux=dus; %d1(ius,:); % Inputs
stu=ux(eve:end,2);
sigu=fft(stu);
yx=dys; %d1(iys,:); % Outputs
for yyy=1:2
sty=yx(eve:end,yyy);
sty=sty-ones(size(sty(:,1)))*mean(sty)*1;
omega=fft(sty);
% subplot(2,1,1);plot(stu);subplot(2,1,2);plot(sty);
sigy=ifft(omega./sigu)';
st_u2(:,yyy)=filtfilt(zfl,sigy(:,:)');
end
%% 
%Datos modulando la entrada tres
tt=Zon3_PL7.Time;
Vy=(Zon3_PL7.Data(:,4:5)-0);
Qu=Zon3_PL7.Data(:,1:3);
% subplot(2,1,1);plot(Qu);subplot(2,1,2);plot(Vy);
ts=tt; dys=Vy; dus=Qu;
its=1;eve=20000;
nyus=length(dys(1,:)); 
nus=length(dus(1,:)); %% Cambio
dt=ts(104)-ts(103);
fs=1/dt;nt=length(ts(eve:end));
ts=ts(eve:end,1);
ux=dus; %d1(ius,:); % Inputs
stu=ux(eve:end,3);
sigu=fft(stu);
yx=dys; %d1(iys,:); % Outputs
for yyy=1:2
sty=yx(eve:end,yyy);
sty=sty-ones(size(sty(:,1)))*mean(sty)*1;
omega=fft(sty);
% subplot(2,1,1);plot(stu);subplot(2,1,2);plot(sty);
sigy=ifft(omega./sigu)';
st_u3(:,yyy)=filtfilt(zfl,sigy(:,:)');
end
%% Identification
% st(2,:)=filtfilt(zfl,sigy(2,:));
figure(3);plot(ts,sigy,'k');hold on;plot(ts,st_u1);plot(ts,st_u2);

t=ts'; yg123=[st_u1 st_u2 st_u3];
% figure(6);plot(t,yg12);hold on;
t=ts; yimp=yg123;
it=1:8010; nus=3; nys=2;
t=downsample(t(it,:),10);
sm1=downsample(yimp(it,:),10);
figure(3);plot(t,sm1,'y--')

means=mean(sm1);
sm1=sm1-ones(size(sm1(:,1)))*means*1;
[N,nyus]=size(sm1);

r=400;
c1=0; c2=1;

for m=1:r
    for l=1:r
H0(1+c1*nyus/nus:c2*nyus/nus,1+(l-1)*nus:l*nus)=reshape(sm1(l+c2,:)',nyus/nus,nus);
H1(1+c1*nyus/nus:c2*nyus/nus,1+(l-1)*nus:l*nus)=reshape(sm1(l+c2+1,:)',nyus/nus,nus);
    end
    c1=c1+1;
    c2=c2+1;
end
[U,S,V] = svd(H0);
energ=cumsum(diag(S))./sum(diag(S));
poenerg=find(energ>0.99); nx=poenerg(1); % nx=nx; 
Sn=S(1:nx,1:nx);
Un=U(:,1:nx);
Vt=V';Vn=Vt(1:nx,:);
Ad=(Sn^(-1/2))*Un'*H1*Vn'*(Sn^(-1/2));
Bd=Sn^(1/2)*Vn(1:nx,1:nus);
Cd=Un(1:nys,1:nx)*Sn^(1/2);
sysdisc=ss(Ad,Bd,Cd,0,dt);
syscont=d2c(sysdisc,'zoh');
[a_mat2,b_vr2,c_spd2,d_e2]=ssdata(syscont);
eia=eig(a_mat2); par=[imag(eia)/(2*pi) real(eia) -100*real(eia)./abs(eia)]
TT=0:0.001:0.4;
YY=impulse(ss(a_mat2,b_vr2(:,1),c_spd2,0),TT);
YY=[YY impulse(ss(a_mat2,b_vr2(:,1),c_spd2,0),TT)];
% figure(6);plot(TT,YY*dt+ones(size(YY(:,1)))*means,'k--')
sysdisc2 = c2d(syscont,Ts,'zoh');
[Ad2_3,Bd2_3,Cd2_3,Dd2_3]=ssdata(sysdisc2);
Qd2=eye(size(Ad2_3)); %Cc'*Cc
%Qd2=diag(Qd2); Qd2(1:4)=0.0001; Qd2(7)=1; Qd2([5 6])=20; Qd2=diag(Qd2);
Rd1=eye(3)*0.01;Rd2=eye(2)*0.01;
Kd2_3 = dlqr(Ad2_3,Bd2_3,Qd2,Rd1)'; %(:,1)
Gd2_3 = dlqr(Ad2_3',Cd2_3',Qd2,Rd2)'; %(1,:)



