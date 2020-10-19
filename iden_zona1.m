%clc; clear all; close all;
V=4.16e3;
Ts=50e-6;
load('Datos_Net_123.mat');

%% Chirp data to impulse data
%Datos modulando la entrada uno
load('zon1.mat');
tt=Zon1_PL1.Time;

%plot(Zon1_PL1.Data(:,3:4))

Vy=(Zon1_PL1.Data(:,3:4)-0);
Qu=Zon1_PL1.Data(:,1:2);

%subplot(2,1,1);plot(Qu);subplot(2,1,2);plot(Vy);
ts=tt; dys=Vy; dus=Qu;
its=1;eve=20000;
nyus=length(dys(1,:)); 
nus=length(dus(1,:)); %% Cambio
dt=ts(104)-ts(103);
fs=1/dt;nt=length(ts(eve:end));
ts=ts(eve:end,1);
ux=dus; %d1(ius,:); % Inputs
stu=ux(eve:end,1);
sigu=fft(stu);
yx=dys; %d1(iys,:); % Outputs
sty=yx(eve:end,:);
sty=sty-ones(size(sty(:,1)))*mean(sty)*1;
omega=fft(sty);
% subplot(2,1,1);plot(stu);subplot(2,1,2);plot(sty);
sigy=ifft(omega./sigu)';
zfl=designfilt('lowpassiir','FilterOrder',15, ...
    'HalfPowerFrequency',0.0015,'DesignMethod','butter'); %0.01 0.005
st_u1=filtfilt(zfl,sigy(:,:)');

%% 
%Datos modulando la entrada dos
tt=Zon1_PL2.Time;
Vy=(Zon1_PL2.Data(:,3:4)-0);
Qu=Zon1_PL2.Data(:,1:2);
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
sty=yx(eve:end,:);
sty=sty-ones(size(sty(:,1)))*mean(sty)*1;
omega=fft(sty);
% subplot(2,1,1);plot(stu);subplot(2,1,2);plot(sty);
sigy=ifft(omega./sigu)';
st_u2=filtfilt(zfl,sigy(:,:)');

%% Identification
% st(2,:)=filtfilt(zfl,sigy(2,:));
figure(3);plot(ts,sigy,'k');hold on;plot(ts,st_u1);plot(ts,st_u2);

t=ts'; yg12=[st_u1 st_u2];
% figure(6);plot(t,yg12);hold on;
t=ts; yimp=yg12;
it=1:8010; nus=2; nys=2;
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
[Ad2_1,Bd2_1,Cd2_1,Dd2_1]=ssdata(sysdisc2);
Qd2=eye(size(Ad2_1))*1; %Cc'*Cc
Qd2=diag(Qd2); Qd2(1:4)=0.0001; Qd2(7)=1; Qd2([5 6])=20; Qd2=diag(Qd2);
Rd2=eye(1)*0.01;
Kd2_1 = dlqr(Ad2_1,Bd2_1,Qd2,Rd2)'; %(:,1)
Gd2_1 = dlqr((Ad2_1)',(Cd2_1)',Qd2,Rd2)'; %(1,:)

% Ts_Power = 1.0000e-05;
% Ts_Control = 1.0000e-05;
% 
% Tt=V_1.Time;
% 
% figure
% subplot(3,1,1)
% plot(Tt,V_1.Data,Tt,V_2.Data,Tt,V_3.Data,Tt,V_4.Data,Tt,V_5.Data,Tt,V_6.Data,Tt,V_7.Data,Tt,V_8.Data,Tt,V_9.Data,Tt,V_10.Data,Tt,V_11.Data,Tt,V_12.Data,Tt,V_13.Data)
% subplot(3,1,2)
% plot(Tt,PQ_1.Data(:,2),Tt,PQ_2.Data(:,2),Tt,PQ_3.Data(:,2),Tt,PQ_4.Data(:,2),Tt,PQ_5.Data(:,2),Tt,PQ_6.Data(:,2),Tt,PQ_7.Data(:,2),Tt,PQ_8.Data(:,2),Tt,PQ_9.Data(:,2),Tt,PQ_10.Data(:,2),Tt,PQ_11.Data(:,2),Tt,PQ_12.Data(:,2),Tt,PQ_13.Data(:,2))
% subplot(3,1,3)
% plot(Tt,PQ_1.Data(:,1),Tt,PQ_2.Data(:,1),Tt,PQ_3.Data(:,1),Tt,PQ_4.Data(:,1),Tt,PQ_5.Data(:,1),Tt,PQ_6.Data(:,1),Tt,PQ_7.Data(:,1),Tt,PQ_8.Data(:,1),Tt,PQ_9.Data(:,1),Tt,PQ_10.Data(:,1),Tt,PQ_11.Data(:,1),Tt,PQ_12.Data(:,1),Tt,PQ_13.Data(:,1))
% 
% figure
% plot(Tt,V_8.Data,Tt,V_7.Data,Tt,V_9.Data,Tt,V_12.Data)
% 
% size(V_8.Data)
% 
% for i = 25000:+1:32000
%    V_8.Data(i)=V_7.Data(i);
%    V_12.Data(i)=V_7.Data(i);
% end
% 
% V_8.Data(25000)=555555;


