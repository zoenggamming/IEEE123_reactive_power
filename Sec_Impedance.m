clc; clear all; close all;

Ts_Powerv2 = 1.0000e-05;

Ts_Power = 5.0000e-05;
Ts_Control = 5.0000e-05;

parameters = xlsread('lines_parameters.xlsx','digsilent');
%load('lines_parameters.mat');

Ainv=(1/3)*[1, 1, 1; 1, (-0.5+1i*0.866), (-0.5-1i*0.866);...
    1, (-0.5-1i*0.866), (-0.5+1i*0.866)];
A=inv(Ainv);

for h = 1:11
   Z0= parameters(h,4)+1i*parameters(h,5); %Z0=R0+1i*X0;
   Z1= parameters(h,2)+1i*parameters(h,3);%Z1=R1+1i*X1;
   Z2= Z1;%Z2=Z1
   
   Z012 =[Z0 0 0; 0 Z1 0; 0 0 Z2];
   
   line_name = strcat('line',num2str(h));
   Zabc.(line_name)= A*Z012*Ainv; %impedance sequence matriz
   line_name2 = strcat('line',num2str(h));
   %All dates of Zabc diagonal are equals, so all R and L in a,b,c phases are iquals to. I will consider only Zabc(1,1) to obtain R and L for each line. 
   R.(line_name2)=real(Zabc.(line_name)(1,1));
   L.(line_name2)=(imag(Zabc.(line_name)(1,1)))/(2*pi*60);
end

% figure (1)
% subplot(2,1,1);
% plot(vpu);hold on;
% plot(vpu_closed_loop);
% xlabel('Time (s)');
% ylabel('Voltage (PU)');
% lgnd=legend({'Voltage open loop','Voltage closed loop'});
% set(0,'defaultLineLineWidth', 4);
% set(0,'defaultAxesFontName', 'Times New Roman');
% set(0,'DefaultAxesFontSize',14);
% grid on;

% subplot(2,1,2);
% plot(Q);hold on;
% xlabel('Time (s)');
% ylabel('Q (kVA');
% lgnd=legend({'Reactive power injected into the grid'});
% set(0,'defaultLineLineWidth', 4);
% set(0,'defaultAxesFontName', 'Times New Roman');
% set(0,'DefaultAxesFontSize',14);
% grid on;



% Control
load('Caso_03_t5em5_ux10.mat');
tt=QV_ident.Time;
Vy=(QV_ident.Data(:,1)-0)*1e1;
Qu=QV_ident.Data(:,2);

% tt=Con_LQG.Time;
% Vy=(Con_LQG.Data(:,1)-0)*1e1;
% Qu=Con_LQG.Data(:,2);


% subplot(2,1,1);plot(Qu);subplot(2,1,2);plot(Vy);
ts=tt; dys=Vy; dus=Qu;
its=1;eve=27521;
nys=length(dys(1,:));
nus=length(dus(1,:)); %% Cambio
dt=ts(104)-ts(103);
fs=1/dt;nt=length(ts(eve:end));
ts=ts(eve:end,1);
ux=dus; %d1(ius,:); % Inputs
stu=ux(eve:end,1);
sigu=fft(stu);
yx=dys; %d1(iys,:); % Outputs
sty=yx(eve:end,1);
sty=sty-mean(sty);
omega=fft(sty);
% subplot(2,1,1);plot(stu);subplot(2,1,2);plot(sty);
sigy=ifft(omega./sigu)';
zfl=designfilt('lowpassiir','FilterOrder',15, ...
    'HalfPowerFrequency',0.0015,'DesignMethod','butter'); %0.01 0.005
st(1,:)=filtfilt(zfl,sigy(1,:));
figure(3);plot(ts,sigy,'k'); hold on; plot(ts,st);
t=ts'; yg1=st';
% figure;plot(t,yg1)
t=ts; yimp=yg1;
it=1:8010; nus=1; nys=1;
t1=downsample(t(it,:),10);
st1=downsample(yimp(it,:),10);
plot(t1,st1)
[N,ns]=size(st1);
dt=t1(round(ns/2,0)+1)-t1(round(ns/2,0)); 
fs=1/dt;means=mean(st1);
for k=1:ns
    st1(:,k)=st1(:,k)-means(k);
end
st=st1;
r=400;
c1=0; c2=1;
for m=1:r
    for l=1:r
        H0(1+c1*ns:c2*ns,l)=st(l+c2,:)';
        H1(1+c1*ns:c2*ns,l)=st(l+c2+1,:)';
    end
    c1=c1+1;
    c2=c2+1;
end
[U,S,V] = svd(H0);
energ=cumsum(diag(S))./sum(diag(S));
poenerg=find(energ>0.999); nx=poenerg(1); % nx=nx; 
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
% figure(7);subplot(1,2,2);plot(TT,(impulse(ss(a_mat2,b_vr2,c_spd2,0),TT))*dt+means)
sysdisc2 = c2d(syscont,Ts_Power,'zoh');
[Ad2,Bd2,Cd2,Dd2]=ssdata(sysdisc2);
Qd2=eye(size(Ad2))*1; %Cc'*Cc
% Qd2=diag(Qd2); Qd2(1:4)=0.01 ; Qd2=diag(Qd2)
Rd2=eye(1)*1;
Kd2 = dlqr(Ad2,Bd2,Qd2,Rd2)'; %(:,1)
Gd2 = dlqr(Ad2',Cd2',Qd2,Rd2)'; %(1,:)

% figure (1999);
% plot(vpu_closed_loop_pi);hold on;
% plot(vpu_closed_loop_lqg);hold on;
% plot(vpu_open_loop);
% xlabel('Time (s)', 'FontSize',12, 'FontName', 'Times New Roman');
% ylabel('V (PU)', 'FontSize',12, 'FontName', 'Times New Roman');
% lgnd=legend({'PI','LQG','LQG+INT', 'Open-loop'})
% set(lgnd,'FontSize',12);
% set(lgnd,'FontName','Times New Roman');
% set(gca,'FontSize',12,'FontName','Times New Roman');

delay_length1=0;
delay_length2=0;  
i=1;

sim("unbal_13nodes_case1_wlatency.slx")
V_out_test(:,i)=vpu_closed_loop_pi.Data;
Q_out_test(:,i)=Q_node5.Data;
P_out_test(:,i)=P_node5.Data;

tt=vpu_open_loop.Time;

% figure
% subplot(2,1,1)
% plot(tt,V_out_test(:,1),tt,V_out_test(:,2),tt,V_out_test(:,3));
% subplot(2,1,2)
% plot(tt,Q_out_test(:,1),tt,Q_out_test(:,2),tt,Q_out_test(:,3));
% 
% 
% 
% figure
% ttt=tt(1:20001);
% subplot(2,1,1)
% plot(ttt,V_out_test(10000:30000,1),ttt,V_out_test(10000:30000,2),ttt,V_out_test(10000:30000,3));
% subplot(2,1,2)
% plot(ttt,Q_out_test(10000:30000,1),ttt,Q_out_test(10000:30000,2),ttt,Q_out_test(10000:30000,3));
% 
% figure
% plot(tt, Q_out_test(:,2))
% 
% figure 
% plot(tt, V_out_test(:,2))
