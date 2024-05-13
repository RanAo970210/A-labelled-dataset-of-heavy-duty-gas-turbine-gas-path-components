% Dataalgorithm
% Project: A-labelled-dataset-of-heavy-duty-gas-turbine-gas-path-components

%% Real unit operating data
%%
%%
%% import data
clc
clear
load Realunitoperatingdata.mat;

data_all_zz=Realunitoperatingdata;
% Hour level data
data_all=data_all_zz([1:60:length(data_all_zz)],:);
Date=data_all(:,1);
% unit conversion
P2_zs=1000.*data_all(:,2);%KPa-Pa
T2_zs=data_all(:,3)+273;%℃-K
T3_zs=data_all(:,4)+273;%℃-K
T4_zs=data_all(:,5)+273;%℃-K
P_zh=data_all(:,end);
DATA=[Date,P2_zs,T2_zs,T3_zs,T4_zs,P_zh];
%% Data processing
% Steady-state screening
%Non-start-stop process
k=find(DATA(:,end)>=0.7);
gonglv_fqtj=DATA(k,end);
DATA_fqtj=DATA(k,:);

% historical data steady-state screening 
% sliding window length = 10 ，significance level a = 0.05 characterized by operating conditions 
% t test ( t-test ) critical value table ( critical confidence level ) 
% The sliding window length is set to the difference of 10，working conditions as the feature. 
gonglvchazhi=gonglv_fqtj(2:length(gonglv_fqtj))-gonglv_fqtj(1:length(gonglv_fqtj)-1);
for ii = 1 : length(gonglvchazhi)-10
   a = gonglvchazhi(ii : ii + 9);
   aa(:,ii)=a;%aa is the sliding window data of power difference.
end
%t_0.025(n-1)=t_0.025(9)=2.262
for i=1:length(aa)
gonglv_pj=mean(aa(:,i));%mean value
gonglv_s=sqrt(sum((aa(:,i)-gonglv_pj).^2)/(10-1));%S
%confidence interval
zxqj_z=gonglv_pj-2.262*gonglv_s/sqrt(10);
zxqj_y=gonglv_pj+2.262*gonglv_s/sqrt(10);

m1(i)=zxqj_z;m2(i)=zxqj_y;
end
m=[m1;m2];%confidence interval
zxqjsy=find(m1<=0&m2>0);%Index of steady-state sliding window

DATA_zxqjsy=DATA_fqtj(zxqjsy,:);%DATA_zxqjsy is steady-state data
%plot(k(zxqjsy),DATA_zxqjsy(:,3),'b.')

%% The health value model considers the data of the first 30 hours as the health state, and constructs the health value model.
niheP2=polyfit(DATA_zxqjsy(1:30,end),DATA_zxqjsy(1:30,2),1);
fP2=polyval(niheP2,DATA_zxqjsy(:,end));

niheT2=polyfit(DATA_zxqjsy(1:30,end),DATA_zxqjsy(1:30,3),1);
fT2=polyval(niheT2,DATA_zxqjsy(:,end));

niheT3=polyfit(DATA_zxqjsy(1:30,end),DATA_zxqjsy(1:30,4),1);
fT3=polyval(niheT3,DATA_zxqjsy(:,end));

niheT4=polyfit(DATA_zxqjsy(1:30,end),DATA_zxqjsy(1:30,5),1);
fT4=polyval(niheT4,DATA_zxqjsy(:,end));

% Deviation
Delta_P2=100.*(DATA_zxqjsy(:,2)-fP2)./fP2;
Delta_T2=100.*(DATA_zxqjsy(:,3)-fT2)./fT2;
Delta_T3=100.*(DATA_zxqjsy(:,4)-fT3)./fT3;
Delta_T4=100.*(DATA_zxqjsy(:,5)-fT4)./fT4;

Delta_sc=[Delta_P2,Delta_T2,Delta_T3,Delta_T4,DATA_zxqjsy(:,end)];%The deviation of real unit operation data
Delta_sc_bq=DATA_zxqjsy(:,1);%Date label of deviation
save('Delta_sc.mat','Delta_sc');
save('Delta_sc_bq.mat','Delta_sc_bq');
%% % BPtraining
% % Simulation data
% clc;clear all;load CF.mat;load CFTF.mat;load CFTE.mat;load TF.mat;load TE.mat;
% data=[CF;TF;TE;CFTF;CFTE];
% input=data(:,5:end);
% biaoqian=data(:,[1:4]);
% input_train=input';
% output_train=biaoqian';
% % Input data normalization
% [inputn,minp,maxp,outputn,mint,maxt]=premnmx(input_train,output_train);%[inputn,inputps]=mapminmax(input_train);
% % BP network training 
% % % Initialize the network structure
% net=newff(inputn,output_train,13);%Number of nodes 13
% net.trainParam.epochs=1000;
% net.trainParam.lr=0.1;
% net.trainParam.goal=0.0000004;
% % network training
% net=train(net,inputn,output_train);
% % training set results
% BPoutput_train=sim(net,inputn);
% figure(1)
% plot(BPoutput_train(1,:),'r.')
% hold on
% plot(output_train(1,:),'b.')
% 
% save my_BP net;%Save the trained BP neural network net in my_BP net
% save my_mint mint;%Save the normalized minimum value of the trained BP neural network in my_mint mint 
% save my_maxt maxt;%Save the maximum normalized value of the trained BP neural network in my_maxt maxt
%% BP test
clc
clear
load my_BP.mat net;
load my_maxt;
load my_mint;
load Delta_sc.mat;
load Delta_sc_bq.mat;

% Test
input_test=Delta_sc';
inputn=premnmx(input_test);
an=sim(net,inputn);

figure(1)
plot(an(1,:),'r-')
hold on
x1_1=1:1789; y1_1=an(1,x1_1);
%cftool(x1_1,y1_1)
a1_1=0.9557;
b1_1=-1.122e-05;
nh1_1= a1_1*exp(b1_1*x1_1);
plot(x1_1,nh1_1,'b','linewidth',2)
set(gcf,'unit','centimeters','position',[7 7 7 5]);
hold on
x2_1=1790:3446; y2_1=an(1,x2_1);
%cftool(x2_1,y2_1)
a2_1=1.041;
b2_1= -2.853e-05;
nh2_1=a2_1*exp(b2_1*x2_1);
plot(x2_1,nh2_1,'b','linewidth',2)
hold on
x3_1=3447:4511; y3_1=an(1,x3_1);
%cftool(x3_1,y3_1)
a3_1=1.046;
b3_1=-1.967e-05;
nh3_1=a3_1*exp(b3_1*x3_1);
plot(x3_1,nh3_1,'b','linewidth',2)
hold on
x4_1=4512:length(Delta_sc); y4_1=an(1,x4_1);
%cftool(x4_1,y4_1)
a4_1=1.325;
b4_1=-6.31e-05;
nh4_1=a4_1*exp(b4_1*x4_1);
plot(x4_1,nh4_1,'b','linewidth',2)
xlabel('Operating time(hours)')
ylabel('σ_G_,_c')
ylim([ 0.9 1.1])


figure(2)
plot(an(2,:),'r-')
hold on
x1_2=1:1789; y1_2=an(2,x1_2);
%cftool(x1_2,y1_2)
a1_2=0.9778;
b1_2=-5.439e-06;
nh1_2= a1_2*exp(b1_2*x1_2);
plot(x1_2,nh1_2,'b','linewidth',2)
set(gcf,'unit','centimeters','position',[7 7 7 5]);
hold on
x2_2=1790:3446; y2_2=an(2,x2_2);
%cftool(x2_2,y2_2)
a2_2=1.02;
b2_2=-1.401e-05;
nh2_2=a2_2*exp(b2_2*x2_2);
plot(x2_2,nh2_2,'b','linewidth',2)
hold on
x3_2=3447:4511; y3_2=an(2,x3_2);
%cftool(x3_2,y3_2)
a3_2=1.022;
b3_2=-9.657e-06;
nh3_2=a3_2*exp(b3_2*x3_2);
plot(x3_2,nh3_2,'b','linewidth',2)
hold on
x4_2=4512:length(Delta_sc); y4_2=an(2,x4_2);
%cftool(x4_2,y4_2)
a4_2=1.148 ;
b4_2=-3.096e-05;
nh4_2=a4_2*exp(b4_2*x4_2);
plot(x4_2,nh4_2,'b','linewidth',2)
xlabel('Operating time(hours)')
ylabel('σ_η_,_c')
ylim([ 0.96 1.05])

figure(3)
plot(an(3,:),'r-')
hold on
x1_3=1:1789; y1_3=an(3,x1_3);
%cftool(x1_3,y1_3)
a1_3=1.041;
b1_3=1.219e-05;
nh1_3= a1_3*exp(b1_3*x1_3);
plot(x1_3,nh1_3,'b','linewidth',2)
set(gcf,'unit','centimeters','position',[7 7 7 5]);
hold on
x2_3=1790:3446; y2_3=an(3,x2_3);
%cftool(x2_3,y2_3)
a2_3= 0.9894;
b2_3=1.398e-05;
nh2_3=a2_3*exp(b2_3*x2_3);
plot(x2_3,nh2_3,'b','linewidth',2)
hold on
x3_3=3447:4511; y3_3=an(3,x3_3);
%cftool(x3_3,y3_3)
a3_3=0.9332;
b3_3=2.417e-05;
nh3_3=a3_3*exp(b3_3*x3_3);
plot(x3_3,nh3_3,'b','linewidth',2)
hold on
x4_3=4512:length(Delta_sc); y4_3=an(3,x4_3);
%cftool(x4_3,y4_3)
a4_3=0.6754;
b4_3=8.111e-05;
nh4_3=a4_3*exp(b4_3*x4_3);
plot(x4_3,nh4_3,'b','linewidth',2)
xlabel('Operating time(hours)')
ylabel('σ_G_,_T')
ylim([ 0.65 1.1])

 figure(4)
plot(an(4,:),'r-')
hold on
x1_4=1:1789; y1_4=an(4,x1_4);
%cftool(x1_4,y1_4)
a1_4=0.9783;
b1_4= -9.014e-06;
nh1_4= a1_4*exp(b1_4*x1_4);
plot(x1_4,nh1_4,'b','linewidth',2)
set(gcf,'unit','centimeters','position',[7 7 7 5]);
hold on
x2_4=1790:3446; y2_4=an(4,x2_4);
%cftool(x2_4,y2_4)
a2_4=1.025;
b2_4= -1.632e-05 ;
nh2_4=a2_4*exp(b2_4*x2_4);
plot(x2_4,nh2_4,'b','linewidth',2)
hold on
x3_4=3447:4511; y3_4=an(4,x3_4);
%cftool(x3_4,y3_4)
a3_4=1.042;
b3_4= -1.464e-05;
nh3_4=a3_4*exp(b3_4*x3_4);
plot(x3_4,nh3_4,'b','linewidth',2)
hold on
x4_4=4512:length(Delta_sc); y4_4=an(4,x4_4);
%cftool(x4_4,y4_4)
a4_4=1.211 ;
b4_4= -4.01e-05;
nh4_4=a4_4*exp(b4_4*x4_4);
plot(x4_4,nh4_4,'b','linewidth',2)
xlabel('Operating time(hours)')
ylabel('σ_η_,_T')   
ylim([ 0.95 1.15])