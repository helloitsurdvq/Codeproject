    

close all
clear
clc
K=7;
Rmin_table = 0.2:0.2:2.4;  %[0.5 1 1.5 2 2.5 3];
ProbTable_ZF = zeros(length(Rmin_table), 4);
SumrateTable_ZF = zeros(length(Rmin_table),4);
Jain_WF_ZF = zeros(length(Rmin_table),200);
Jain_alg_ZF = zeros(length(Rmin_table),200);
Jain_equal_ZF = zeros(length(Rmin_table),200);
Jain_maxQ_ZF = zeros(length(Rmin_table),200);

ProbTable_RZF = zeros(length(Rmin_table), 4);
SumrateTable_RZF = zeros(length(Rmin_table),4);
Jain_WF_RZF = zeros(length(Rmin_table),200);
Jain_alg_RZF = zeros(length(Rmin_table),200);
Jain_equal_RZF = zeros(length(Rmin_table),200);
Jain_maxQ_RZF = zeros(length(Rmin_table),200);

ProbTable_MRT = zeros(length(Rmin_table), 4);
SumrateTable_MRT = zeros(length(Rmin_table),4);
Jain_WF_MRT = zeros(length(Rmin_table),200);
Jain_alg_MRT = zeros(length(Rmin_table),200);
Jain_equal_MRT = zeros(length(Rmin_table),200);
Jain_maxQ_MRT = zeros(length(Rmin_table),200);

for i = 1:length(Rmin_table)
    %% MRT
    Rmin=Rmin_table(i);
%     filename_ZF = ['Rmin_RZF2_',num2str(Rmin),'_200Samples.mat'];
    filename_MRT = ['Rmin_MRT2_',num2str(Rmin),'_200Samples.mat'];
    Rmin = Rmin-1e-5;
%     filename = ['Rmin_',num2str(Rmin),'.mat'];
    load(filename_MRT);
    ProbTable_MRT(i,1) = length(find(rates0>=Rmin))/(1400);
    ProbTable_MRT(i,2) = length(find(rates1>=Rmin))/(1400);
    ProbTable_MRT(i,3) = length(find(rates3>=Rmin))/(1400);
    ProbTable_MRT(i,4) = length(find(rates_maxQ>=Rmin))/(1400);
    
    SumrateTable_MRT(i,1) = mean(sumrates0);
    SumrateTable_MRT(i,2) = mean(sumrates1);
    SumrateTable_MRT(i,3) = mean(sumrates3);
    SumrateTable_MRT(i,4) = mean(sumrates_maxQ);
   
    u_WF = rates1./Rmin;
    Jain_WF_MRT(i,:) = (sum(rates1,2).^2./(sum(rates1.^2,2)*K))';
    
    u_alg = rates3./Rmin;
    Jain_alg_MRT(i,:) = (sum(rates3,2).^2./(sum(rates3.^2,2)*K))';
    
    u_equal = rates0./Rmin;
    Jain_equal_MRT(i,:) = (sum(u_equal,2).^2./(sum(u_equal.^2,2)*K))';
    
    u_maxQ = rates_maxQ./Rmin;
    Jain_maxQ_MRT(i,:) = (sum(u_maxQ,2).^2./(sum(u_maxQ.^2,2)*K))';
    
    %% ZF
    Rmin=Rmin_table(i);
    filename_ZF = ['Rmin_ZF_trueSRM_',num2str(Rmin),'_200Samples_PQ.mat'];
    Rmin = Rmin-1e-5;
%     filename = ['Rmin_',num2str(Rmin),'.mat'];
    load(filename_ZF);
    ProbTable_ZF(i,1) = length(find(rates0>=Rmin))/(1400);
    ProbTable_ZF(i,2) = length(find(rates1>=Rmin))/(1400);
    ProbTable_ZF(i,3) = length(find(rates3>=Rmin))/(1400);
    ProbTable_ZF(i,4) = length(find(rates_maxQ>=Rmin))/(1400);
    
    SumrateTable_ZF(i,1) = mean(sumrates0);
    SumrateTable_ZF(i,2) = mean(sumrates1);
    SumrateTable_ZF(i,3) = mean(sumrates3);
    SumrateTable_ZF(i,4) = mean(sumrates_maxQ);
   
    u_WF = rates1./Rmin;
    Jain_WF_ZF(i,:) = (sum(rates1,2).^2./(sum(rates1.^2,2)*K))';
    
    u_alg = rates3./Rmin;
    Jain_alg_ZF(i,:) = (sum(rates3,2).^2./(sum(rates3.^2,2)*K))';
    
    u_equal = rates0./Rmin;
    Jain_equal_ZF(i,:) = (sum(u_equal,2).^2./(sum(u_equal.^2,2)*K))';
    
    u_maxQ = rates_maxQ./Rmin;
    Jain_maxQ_ZF(i,:) = (sum(u_maxQ,2).^2./(sum(u_maxQ.^2,2)*K))';
    
    %% ZRF    
    Rmin=Rmin_table(i);
    filename_RZF = ['Rmin_RZF2_',num2str(Rmin),'_200Samples.mat'];
    Rmin = Rmin-1e-5;
    load(filename_RZF);
    ProbTable_RZF(i,1) = length(find(rates0>=Rmin))/(1400);
    ProbTable_RZF(i,2) = length(find(rates1>=Rmin))/(1400);
    ProbTable_RZF(i,3) = length(find(rates3>=Rmin))/(1400);
    ProbTable_RZF(i,4) = length(find(rates_maxQ>=Rmin))/(1400);
    
    SumrateTable_RZF(i,1) = mean(sumrates0);
    SumrateTable_RZF(i,2) = mean(sumrates1);
    SumrateTable_RZF(i,3) = mean(sumrates3);
    SumrateTable_RZF(i,4) = mean(sumrates_maxQ);
   
    u_WF = rates1./Rmin;
    Jain_WF_RZF(i,:) = (sum(u_WF,2).^2./(sum(u_WF.^2,2)*K))';
    
    u_alg = rates3./Rmin;
    Jain_alg_RZF(i,:) = (sum(u_alg,2).^2./(sum(u_alg.^2,2)*K))';
    
    u_equal = rates0./Rmin;
    Jain_equal_RZF(i,:) = (sum(u_equal,2).^2./(sum(u_equal.^2,2)*K))';
    
    u_maxQ = rates_maxQ./Rmin;
    Jain_maxQ_RZF(i,:) = (sum(u_maxQ,2).^2./(sum(u_maxQ.^2,2)*K))';
end

Rmin_table = Rmin_table*500;
%% Plot
figure(1)
hold on

plot(Rmin_table, ProbTable_ZF(:,1),'b--', 'linewidth', 1, 'DisplayName', 'ZF-EqualPower')
plot(Rmin_table, ProbTable_ZF(:,2),'b-', 'linewidth', 1, 'DisplayName', 'ZF-SumOpt')
plot(Rmin_table, ProbTable_ZF(:,4),'b-o', 'linewidth', 1, 'DisplayName', 'ZF-SatisSetOpt')
plot(Rmin_table, ProbTable_ZF(:,3),'b-*', 'linewidth', 1, 'DisplayName', 'ZF-JointOpt')

plot(Rmin_table, ProbTable_RZF(:,1),'r--', 'linewidth', 1, 'DisplayName', 'RZF-EqualPower')
plot(Rmin_table, ProbTable_RZF(:,2),'r-', 'linewidth', 1, 'DisplayName', 'RZF-SumOpt')
plot(Rmin_table, ProbTable_RZF(:,4),'r-o', 'linewidth', 1, 'DisplayName', 'RZF-SatisSetOpt')
plot(Rmin_table, ProbTable_RZF(:,3),'r-*', 'linewidth', 1, 'DisplayName', 'RZF-JointOpt')
grid on
% grid minor
xlim([Rmin_table(1) Rmin_table(end)])
xlabel('QoS requirement [Mbps]')
ylabel('Probility of satisfied demand-based constraints')
lgd = legend;
lgd.NumColumns = 2;


%%
figure(2)
hold on

plot(Rmin_table, SumrateTable_ZF(:,1).*500, 'b--', 'linewidth', 1, 'DisplayName', 'ZF-EqualPower')
plot(Rmin_table, SumrateTable_ZF(:,2).*500, 'b-', 'linewidth', 1, 'DisplayName', 'ZF-SumOpt')
plot(Rmin_table, SumrateTable_ZF(:,4).*500, 'b-o', 'linewidth', 1, 'DisplayName', 'ZF-SatisSetOpt')
plot(Rmin_table, SumrateTable_ZF(:,3).*500, 'b-*', 'linewidth', 1, 'DisplayName', 'ZF-JointOpt')

plot(Rmin_table, SumrateTable_RZF(:,1).*500, 'r--', 'linewidth', 1, 'DisplayName', 'RZF-EqualPower')
plot(Rmin_table, SumrateTable_RZF(:,2).*500, 'r-', 'linewidth', 1, 'DisplayName', 'RZF-SumOpt')
plot(Rmin_table, SumrateTable_RZF(:,4).*500, 'r-o', 'linewidth', 1, 'DisplayName', 'RZF-SatisSetOpt')
plot(Rmin_table, SumrateTable_RZF(:,3).*500, 'r-*', 'linewidth', 1, 'DisplayName', 'RZF-JointOpt')

grid on
% grid minor
xlabel('QoS requirement [Mbps]')
ylabel('Sum rate [Mbps]')
xlim([Rmin_table(1) Rmin_table(end)])

%%
figure(3); hold on;

f1=cdfplot(Jain_alg_ZF(1,:)); set(f1,'color','b');
f3=cdfplot(Jain_alg_ZF(5,:)); set(f3,'color','r');
f4=cdfplot(Jain_alg_ZF(9,:)); set(f4,'color','g');
f2=cdfplot(Jain_WF_ZF(9,:)); set(f2,'color','k');
legend('ZF - \xi = 100','ZF - \xi = 300','ZF - \xi = 500','WF');

%%
figure(4);
hold on
a= Rmin_table; % 0.5:0.5:3;
plot(a, mean(Jain_equal_ZF,2),'b--','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-EqualPower')
plot(a, mean(Jain_WF_ZF,2),'b','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-SumOpt')
plot(a, mean(Jain_maxQ_ZF,2),'b-o','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-SatisSetOpt')
plot(a, mean(Jain_alg_ZF,2),'b-*','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-JointOpt')

plot(a, mean(Jain_equal_RZF,2),'r--','linewidth',1,'MarkerSize',4, 'DisplayName', 'RZF-EqualPower')
plot(a, mean(Jain_WF_RZF,2),'r-','linewidth',1,'MarkerSize',4, 'DisplayName', 'RZF-SumOpt')
plot(a, mean(Jain_maxQ_RZF,2),'r-o','linewidth',1,'MarkerSize',4, 'DisplayName', 'RZF-SatisSetOpt')
plot(a, mean(Jain_alg_RZF,2),'r-*','linewidth',1,'MarkerSize',4,'DisplayName', 'RZF-JointOpt')

xlim([Rmin_table(1) Rmin_table(end)])

xlabel('QoS requirement [Mbps]')
ylabel('Jain index')
lgd = legend;
lgd.NumColumns = 2;
grid on

%%
figure
hold on
plot(Rmin_table, (ProbTable_ZF(:,1) + SumrateTable_ZF(:,1)./SumrateTable_ZF(:,2))/2,'b--','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-EqualPower')
plot(Rmin_table, (ProbTable_ZF(:,2) + SumrateTable_ZF(:,2)./SumrateTable_ZF(:,2))/2,'b','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-SumOpt')
plot(Rmin_table, (ProbTable_ZF(:,4) + SumrateTable_ZF(:,4)./SumrateTable_ZF(:,2))/2,'b-o','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-SatisSetOpt')
plot(Rmin_table, (ProbTable_ZF(:,3) + SumrateTable_ZF(:,3)./SumrateTable_ZF(:,2))/2,'b-*','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-JointOpt')

plot(Rmin_table, (ProbTable_RZF(:,1) + SumrateTable_RZF(:,1)./SumrateTable_RZF(:,2))/2,'r--','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-EqualPower')
plot(Rmin_table, (ProbTable_RZF(:,2) + SumrateTable_RZF(:,2)./SumrateTable_RZF(:,2))/2,'r','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-SumOpt')
plot(Rmin_table, (ProbTable_RZF(:,4) + SumrateTable_RZF(:,4)./SumrateTable_RZF(:,2))/2,'r-o','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-SatisSetOpt')
plot(Rmin_table, (ProbTable_RZF(:,3) + SumrateTable_RZF(:,3)./SumrateTable_RZF(:,2))/2,'r-*','linewidth',1,'MarkerSize',4, 'DisplayName', 'ZF-JointOpt')

xlim([Rmin_table(1) Rmin_table(end)])
xlabel('QoS requirement [Mbps]')
ylabel('Normalized objective function')
lgd = legend;
lgd.NumColumns = 2;
grid on
