
close all;
clear all;
clc

%% Parameters
rng('shuffle');
BW=500e6;                         % Bandwidth user link [Hz]
BW_Margin=10;                     % in [%]
Useful_BW=BW*(1-BW_Margin/100);   % Useful Bandwidth [Hz]
rolloff = 5/100;                 % Roll-off of the shaping filter, latest value is 5%
BW_with_RO = Useful_BW/(1+rolloff);
T=50;                             % Noise temperature at the antenna, 50 degree Kelvin
Tref=290;                         % Reference Temperature 290 degree kelvin
k_B = 1.3806503e-23;              % Boltzmann constant
NF=10^(2.278/10);                 % Noise Figure of the UT in linear scale
sigma = k_B*((NF-1)*Tref + T)*BW_with_RO;    % Noise Power per User


N = 7; %Number of transmit antennas
K = 7; %Number of users

Rmin_table = 0.2:0.2:2.4;

%Number of realizations in the Monte Carlo simulations
nbrOfMonteCarloRealizations = 200; %100;



PdB = 14.9202;
P = 10.^(PdB/10)*7; %Linear scale


%% Channel creating
% % Radius = 100;
% % sigma = 1e-11;
% % H = CreateChannel(K,N, Radius);
% % H = H./sqrt(sigma);
% % % load H_sat % 1 case tot

% Hcell = cell(nbrOfMonteCarloRealizations,1);
% load('True_Channel_Matrix_175Users_7Beams_200Realizations.mat')
% for m=1:nbrOfMonteCarloRealizations
%     H=[];
%     for k=1:K
%         H = [H; True_Channel_Matrix(ceil(rand(1)*25)+25*(k-1),:,1)];
%     end
%     H = H./sqrt(sigma);
%     Hcell{m}  = H;
% end
% save('Channel_sat_cell_40000','Hcell')
% dfdfdfd

load Channel_sat_cell_40000

rates0 = zeros(nbrOfMonteCarloRealizations,K);
sumrates0 = zeros(nbrOfMonteCarloRealizations,1);
powerAllocation0 = zeros(nbrOfMonteCarloRealizations,K);

rates1 = zeros(nbrOfMonteCarloRealizations,K);
sumrates1 = zeros(nbrOfMonteCarloRealizations,1);
powerAllocation1 = zeros(nbrOfMonteCarloRealizations,K);

rates2 = zeros(nbrOfMonteCarloRealizations,K);
sumrates2 = zeros(nbrOfMonteCarloRealizations,1);
powerAllocation2 = zeros(nbrOfMonteCarloRealizations,K);

rates3 = zeros(nbrOfMonteCarloRealizations,K);
sumrates3 = zeros(nbrOfMonteCarloRealizations,1);
powerAllocation3 = zeros(nbrOfMonteCarloRealizations,K);

rates_maxQ = zeros(nbrOfMonteCarloRealizations,K);
sumrates_maxQ = zeros(nbrOfMonteCarloRealizations,1);
powerAllocation_maxQ = zeros(nbrOfMonteCarloRealizations,K);

% kenh 122,200 dang co van de
%% Go through all channel realizations
for Rmin = Rmin_table
    for m =  1 :nbrOfMonteCarloRealizations
            H = Hcell{m};
            %% Step 1: Calculate waterfilling

            %Compute normalized beamforming vectors for ZFBF
            wZFBF = functionZFBF(H);

            csi = real(diag(H*wZFBF)).^2';
            
            % Equal power
            powerAllocationwEqual = P/K;
            signalGains = (powerAllocationwEqual.*csi)';
            rates = log2(1+signalGains)';
            rates0(m,:) = rates;
            sumrates0(m,:) = sum(rates);
            powerAllocation0(m,:) = powerAllocationwEqual;
            
            % Warter filling
            powerAllocationwZFBF = functionHeuristicPowerAllocation(csi,P,ones(length(csi),1));
            signalGains = (powerAllocationwZFBF.*csi)';
            rates = log2(1+signalGains)';
            rates1(m,:) = rates;
            sumrates1(m,:) = sum(rates);
            powerAllocation1(m,:) = powerAllocationwZFBF;
            
            
            %% Step 2: Processing:
            %% tinh PQ
            Q = abs(H*wZFBF).^2;
            alpha = 2^(Rmin) -1;
            v = zeros(K,1);
            for k=1:K
                v(k) = alpha/((alpha+1)*Q(k,k));
            end
            R = diag(v);
            R = R./sigma;
            Q =Q.*sigma;
            Pmin=(inv(eye(K)-R*Q)*v)';
            
            if(min(Pmin)<P)
                if max(abs(eig(R*Q)))<1 && sum(inv(eye(K)-R*Q)*v)<P %% toan bo user thoa man QoS
                    PtotalLeft = P - sum(Pmin);
                    powerAllocationwZFBF_left = functionHeuristicPowerAllocation(csi,PtotalLeft,ones(length(csi),1));
                    
                    powerAllocationwZFBF_maxQ = Pmin+PtotalLeft/K;
                    W_maxQ = kron(sqrt(powerAllocationwZFBF_maxQ),ones(N,1)).*wZFBF;
                    channelGains_maxQ = abs(H*W_maxQ).^2;
                    signalGains_maxQ = diag(channelGains_maxQ);
                    interferenceGains_maxQ = sum(channelGains_maxQ,2)-signalGains_maxQ;
                    rates_maxQ_temp = log2(1+signalGains_maxQ./(interferenceGains_maxQ+1))';

                    rates_maxQ(m,:) = rates_maxQ_temp;
                    sumrates_maxQ(m) = sum(rates_maxQ(m,:));
                    powerAllocation_maxQ(m,:) = powerAllocationwZFBF_maxQ;


                    %
                    powerAllocationwZFBF_final = Pmin + powerAllocationwZFBF_left;
                    W = kron(sqrt(powerAllocationwZFBF_final),ones(N,1)).*wZFBF;
                    channelGains = abs(H*W).^2;
                    signalGains = diag(channelGains);
                    interferenceGains = sum(channelGains,2)-signalGains;
                    rates_final = log2(1+signalGains./(interferenceGains+1))';
                    
                else
                    [Psort, index_Psort] = sort(Pmin);
                    csi_sort = csi(index_Psort);
                    for i = 1:length(csi)
                        if(sum(Psort(1:i)) > P)
                            selected_number = i-1;
                            break;
                        end
                    end
                    csi_left = csi_sort(selected_number+1:end);
                    PtotalLeft = P - sum(Psort(1:selected_number));
                    powerAllocationwZFBF_left = functionHeuristicPowerAllocation(csi_left,PtotalLeft,ones(length(csi_left),1));
                    signalGains_left = (powerAllocationwZFBF_left.*csi_left)';
                    rates_left = log2(1+signalGains_left)';
                    
                    rates_final = zeros(1,7); powerAllocationwZFBF_final= zeros(1,7);
                    rates_final(index_Psort(1:selected_number)) = Rmin;
                    rates_final(index_Psort(selected_number+1:end))=rates_left;
                    powerAllocationwZFBF_final(index_Psort(1:selected_number)) = Pmin(index_Psort(1:selected_number));
                    powerAllocationwZFBF_final(index_Psort(selected_number+1:end)) = powerAllocationwZFBF_left;
                    
                    rates_maxQ(m,:) = rates_final;
                    sumrates_maxQ(m) = sum(rates_final);
                    powerAllocation_maxQ(m,:) = powerAllocationwZFBF_final;
                end
                rates3(m,:) = rates_final;
                sumrates3(m,:) = sum(rates_final);
                powerAllocation3(m,:) = powerAllocationwZFBF_final;
                
            else
                rates_maxQ(m,:) = rates1(m,:);
                sumrates_maxQ(m) = sum(rates1(m,:));
                powerAllocation_maxQ(m,:) = powerAllocation1(m,:);
                
                rates3(m,:) = rates1(m,:);
                sumrates3(m,:) = sum(rates1(m,:));
                powerAllocation3(m,:) = powerAllocation1(m,:);
            end
        end
      
    
    
    
    
    filename=['Results/Rmin_ZF_trueSRM_',num2str(Rmin),'_200Samples_PQ.mat'];
%     save(filename,'rates0','rates1','rates2','rates3','sumrates0','sumrates1','sumrates2','sumrates3','powerAllocation0','powerAllocation1','powerAllocation2','powerAllocation3','rates_maxQ','sumrates_maxQ','powerAllocation_maxQ')
end

disp('Finish')

