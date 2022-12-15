%implementation that receives the cks (the channel noise ratio)
%or the channel state info with the power budget pt => compute the pk star
%through non-linear minimization and the pk star equation

function waterfilling(csi, Ptotal)
% csi: vector size 1xn

n = length(csi);
%compute mu star fuction
M = 1e4; %a user specified number of grid points 
mu_axis = linspace(1e-15, 5, M); 
Pn = max(1./mu_axis - (1./csi)', 0);
g = sum(log(1 + Pn.*(repmat(csi', [1, M])))) - mu_axis.*(sum(Pn) - Ptotal);
[min_g, idx] = find(g == min(g)); %min of g
mu = mu_axis(idx);

%optimal power allocation
Pn_opt = max(1./mu - 1./csi, 0);

f1 = figure(1);
clf(f1);
subplot(2, 1, 1) %power allocatied to each user
set(f1, 'Color', [1 1 1])
bar(Pn_opt, 1, 'red')
hold on
xlim = get(gca, 'xlim');
plot(xlim, [1/mu 1/mu], 'k--')
xlabel('Users')
ylabel('Power Allocated')

subplot(2, 1, 2) % plot the csi
set(f1, 'Color', [1 1 1])
bar(csi, 1)
xlabel('Users')
ylabel('CSI')
sgtitle(['Water filling Algorithm with total power = ', num2str(sum(Pn_opt))])
