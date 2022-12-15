% Number of channels
N = 20;
N0 = 1; % Normalized noise level
SNR_dB = 10; % Signal-to-noise ratio in dB
P = 10^(SNR_dB / 10);

g = abs(randn(N, 1));
%g = abs(2 15 3 17 5 5 5 6 4 5 11 11 8 9 10 14 6 2 3 13);
% Bisection search for alpha
alpha_low = min(N0 / g);
alpha_high = (P + sum(N0 / g)) / N; % Initial high (upper bound)

stop_threshold = 1e-5; % Stop threshold

% Iterate while low/high bounds are futher than stop_threshold
while(abs(alpha_low - alpha_high) > stop_threshold)
    alpha = (alpha_low + alpha_high) / 2; % Test value in the middle of low/high

    p = 1 / alpha - N0 / g;

    p(p < 0) = 0; % Consider only positive powers

    if (sum(p) > P)
        alpha_low = alpha;
    else
        alpha_high = alpha;
    end
end

% Print the achievable rate in bits/s
rates = log2(1 + g*p/N0);
disp(g)
disp(['Achievable rate ', num2str(sum(rates(:))), ' bits/s'])