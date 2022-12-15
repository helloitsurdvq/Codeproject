function final_pwr_alloc_vec = waterfilling2(p, n_vec)
num_channels = numel(n_vec);
[sorted_noise, sort_idx] = sort(n_vec);
[unq_sorted_noise, last_unq_idx] = unique(sorted_noise);
num_unq = numel(unq_sorted_noise);
[num_reps, ~] = histogram(sorted_noise, unq_sorted_noise) ;
tot_pow_steps = [unq_sorted_noise(2:end).*cumsum(num_reps(1:num_unq-1)) - cumsum(unq_sorted_noise(1:end-1).*num_reps(1:end-1)) 0];
first_overflow_bin = find((p - tot_pow_steps)<0, 1);
if(~isempty(first_overflow_bin))
    init_nu = unq_sorted_noise(first_overflow_bin);
    init_pwr_alloc_vec = max(0, init_nu - sorted_noise);
    tot_pwr_allocated = sum(init_pwr_alloc_vec);
    res_pwr = p - tot_pwr_allocated;
    res_pwr_per_bin = res_pwr/last_unq_idx(first_overflow_bin);
    num_res_allocs = last_unq_idx(first_overflow_bin);
    final_sorted_pwr_alloc_vec = init_pwr_alloc_vec + res_pwr_per_bin*[ones(1, num_res_allocs) zeros(1, num_channels-num_res_allocs )];
    final_pwr_alloc_vec = zeros(1, num_channels);
    final_pwr_alloc_vec(sort_idx) = final_sorted_pwr_alloc_vec;
else
    final_pwr_alloc_vec = p - n_vec;
end
end