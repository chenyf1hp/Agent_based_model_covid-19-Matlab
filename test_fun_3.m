exp_indx = [];
for exp_count = 1:3
    exp_index = [exp_indx find(people(:,1)==exp_uid(exp_count))];
end
% disp(exp_index)
% latent_time = truncated_poisson(epsilon,1,exp_n);
% pre_time = truncated_poisson(gamma,1,exp_n);
% inf_time = truncated_poisson(mu,1,exp_n);
% 
% people(exp_index,4) = 0;
% people(exp_index,5) = 1;
% people(exp_index,11) = latent_time;
% people(exp_index,12) = pre_time;
% people(exp_index,13) = inf_time;
% people(exp_index,14) = n_day;
% people(exp_index,15) = n_day + latent_time;
% people(exp_index,16) = n_day + latent_time + pre_time;
% people(exp_index,17) = n_day + latent_time + pre_time + inf_time;