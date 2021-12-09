num_dis=6;
cdf_mu =2;
    dis = 1:(num_dis+1);
    p = expcdf(dis,cdf_mu);
    diff_p = [p(1) diff(p)];
    diff_p_n = diff_p + diff_p * (1 - max(p));
    p_n = cumsum(diff_p_n);
    p_n(num_dis+1)=1;
