day = 1:n_days;


plot(day,S_ini_num-S_n(day),'k','LineWidth',2);
hold on;

plot(day,S_ini_num-S_d(day),'k--','LineWidth',2);


grid on
legend('cum_I_s','cum_I_d')