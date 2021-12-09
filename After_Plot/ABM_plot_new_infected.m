day = 1:n_days;


plot(day,-dS_s(day),'k','LineWidth',2);
hold on;

plot(day,-dS_d(day),'k--','LineWidth',2);


grid on
legend('new_I_s','new_I_d')