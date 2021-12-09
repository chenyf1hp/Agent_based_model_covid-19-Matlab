day = 1:n_days;


plot(day,I_n(day),'k','LineWidth',2);
hold on;

plot(day,I_d(day),'k--','LineWidth',2);


grid on
legend('I_s','I_d')