close all
day = 1:n_days;

%% SIR trend
figure(1)
plot(day,S_n(day),'red','LineWidth',2);
hold on;
plot(day,I_n(day),'k','LineWidth',2);
hold on;
plot(day,R_n(day),'b','LineWidth',2)
hold on;

plot(day,S_d(day),'red--','LineWidth',2);
hold on;
plot(day,I_d(day),'k--','LineWidth',2);
hold on;
plot(day,R_d(day),'b--','LineWidth',2)
hold on;

grid on
legend('S_s','I_s','R_s','S_d','I_d','R_d')

%% New_infected
figure(2)
plot(day,-dS_s(day),'k','LineWidth',2);
hold on;

plot(day,-dS_d(day),'k--','LineWidth',2);


grid on
legend('new_I_s','new_I_d')

%% Cumulative infected
figure(3)
plot(day,S_ini_num-S_n(day),'k','LineWidth',2);
hold on;

plot(day,S_ini_num-S_d(day),'k--','LineWidth',2);


grid on
legend('cum_I_s','cum_I_d')