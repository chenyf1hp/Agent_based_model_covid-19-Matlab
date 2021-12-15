close all
day = 1:(n_days+1);

%% SIR trend
figure(1)
plot(day,S_n(day),'red','LineWidth',2);
hold on;
plot(day,Exp_n(day),'y','LineWidth',2);
hold on;
plot(day,I_n(day),'k','LineWidth',2);
hold on;
plot(day,R_n(day),'b','LineWidth',2);

grid on
legend('S','E','I','R')

%% New_infected
figure(2)
plot(day,dE(day),'k','LineWidth',2);

grid on
legend('New infected')

%% Cumulative infected
figure(3)
plot(day,Infected_cum(day),'k','LineWidth',2);

grid on
legend('Cumulative Infected')