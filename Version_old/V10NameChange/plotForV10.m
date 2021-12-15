close all
day = 1:(nDays+1);

%% SIR trend
figure(1)
plot(day,nS(day),'red','LineWidth',2);
hold on;
plot(day,nE(day),'y','LineWidth',2);
hold on;
plot(day,nInfectious(day),'k','LineWidth',2);
hold on;
plot(day,nRe(day),'b','LineWidth',2);

grid on
legend('S','E','I','R')

%% New_infected
figure(2)
plot(day,nEPerDay(day),'k','LineWidth',2);

grid on
legend('New infected')

%% Cumulative infected
figure(3)
plot(day,nOfCumInfected(day),'k','LineWidth',2);

grid on
legend('Cumulative Infected')