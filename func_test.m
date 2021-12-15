minDelayC = 3.89;  maxDelayC = 8.93;
minDelayT = 1;     maxDelayT = 5.53;
kappaC = 9.5e-2;    thldC = 20.19;
kappaT = 4.2e-3;    thldT = 19.85;
p2 = 24.4/100;   %probability of back-traced
q = 11.4/100;    %quarantined to unquarantined
startDayOfContactTracing = 5;

if day >= startDayOfContactTracing
    
end




function timeOfContracing = f(kappaC,minDelayC,maxDelayC,nCaseC,thldC)
    timeOfContracing = minDelayC + (maxDelayC - minDelayC)*(1 - exp(-(nCaseC - thldC)/kappaC));
end

function timeOfTestDelay = g(kappaT,minDelayT,maxDelayT,nCaseT,thldT)
    timeOfTestDelay = minDelayT + (maxDelayT - minDelayT)*(1 - exp(-(nCaseT - thldT)/kappaT));
end




