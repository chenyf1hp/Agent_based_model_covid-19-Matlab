clear all;
clc
close all;
% rng(1)
%% Default parameters 
%Hong Kong population_size : 7428300
pop_size = 150;
R0 = 2.5;           %basic reproduction number
mu = 1/2.9;         %a rate of moving from I to R
beta = R0 * mu;     %a rate of moving from S to I
contact_rate =1/pop_size;
days =30;

%% Initialization
% S = Susceptible       I = Infectious         R = Recovered 
% Region = Reg          Asy = Asymptomatic     Dur = Duration time
uID = (1:pop_size)';
S = zeros(pop_size,1);
I = zeros(pop_size,1);
R = zeros(pop_size,1);
Reg = zeros(pop_size,1);
Dur = zeros(pop_size,1);
people = [uID S I R Reg Dur];

%Defined by users:
S_ini_num = round(0.9*pop_size);
I_ini_num = round(0.1*pop_size);
R_ini_num = pop_size - S_ini_num - I_ini_num;

%%Generation Cases
rand_pop = randperm(pop_size);
S_ini_id = rand_pop(1:S_ini_num);
I_ini_id = rand_pop((S_ini_num + 1) : (S_ini_num + I_ini_num));
R_ini_id = rand_pop((S_ini_num + I_ini_num + 1): pop_size);
Inf_time_ini = poissrnd(1/mu,I_ini_num,1);

people(S_ini_id,2)=1;

people(I_ini_id,3)=1;
people(I_ini_id,6)=Inf_time_ini;

people(R_ini_id,4)=1;

S_n = count_s(people);
I_n = count_i(people);
R_n = count_r(people);

dS_s = zeros(days,1);
dI_s = zeros(days,1);
dR_s = zeros(days,1);
%% Run
for day = 1:days

    S_id = get_s_id(people);    
    I_id = get_i_id(people);
    I2R_id = check_i2r(people,I_id);
    I_keep_id = setdiff(I_id,I2R_id);
    
    dS_s(day) = - binornd(S_n(day), beta*I_n(day)/pop_size, 1,1);
    if abs(dS_s(day)) > S_n(day)
        dS_s(day) = -S_n(day);
    end
    dI_s(day) = -dS_s(day) - length(I2R_id);
    dR_s(day) = length(I2R_id);
    r = randperm(S_n(day));
    s_new_rand_id = S_id(r);
    i_new_id = s_new_rand_id(1:(-dS_s(day)));
    Inf_time = poissrnd(1/mu,(-dS_s(day)),1);

    people(i_new_id,2) = 0;
    people(i_new_id,3) = 1;
    people(i_new_id,6) = Inf_time;

    people(I2R_id,3) = 0;
    people(I2R_id,4) = 1;
    people(I_keep_id,6) = people(I_keep_id,6)-1;
    
    S_n = [S_n count_s(people)];
    I_n = [I_n count_i(people)];
    R_n = [R_n count_r(people)];

end

S_d = round(0.9*pop_size);
I_d = round(0.1*pop_size);
R_d = pop_size - S_ini_num - I_ini_num;
dS_d = [];
dI_d = [];
dR_d = [];

%% Deterministic model
for k = 1:days
    dS_d(k) = - round(beta * S_d(k) * I_d(k) / pop_size);
    dI_d(k) = round(beta * S_d(k) * I_d(k) / pop_size) -round( mu * I_d(k));
    dR_d(k) =round( mu * I_d(k));    
    S_d = [S_d S_d(k)+dS_d(k)];
    I_d = [I_d I_d(k)+dI_d(k)];
    R_d = [R_d R_d(k)+dR_d(k)];
end



%% Plot
% 
% day = 1:14;
% m = [S_n(day); I_n(day); R_n(day)];
% plot(day,m);
% x = -1;
% axis([x,x + 17,-1,pop_size+5]);
% grid on
% while 1
%     if x>max(day)
%         break;
%     end
%     x = x+0.1;
%     axis([x,x + 17,-1,pop_size+5]); %ÒÆ¶¯×ø±êÏµ
%     pause(0.1);
% end

day = 1:days;

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










%% Definition of function

% count Number
function s_number = count_s(x)
    s_number = sum(x(:,2)); 
end

function i_number = count_i(x)
    i_number = sum(x(:,3)); 
end

function r_number = count_r(x)
    r_number = sum(x(:,4)); 
end

function s_id = get_s_id(x)
    s_id =  find(x(:,2)); 
end

function i_id = get_i_id(x)
    i_id =  find(x(:,3)); 
end

function r_id = get_r_id(x)
    r_id =  find(x(:,4)); 
end



function i2r_id = check_i2r(x,i_id)
    dur_time = x(i_id,6);
    index_id = find(dur_time == 0);
    i2r_id = i_id(index_id);
end




