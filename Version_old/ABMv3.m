clear all;
clc
close all;
rng(1)

%% Global varibables

global people 


%% Default parameters 
% Population parameters
pop_size = 150;               %Hong Kong population_size : 7396000   https://www.censtatd.gov.hk/tc/web_table.html?id=216#
pop_infected = 2;           %round(0.1*pop_size)    Initial infected people

% Simulation parameters
start_day = datetime('today');
end_day = 0;
n_days = 30;

%NPI_impact factor

npi_factor = 1;


% Basic disease transmission parameters
R0 = 2.5;           %basic reproduction number
mu = 1/2.9;         %2.9 a rate of moving from I to R
beta = npi_factor * R0 * mu;     %a rate of moving from S to I
contact_rate =1/pop_size;
asymp_factor = 1;

% Efficacy of protection measures
iso_factor = 0;
quar_factor = 0;

%% Initialization Peoplemeta
% S = Susceptible       I = Infectious         R = Recovered 
% Region = Reg          Asy = Asymptomatic     Dur = Duration time
uid = randperm(pop_size)';                      % 1st
age = zeros(pop_size,1);                  % 2nd    https://www.censtatd.gov.hk/tc/web_table.html?id=216#
sex = zeros(pop_size,1);                  % 3rd
states_S = zeros(pop_size,1);             % S = Susceptible  4th
states_I = zeros(pop_size,1);             % I = Infectious   5th
states_R = zeros(pop_size,1);             % R = Recovered    6th
reg = zeros(pop_size,1);                  % Reg = region     7th
dur_inf2rec = zeros(pop_size,1);          % Duration for people with infectious to recover  8th
inf_time = NaN(pop_size,1);               % inf_time = infected time  9th
rec_time = NaN(pop_size,1);               % rec_time = recovery time  10th

people = [uid age sex states_S states_I states_R reg dur_inf2rec inf_time rec_time];

%Defined by users:
S_ini_num = pop_size - pop_infected;
I_ini_num = pop_infected;
R_ini_num = pop_size - S_ini_num - I_ini_num;

%%Generation Cases
S_ini_id = uid(1:S_ini_num);
I_ini_id = uid((S_ini_num + 1) : (S_ini_num + I_ini_num));
R_ini_id = uid((S_ini_num + I_ini_num + 1): pop_size);
Inf_time_ini = truncated_poisson(mu,1,I_ini_num);

people(S_ini_id,4) = 1;
% people(I_ini_id,8) = Inf_time_ini;
people = infect_happen(people,I_ini_id,Inf_time_ini,0);
people(R_ini_id,6) = 1;

S_n = count_s(people);
I_n = count_i(people);
R_n = count_r(people);

dS_s = zeros(n_days,1);
dI_s = zeros(n_days,1);
dR_s = zeros(n_days,1);


%% Generate pop distribution
pop_cdc_stat = csvread('numeric_data.csv');    
pop_proportion = round(pop_cdc_stat(:,5) * pop_size);
sum_generation = sum(pop_proportion);
if sum_generation ~= pop_size
    randsample = randi([1,length(pop_cdc_stat(:,5))],abs(pop_size-sum_generation),1);
    if sum_generation < pop_size
        pop_proportion(randsample) = pop_proportion(randsample) + 1;
    else
        pop_proportion(randsample) = pop_proportion(randsample) - 1;
    end
end

for pop_dis_id = 1:length(pop_cdc_stat(:,5))
    if pop_proportion(pop_dis_id) == 0
        continue;
    else
        up2now = sum(pop_proportion(1:pop_dis_id));
        start_index = up2now - pop_proportion(pop_dis_id) + 1;
        end_index = start_index + pop_proportion(pop_dis_id) - 1;
        people(start_index:end_index,2) = pop_cdc_stat(pop_dis_id,3);
        people(start_index:end_index,3) = pop_cdc_stat(pop_dis_id,2);
        people(start_index:end_index,7) = pop_cdc_stat(pop_dis_id,1);
    end
end


%% Run
for day = 1:n_days
%     tic
    S_id = get_s_id(people);    
    I_id = get_i_id(people);
    I2R_id = check_i2r(people,I_id);
    I_keep_id = setdiff(I_id,I2R_id);
    
    dS_s(day) = - binornd(S_n(day), beta * I_n(day) / pop_size, 1,1);
    if abs(dS_s(day)) > S_n(day)
        dS_s(day) = -S_n(day);
    end
    dI_s(day) = -dS_s(day) - length(I2R_id);
    dR_s(day) = length(I2R_id);
    r = randperm(S_n(day));
    s_new_rand_id = S_id(r);
    i_new_id = s_new_rand_id(1:(-dS_s(day)));
    Inf_time = truncated_poisson(mu,1,(-dS_s(day)));

    people(i_new_id,4) = 0;
    people = infect_happen(people,i_new_id,Inf_time,day);

    people(I2R_id,5) = 0;
    people(I2R_id,6) = 1;
    people(I_keep_id,8) = people(I_keep_id,8)-1;
    
    S_n = [S_n count_s(people)];
    I_n = [I_n count_i(people)];
    R_n = [R_n count_r(people)];
%     toc
end

S_d = pop_size - pop_infected;
I_d = pop_infected;
R_d = pop_size - S_d - I_d;
dS_d = [];
dI_d = [];
dR_d = [];

%% Deterministic model
for k = 1:n_days
    dS_d(k) = - round(beta * S_d(k) * I_d(k) / pop_size);
    dI_d(k) = round(beta * S_d(k) * I_d(k) / pop_size) - round( mu * I_d(k));
    dR_d(k) = round( mu * I_d(k));    
    S_d = [S_d S_d(k)+dS_d(k)];
    I_d = [I_d I_d(k)+dI_d(k)];
    R_d = [R_d R_d(k)+dR_d(k)];
end

A = get_state(people,7);


%% Definition of function

% Count Number
function s_number = count_s(x)
    s_number = sum(x(:,4)); 
end

function i_number = count_i(x)
    i_number = sum(x(:,5)); 
end

function r_number = count_r(x)
    r_number = sum(x(:,6)); 
end

function s_id = get_s_id(x)
    s_id =  find(x(:,4)); 
end

function i_id = get_i_id(x)
    i_id =  find(x(:,5)); 
end

function r_id = get_r_id(x)
    r_id =  find(x(:,6)); 
end

function whos_state = get_state(people,person_id)
    person_index = find(people(:,1) == person_id);
    whos_state = people(person_index,:); 
end


function i2r_id = check_i2r(x,i_id)
    dur_time = x(i_id,8);
    index_id = find(dur_time == 0);
    i2r_id = i_id(index_id);
end

function truncated_inf_time = truncated_poisson(mu,minvalue,size)
    truncated_inf_time = poissrnd(1/mu,size,1);
    truncated_inf_time(truncated_inf_time < minvalue) = minvalue;
end


function people = infect_happen(people,infected_id,i2r_time,n_day)
    people(infected_id,5) = 1;
    people(infected_id,8) = i2r_time;
    people(infected_id,9) = n_day;
    people(infected_id,10)= n_day + i2r_time;
end



