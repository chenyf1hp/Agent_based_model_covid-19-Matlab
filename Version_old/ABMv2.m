clear all;
clc
close all;
rng(1)
%% Default parameters 
% Population parameters
pop_size =1500;               %Hong Kong population_size : 7428300
pop_infected = 2;   %round(0.1*pop_size)
% Simulation parameters
start_day = datetime('today');
end_day = 0;
n_days =30;

%NPI_impact factor

npi_factor = 1;


% Basic disease transmission parameters
R0 = 2.5;           %basic reproduction number
mu = 1/2.9;         %a rate of moving from I to R
beta = npi_factor* R0 * mu;     %a rate of moving from S to I
contact_rate =1/pop_size;
asymp_factor = 1;

% Efficacy of protection measures
iso_factor = 0;
quar_factor = 0;

%% Initialization Peoplemeta
% S = Susceptible       I = Infectious         R = Recovered 
% Region = Reg          Asy = Asymptomatic     Dur = Duration time
uid = (1:pop_size)';
age = zeros(pop_size,1);
sex = zeros(pop_size,1);
states_S = zeros(pop_size,1);             % S = Susceptible 
states_I = zeros(pop_size,1);             % I = Infectious
states_R = zeros(pop_size,1);             % R = Recovered
reg = zeros(pop_size,1);                  % Reg = region  
dur_inf2rec = zeros(pop_size,1);               % Duration for people with infectious to recover
inf_time = NaN(pop_size,1);             % inf_time = infected time
rec_time = NaN(pop_size,1);             % rec_time = recovery time
people = [uid age sex states_S states_I states_R reg dur_inf2rec inf_time rec_time];

%Defined by users:
S_ini_num = pop_size - pop_infected;
I_ini_num = pop_infected;
R_ini_num = pop_size - S_ini_num - I_ini_num;

%%Generation Cases
rand_pop = randperm(pop_size);
S_ini_id = rand_pop(1:S_ini_num);
I_ini_id = rand_pop((S_ini_num + 1) : (S_ini_num + I_ini_num));
R_ini_id = rand_pop((S_ini_num + I_ini_num + 1): pop_size);
Inf_time_ini = poissrnd(1/mu,I_ini_num,1);

people(S_ini_id,4)=1;

people(I_ini_id,5)=1;
people(I_ini_id,8)=Inf_time_ini;
people(I_ini_id,9)=Inf_time_ini;

people(R_ini_id,6)=1;

S_n = count_s(people);
I_n = count_i(people);
R_n = count_r(people);

dS_s = zeros(n_days,1);
dI_s = zeros(n_days,1);
dR_s = zeros(n_days,1);
%% Run
for day = 1:n_days

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

    people(i_new_id,4) = 0;
    people(i_new_id,5) = 1;
    people(i_new_id,8) = Inf_time;

    people(I2R_id,5) = 0;
    people(I2R_id,6) = 1;
    people(I_keep_id,8) = people(I_keep_id,8)-1;
    
    S_n = [S_n count_s(people)];
    I_n = [I_n count_i(people)];
    R_n = [R_n count_r(people)];

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
    dI_d(k) = round(beta * S_d(k) * I_d(k) / pop_size) -round( mu * I_d(k));
    dR_d(k) =round( mu * I_d(k));    
    S_d = [S_d S_d(k)+dS_d(k)];
    I_d = [I_d I_d(k)+dI_d(k)];
    R_d = [R_d R_d(k)+dR_d(k)];
end




%% Definition of function

% count Number
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



function i2r_id = check_i2r(x,i_id)
    dur_time = x(i_id,8);
    index_id = find(dur_time == 0);
    i2r_id = i_id(index_id);
end




