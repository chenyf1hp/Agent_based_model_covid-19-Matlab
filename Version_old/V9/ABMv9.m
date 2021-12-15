clear all;
clc
close all;
rng(1)


%% load population
load("Population_Generation.mat")


%% Default parameters 
% Population parameters   set by users or import
exp_initial = 2;
pre_initial = 1;
asym_initial = 2;              
sym_initial = 2;
initial_infected = exp_initial + pre_initial + asym_initial + sym_initial;

% Simulation parameters
n_days = 15;

% Basic disease transmission parameters
R0 = 2.5;           %basic reproduction number
mu = 1/2.9;         %2.9 a rate of moving from I to R
r = 0.55;          %ratio of transmission rate for asymptomatic over symptomatic cases
epsilon = 1/2.9;   %a rate of moving from E to P
gamma = 1/2.3;     %a rate of moving from P to I
p = 30.8/100;      %Asymptomatic rate
beta = R0*mu/(mu*r/gamma+p*r+1-p);   %a rate of moving from S to E  0.6644


N_mean_contact = 10;
Var_contact = 0.6;


% Defined by users:
S_ini_num = pop_size - exp_initial - pre_initial - asym_initial - sym_initial;
I_ini_num = exp_initial + pre_initial + asym_initial + sym_initial;    %infected
R_ini_num = pop_size - S_ini_num - I_ini_num;

%% Generation Cases
S_ini_id = uid(1:S_ini_num);
Exp_ini_id = uid((S_ini_num + 1) : (S_ini_num + exp_initial));
Pre_ini_id = uid((S_ini_num + exp_initial + 1) : (S_ini_num + exp_initial + pre_initial));
Asym_ini_id = uid((S_ini_num + exp_initial + pre_initial + 1) : (S_ini_num + exp_initial + pre_initial + asym_initial));
Sym_ini_id = uid((S_ini_num + exp_initial + pre_initial + asym_initial + 1) : (S_ini_num + I_ini_num));
R_ini_id = uid((S_ini_num + I_ini_num + 1): pop_size);
latent_time_ini = truncated_poisson(epsilon,1,exp_initial);
pre_time_ini = truncated_poisson(gamma,1,exp_initial+pre_initial);
inf_time_ini = truncated_poisson(mu,1,I_ini_num);

people(S_ini_id,4) = 1;
people = infected_happen(people,Exp_ini_id,Pre_ini_id,Asym_ini_id,Sym_ini_id,latent_time_ini,pre_time_ini,inf_time_ini,0);
people(R_ini_id,9) = 1;

S_n = count_s(people);
Exp_n = count_exp(people);
Pre_n = count_pre(people);
Asym_n = count_asym(people);
Sym_n = count_sym(people);
R_n = count_r(people);
dE = 0;

s_color = [55,126,184]/256;
exp_color = [152,78,163]/256;
pre_color = [228,26,28]/256;
i_color = [228,26,28]/256;
r_color = [77,175,74]/256;


transmission_source = [];
transmission_chain =[];



%% Import cases
% import_cases = csvread('Input_cases.csv');    % 1st col distriction  2nd col infected_start_day
% import_cases_inf_start = import_cases(:,2);


%% Run
% figure
% mapshow(S, 'FaceColor', 'white');
% % mapshow(coordinate_points_set(S_ini_id), 'FaceColor', s_color);
% mapshow(coordinate_points_set(Exp_ini_id), 'FaceColor', exp_color);
% mapshow(coordinate_points_set(Pre_ini_id), 'FaceColor', pre_color);
% mapshow(coordinate_points_set(Asym_ini_id), 'FaceColor', i_color);
% mapshow(coordinate_points_set(Sym_ini_id), 'FaceColor', i_color);
% mapshow(coordinate_points_set(R_ini_id), 'FaceColor', r_color);
chain_count = 1;

for day = 1:n_days
    tic
    dE_day = 0;
    exp_new_index =[];
    pre_new_index =[];
    asym_new_index =[];
    sym_new_index =[];

    s_index = get_s_id(people);
    exp_index = get_exp_id(people);
    pre_index = get_pre_id(people);
    asym_index = get_asym_id(people);    %index
    sym_index = get_sym_id(people);
    r_index = get_r_id(people);
    
    e2p_index = check_e2p(people,exp_index);
    p2i_index = check_p2i(people,pre_index);
    asym2r_index = check_asym2r(people,asym_index);
    sym2r_index = check_sym2r(people,sym_index);
    exp_keep_index = setdiff(exp_index,e2p_index);
    pre_keep_index = setdiff(pre_index,p2i_index);
    asym_keep_index = setdiff(asym_index,asym2r_index);
    sym_keep_index = setdiff(sym_index,sym2r_index);
    
    N_contact = round(lognrnd(log(N_mean_contact),Var_contact,Pre_n(day)+Asym_n(day)+Sym_n(day),1));
    infectious_count = 1;
    for sym_count = 1:length(sym_index)
        infectious_flag = 3;
        origin_dis = people(sym_index(sym_count),10);
        sym_contact_n = N_contact(infectious_count);
        source_uid = people(sym_index(sym_count),1);
        dist_u = rand(sym_contact_n,1);
        [contact_index, contact_uid, infected_index, infected_uid]  = get_contact_uid(people,cum_p_contact,dis_index_pool,dist_u,sym_contact_n,origin_dis,infectious_flag,r,beta);
        sym_source = [infected_uid,repmat(source_uid,length(infected_uid),1)];
        transmission_source = [transmission_source;sym_source];
        if isempty(transmission_chain)
            transmission_chain(chain_count).source_uid = source_uid;
            transmission_chain(chain_count).contact_uid = contact_uid;
            transmission_chain(chain_count).infected_uid = infected_uid;
            chain_count = chain_count + 1;
        else
            exist_source_uid = cat(1,transmission_chain.source_uid);
            if find(exist_source_uid==source_uid)
                source_index = find(exist_source_uid==source_uid);
                chain_contact_uid = [transmission_chain(source_index).contact_uid;contact_uid];
                chain_contact_uid = unique(chain_contact_uid);
                transmission_chain(source_index).contact_uid = chain_contact_uid;
                chain_infected_uid = [transmission_chain(source_index).infected_uid;infected_uid];
                transmission_chain(source_index).infected_uid = chain_infected_uid;
            else
                transmission_chain(chain_count).source_uid = source_uid;
                transmission_chain(chain_count).contact_uid = contact_uid;
                transmission_chain(chain_count).infected_uid = infected_uid;   
                chain_count = chain_count + 1;
            end
        end
        dE_day = dE_day + length(infected_uid);
        people = exp_happen(people,infected_uid,length(infected_uid),day,epsilon,gamma,mu);
        infectious_count = infectious_count + 1;
        
    end
    if length(unique(transmission_source(:,1)))~=size(transmission_source,1)
        disp(day)
        disp('symend')
        break
    end    
    for asym_count = 1:length(asym_index)
        infectious_flag = 2;
        origin_dis = people(asym_index(asym_count),10);
        asym_contact_n = N_contact(infectious_count);
        source_uid = people(asym_index(asym_count),1);
        dist_u = rand(asym_contact_n,1);
        [contact_index, contact_uid, infected_index, infected_uid]  = get_contact_uid(people,cum_p_contact,dis_index_pool,dist_u,asym_contact_n,origin_dis,infectious_flag,r,beta);
        asym_source = [infected_uid,repmat(source_uid,length(infected_uid),1)];
        transmission_source = [transmission_source;asym_source];
        if isempty(transmission_chain)
            transmission_chain(chain_count).source_uid = source_uid;
            transmission_chain(chain_count).contact_uid = contact_uid;
            transmission_chain(chain_count).infected_uid = infected_uid;
            chain_count = chain_count + 1;
        else
            exist_source_uid = cat(1,transmission_chain.source_uid);
            if find(exist_source_uid==source_uid)
                source_index = find(exist_source_uid==source_uid);
                chain_contact_uid = [transmission_chain(source_index).contact_uid;contact_uid];
                chain_contact_uid = unique(chain_contact_uid);
                transmission_chain(source_index).contact_uid = chain_contact_uid;
                chain_infected_uid = [transmission_chain(source_index).infected_uid;infected_uid];
                transmission_chain(source_index).infected_uid = chain_infected_uid;
            else
                transmission_chain(chain_count).source_uid = source_uid;
                transmission_chain(chain_count).contact_uid = contact_uid;
                transmission_chain(chain_count).infected_uid = infected_uid;
                chain_count = chain_count + 1;
            end
        end
        dE_day = dE_day + length(infected_uid);
        people = exp_happen(people,infected_uid,length(infected_uid),day,epsilon,gamma,mu);
        infectious_count = infectious_count + 1;
    end
    if length(unique(transmission_source(:,1)))~=size(transmission_source,1)
        disp(day)
        disp('asymend')
        break
    end 
    for pre_count = 1:length(pre_index)
        infectious_flag = 1;
        origin_dis = people(pre_index(pre_count),10);
        pre_contact_n = N_contact(infectious_count);
        source_uid = people(pre_index(pre_count),1);
        dist_u = rand(pre_contact_n,1);
        [contact_index, contact_uid, infected_index, infected_uid]  = get_contact_uid(people,cum_p_contact,dis_index_pool,dist_u,pre_contact_n,origin_dis,infectious_flag,r,beta);
        pre_source = [infected_uid,repmat(source_uid,length(infected_uid),1)];
        transmission_source = [transmission_source;pre_source];
        if isempty(transmission_chain)
            transmission_chain(chain_count).source_uid = source_uid;
            transmission_chain(chain_count).contact_uid = contact_uid;
            transmission_chain(chain_count).infected_uid = infected_uid;
            chain_count = chain_count + 1;
        else
            exist_source_uid = cat(1,transmission_chain.source_uid);
            if find(exist_source_uid==source_uid)
                source_index = find(exist_source_uid==source_uid);
                chain_contact_uid = [transmission_chain(source_index).contact_uid;contact_uid];
                chain_contact_uid = unique(chain_contact_uid);
                transmission_chain(source_index).contact_uid = chain_contact_uid;
                chain_infected_uid = [transmission_chain(source_index).infected_uid;infected_uid];
                transmission_chain(source_index).infected_uid = chain_infected_uid;
            else
                transmission_chain(chain_count).source_uid = source_uid;
                transmission_chain(chain_count).contact_uid = contact_uid;
                transmission_chain(chain_count).infected_uid = infected_uid;
                chain_count = chain_count + 1;
            end
        end
        dE_day = dE_day + length(infected_uid);
        people = exp_happen(people,infected_uid,length(infected_uid),day,epsilon,gamma,mu);
        infectious_count = infectious_count + 1;
    end
    if length(unique(transmission_source(:,1)))~=size(transmission_source,1)
        disp(day)
        disp('preend')
        break
    end 

    pre2i = rand(numel(p2i_index),1);
    pre2i(pre2i>p) = 1;
    pre2i(pre2i<=p) = 0;

    people(e2p_index,5) = 0;
    people(e2p_index,6) = 1;
    people(p2i_index,6) = 0;
    people(p2i_index,7) = 1-pre2i;
    people(p2i_index,8) = pre2i;
    people(asym2r_index,7) = 0;
    people(asym2r_index,9) = 1;
    people(sym2r_index,8) = 0;
    people(sym2r_index,9) = 1;


    people(exp_keep_index,11) = people(exp_keep_index,11)-1;
    people(pre_keep_index,12) = people(pre_keep_index,12)-1;
    people(asym_keep_index,13)= people(asym_keep_index,13)-1;
    people(sym_keep_index,13) = people(sym_keep_index,13)-1;



    S_n = [S_n count_s(people)];
    Exp_n = [Exp_n count_exp(people)];
    Pre_n = [Pre_n count_pre(people)];
    Asym_n =[Asym_n count_asym(people)];
    Sym_n =[Sym_n count_sym(people)];
    R_n = [R_n count_r(people)];
    
    dE = [dE; dE_day];
    

%     s_now_index = get_s_id(people);
%     s2e_index = setdiff(s_index, s_now_index);
%     s2e_coordinate = coordinate_points_set(s2e_index);
%     mapshow(s2e_coordinate,'FaceColor', exp_color); 
% 
%     e2p_coordinate = coordinate_points_set(e2p_index);
%     mapshow(e2p_coordinate,'FaceColor', pre_color); 
% 
%     p2i_coordinate = coordinate_points_set(p2i_index);
%     mapshow(p2i_coordinate,'FaceColor', i_color); 
% 
%     asym2r_coordinate = coordinate_points_set(asym2r_index);
%     mapshow(asym2r_coordinate,'FaceColor', r_color); 
% 
%     sym2r_coordinate = coordinate_points_set(sym2r_index);
%     mapshow(sym2r_coordinate,'FaceColor', r_color); 
%     pause(0.5)
    toc
        

end

Infected_cum = cumsum(dE) + initial_infected;
I_n = Pre_n + Asym_n + Sym_n;


%% Definition of function

% Count Number
function N_number_dis = count_N_dis(x,dis,flag_work)
    if flag_work ==1 
        dis_index = find(x(:,22)==dis);
        N_number_dis = length(dis_index);
    else
        dis_index = find(x(:,10)==dis);
        N_number_dis = length(dis_index);
    end
end

function S_number_dis = count_S_dis(x,dis,flag_work)
    if flag_work ==1 
        dis_index = find(x(:,22)==dis);
        S_number_dis = sum(x(dis_index,4));
    else
        dis_index = find(x(:,10)==dis);
        S_number_dis = sum(x(dis_index,4));
    end
end

function Exp_number_dis = count_Exp_dis(x,dis,flag_work)
    if flag_work ==1 
        dis_index = find(x(:,22)==dis);
        Exp_number_dis = sum(x(dis_index,5));
    else
        dis_index = find(x(:,10)==dis);
        Exp_number_dis = sum(x(dis_index,5));
    end
end

function Pre_number_dis = count_Pre_dis(x,dis,flag_work)
    if flag_work ==1 
        dis_index = find(x(:,22)==dis);
        Pre_number_dis = sum(x(dis_index,6));
    else
        dis_index = find(x(:,10)==dis);
        Pre_number_dis = sum(x(dis_index,6));
    end
end

function Asym_number_dis = count_Asym_dis(x,dis,flag_work)
    if flag_work ==1 
        dis_index = find(x(:,22)==dis);
        Asym_number_dis = sum(x(dis_index,7));
    else
        dis_index = find(x(:,10)==dis);
        Asym_number_dis = sum(x(dis_index,7));
    end
end

function Sym_number_dis = count_Sym_dis(x,dis,flag_work)
    if flag_work ==1 
        dis_index = find(x(:,22)==dis);
        Sym_number_dis = sum(x(dis_index,8));
    else
        dis_index = find(x(:,10)==dis);
        Sym_number_dis = sum(x(dis_index,8));
    end
end

function R_number_dis = count_R_dis(x,dis,flag_work)
    if flag_work ==1 
        dis_index = find(x(:,22)==dis);
        R_number_dis = sum(x(dis_index,9));
    else
        dis_index = find(x(:,10)==dis);
        R_number_dis = sum(x(dis_index,9));
    end
end

function N_number_dis_sunday = count_N_dis_sunday(x,dis)
    dis_index = find(x(:,23)==dis);
	N_number_dis_sunday = length(dis_index);
end

function S_number_dis_sunday = count_S_dis_sunday(x,dis)
    dis_index = find(x(:,23)==dis);
	S_number_dis_sunday = sum(x(dis_index,4));
end

function Exp_number_dis_sunday = count_Exp_dis_sunday(x,dis)
    dis_index = find(x(:,23)==dis);
	Exp_number_dis_sunday = sum(x(dis_index,5));
end

function Pre_number_dis_sunday = count_Pre_dis_sunday(x,dis)
    dis_index = find(x(:,23)==dis);
	Pre_number_dis_sunday = sum(x(dis_index,6));
end

function Asym_number_dis_sunday = count_Asym_dis_sunday(x,dis)
    dis_index = find(x(:,23)==dis);
	Asym_number_dis_sunday = sum(x(dis_index,7));
end

function Sym_number_dis_sunday = count_Sym_dis_sunday(x,dis)
    dis_index = find(x(:,23)==dis);
	Sym_number_dis_sunday = sum(x(dis_index,8));
end

function R_number_dis_sunday = count_R_dis_sunday(x,dis)
    dis_index = find(x(:,23)==dis);
	R_number_dis_sunday = sum(x(dis_index,9));
end


function s_number = count_s(x)
    s_number = sum(x(:,4)); 
end
function dis_people = get_dis_people(x,dis,workflag)
    if workflag==1 
        dis_people = x(x(:,22)==dis,:);
    else
        dis_people = x(x(:,10)==dis,:);
    end
end

function dis_people_sunday = get_dis_people_sunday(x,dis)
    dis_people_sunday = x(x(:,23)==dis,:);
end


function exp_number = count_exp(x)
    exp_number = sum(x(:,5)); 
end
function pre_number = count_pre(x)
    pre_number = sum(x(:,6)); 
end
function asym_number = count_asym(x)
    asym_number = sum(x(:,7)); 
end
function sym_number = count_sym(x)
    sym_number = sum(x(:,8)); 
end

function r_number = count_r(x)
    r_number = sum(x(:,9)); 
end

function s_id = get_s_id(x)
    s_id =  find(x(:,4)); 
end

function exp_id = get_exp_id(x)
    exp_id =  find(x(:,5)); 
end
function pre_id = get_pre_id(x)
    pre_id =  find(x(:,6)); 
end
function asym_id = get_asym_id(x)
    asym_id =  find(x(:,7)); 
end
function sym_id = get_sym_id(x)
    sym_id =  find(x(:,8)); 
end

function r_id = get_r_id(x)
    r_id =  find(x(:,9)); 
end

function whos_state = get_state(people,person_id)
    person_index = people(:,1) == person_id;
    whos_state = people(person_index,:); 
end

function e2p_id = check_e2p(x,exp_id)
    dur_time = x(exp_id,11);
    e2p_id = exp_id(dur_time == 0);
end
function p2i_id = check_p2i(x,pre_id)
    dur_time = x(pre_id,12);
    p2i_id = pre_id(dur_time == 0);
end

function asym2r_id = check_asym2r(x,asym_id)
    dur_time = x(asym_id,13);
    asym2r_id = asym_id(dur_time == 0);
end

function sym2r_id = check_sym2r(x,sym_id)
    dur_time = x(sym_id,13);
    sym2r_id = sym_id(dur_time == 0);
end


function truncated_inf_time = truncated_poisson(mu,minvalue,size)
    truncated_inf_time = poissrnd(1/mu,size,1);
    truncated_inf_time(truncated_inf_time < minvalue) = minvalue;
end

function people = exp_happen(people,exp_uid,exp_n,n_day,epsilon,gamma,mu)
    exp_indx = [];
    for exp_count = 1:exp_n
        exp_indx = [exp_indx find(people(:,1)==exp_uid(exp_count))];
    end
    latent_time = truncated_poisson(epsilon,1,exp_n);
    pre_time = truncated_poisson(gamma,1,exp_n);
    inf_time = truncated_poisson(mu,1,exp_n);
    
    people(exp_indx,4) = 0;
    people(exp_indx,5) = 1;
    people(exp_indx,11) = latent_time;
    people(exp_indx,12) = pre_time;
    people(exp_indx,13) = inf_time;
    people(exp_indx,14) = n_day;
    people(exp_indx,15) = n_day + latent_time;
    people(exp_indx,16) = n_day + latent_time + pre_time;
    people(exp_indx,17) = n_day + latent_time + pre_time + inf_time;
end



function people = infected_happen(people,exp_id,pre_id,asym_id,sym_id,e2p_time,p2i_time,i2r_time,n_day)
    people(exp_id,5) = 1;
    people(exp_id,11) = e2p_time;
    people(exp_id,12) = p2i_time(1:length(exp_id));
    people(exp_id,13) = i2r_time(1:length(exp_id));
    people(exp_id,14) = n_day;
    people(exp_id,15) = n_day + e2p_time;
    people(exp_id,16) = n_day + e2p_time + p2i_time(1:length(exp_id));
    people(exp_id,17) = n_day + e2p_time + p2i_time(1:length(exp_id)) + i2r_time(1:length(exp_id));
    
    people(pre_id,6) = 1;
    people(pre_id,12) = p2i_time(length(exp_id)+1:length(exp_id)+length(pre_id));
    people(pre_id,13) = i2r_time(length(exp_id)+1:length(exp_id)+length(pre_id));
    people(pre_id,15) = n_day;
    people(pre_id,16) = n_day + p2i_time(length(exp_id)+1:length(exp_id)+length(pre_id));
    people(pre_id,17) = n_day + p2i_time(length(exp_id)+1:length(exp_id)+length(pre_id)) + i2r_time(length(exp_id)+1:length(exp_id)+length(pre_id));
    
    people(asym_id,7) = 1;
    people(asym_id,13) = i2r_time(length(exp_id)+length(pre_id)+1:length(exp_id)+length(pre_id)+length(asym_id));
    people(asym_id,16) = n_day;
    people(asym_id,17) = n_day + i2r_time(length(exp_id)+length(pre_id)+1:length(exp_id)+length(pre_id)+length(asym_id));
   
    people(sym_id,8) = 1;
    people(sym_id,13) = i2r_time(length(exp_id)+length(pre_id)+length(asym_id)+1:length(exp_id)+length(pre_id)+length(asym_id)+length(sym_id));
    people(sym_id,16) = n_day;
    people(sym_id,17) = n_day + i2r_time(length(exp_id)+length(pre_id)+length(asym_id)+1:length(exp_id)+length(pre_id)+length(asym_id)+length(sym_id));
end



function district_s_id = get_district_s_id(district,x)
    s_id = get_s_id(x);
    district_s_id = s_id(x(s_id,7)==district);
end


function valid_district = check_valid_dis(s_id,people)
    valid_district = unique(people(s_id,7));
end

function district_visiting = get_visit_district(district_u,dis_p,origin_dis)
    p_dist = dis_p(:,origin_dis);
    x_id = find(p_dist() >= district_u);
    district_visiting = x_id(1);
end


function dis_p = relative_p(d2d_dis,cdf_mu)
    dis_p = zeros(18);
    for i=1:18
        f = exp(-cdf_mu*d2d_dis(:,i));
        p = f/sum(f);
        dis_p(:,i) = cumsum(p);
        dis_p(end,i) = 1; 
    end
end


function [contact_index, contact_uid, infected_index, infected_uid]  = get_contact_uid(people,cum_p_contact,dis_index_pool,dist_u,contact_n,origin_dis,infectious_flag,r,beta)
    contact_index=zeros(contact_n,1); 
    contact_uid=zeros(contact_n,1);
    infected_index=[];
    infected_uid=[];
    for i=1:contact_n
        x_id = find(cum_p_contact(:,origin_dis) >= dist_u(i));
        contact_dis = x_id(1);
        contact_dis_index = dis_index_pool(contact_dis).index;
        contact_dis_uid = dis_index_pool(contact_dis).uid;
        random_index = randperm(length(contact_dis_index),1);
        contact_index(i)=contact_dis_index(random_index);
        contact_uid(i)=contact_dis_uid(random_index);
        if people(contact_index(i),4)==1
            infect_u = rand();
            if infectious_flag ==1
                if infect_u < r*beta
                    infected_index=[infected_index;contact_index(i)];
                    infected_uid=[infected_uid;contact_uid(i)];
                end
            elseif infectious_flag ==2
                if infect_u < r*beta
                    infected_index=[infected_index;contact_index(i)];
                    infected_uid=[infected_uid;contact_uid(i)];
                end                
            else
                if infect_u < beta
                    infected_index=[infected_index;contact_index(i)];
                    infected_uid=[infected_uid;contact_uid(i)];
                end
            end
        end
    end
    contact_index = unique(contact_index);
    contact_uid = unique(contact_uid);
    infected_index = unique(infected_index);
    infected_uid = unique(infected_uid);
end

