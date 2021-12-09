% clear all;
% clc
% close all;
rng(1)

% Population Generation 
%% Global varibables
global people 

%% Default parameters 
% Population parameters
pop_size = 1500;               %Hong Kong population_size : 7396000   https://www.censtatd.gov.hk/tc/web_table.html?id=216#
pop_infected = 2;              %round(0.1*pop_size)    Initial infected people

% Simulation parameters
start_day = datetime('today');
end_day = 0;
n_days = 10;

%NPI_impact factor

npi_factor = 1;


% Basic disease transmission parameters
R0 = 2.5;           %basic reproduction number
mu = 1/2.9;         %2.9 a rate of moving from I to R
beta = npi_factor * R0 * mu;     %a rate of moving from S to I
contact_rate =1/pop_size;
asymp_factor = 1;

dis_mu = 2;

% Efficacy of protection measures
iso_factor = 0;
quar_factor = 0;

%% Initialization Peoplemeta
% S = Susceptible       I = Infectious         R = Recovered 
% Region = Reg          Asy = Asymptomatic     Dur = Duration time
uid = randperm(pop_size)';                      % 1st
age = zeros(pop_size,1);                  % 2nd    https://www.censtatd.gov.hk/tc/web_table.html?id=216#  
%<15 1 15-24 2 25-44 3 45-64 4 65+ 5
sex = zeros(pop_size,1);                  % 3rd male 0 female 1
states_S = zeros(pop_size,1);             % S = Susceptible  4th
states_I = zeros(pop_size,1);             % I = Infectious   5th
states_R = zeros(pop_size,1);             % R = Recovered    6th
reg = zeros(pop_size,1);                  % Reg = region     7th
% https://hkplace.fandom.com/wiki/%E5%8D%81%E5%85%AB%E5%8D%80?file=Map_of_Hong_Kong_District_zh-hant.png
%Central and Western 1 Wan Chai 2 Eastern 3 Southern 4 Yau Tsim Mong 5 Sham Shui Po 6 
%Kowloon City 7 Wong Tai Sin 8 Kwun Tong 9 Kwai Tsing 10 Tsuen Wan 11
%Tuen Mun 12 Yuen Long 13 North 14 Tai Po 15 Sha Tin 16 
%Sai Kung 17 Islands 18

% TPU [121 and 123 - 124] 
dcd_dis = csvread('numeric_DCD_dis.csv'); 

dur_inf2rec = zeros(pop_size,1);          % Duration for people with infectious to recover  8th
inf_time = NaN(pop_size,1);               % inf_time = infected time  9th
rec_time = NaN(pop_size,1);               % rec_time = recovery time  10th
% day_stay_region = NaN(pop_size,1);        % work_region  11th
TPU = zeros(pop_size,1);   % TPU  11th
TPU_edge_index = zeros(pop_size,1);  % TPU_edge_index  12th
Coordinate_X = zeros(pop_size,1);   %Coordinate_X 13th
Coordinate_Y = zeros(pop_size,1);   %Coordinate_Y 14th

people = [uid age sex states_S states_I states_R reg dur_inf2rec inf_time rec_time TPU TPU_edge_index Coordinate_X Coordinate_Y];

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
people = infect_happen(people,I_ini_id,Inf_time_ini,0);
people(R_ini_id,6) = 1;

S_n = count_s(people);
I_n = count_i(people);
R_n = count_r(people);

dS_s = zeros(n_days,1);
dI_s = zeros(n_days,1);
dR_s = zeros(n_days,1);

reg_count = zeros(1,18);


%% Read statistic data

[data,text]  = xlsread('TPU_DCD_STATISTIC.xlsx');
TPU2D  = xlsread('TPU_D.xlsx');
TPU_shp = shaperead('Boundaries_of_TPU_for_2016_Population_By_Census_of_Hong_Kong.shp');

Table_TPU = struct2table(TPU_shp);
sorted_Table_TPU = sortrows(Table_TPU,'TPU');
TPU_shp = table2struct(sorted_Table_TPU);
TPU_edge=[];
TPU_edge_id_all = cat(1,TPU_shp.TPU);

for i = 1:size(TPU2D,1)
    TPU_edge_i_index = find(TPU_edge_id_all==TPU2D(i));
    TPU_edge_count_i = length(TPU_edge_i_index);
    TPU_edge_count_i_start = min(TPU_edge_i_index);
    TPU_edge_count_i_end = max(TPU_edge_i_index);
    TPU_edge_i = [TPU2D(i),TPU_edge_count_i,TPU_edge_count_i_start,TPU_edge_count_i_end];
    TPU_edge = [TPU_edge;TPU_edge_i];
end

%% Generate pop distribution
pop_tpu_stat = data;   
pop_tpu_stat(1:215,2:9) = pop_tpu_stat(1:215,2:9)./pop_tpu_stat(215,2);
pop_proportion = round(pop_tpu_stat(1:214,3:9) * pop_size);
pop_proportion(:,8) = sum(pop_proportion(:,1:2),2);
pop_proportion(:,9) = sum(pop_proportion(:,3:7),2);

sex_age_diff = pop_proportion(:,8)-pop_proportion(:,9);
age_g = find(sex_age_diff<0);
age_l = find(sex_age_diff>0);

for i = 1:length(age_g)
    age_sample = [ones(pop_proportion(age_g(i),3),1)+2; 
            ones(pop_proportion(age_g(i),4),1)+3;
            ones(pop_proportion(age_g(i),5),1)+4;
            ones(pop_proportion(age_g(i),6),1)+5;
            ones(pop_proportion(age_g(i),7),1)+6;
            ];
        
    rand_index = randperm(pop_proportion(age_g(i),9));
    draw_rand_index = rand_index(1:-sex_age_diff(age_g(i)));
    whonum = age_sample(draw_rand_index);
    pop_proportion(age_g(i),3) = pop_proportion(age_g(i),3) - sum(whonum == 3);
    pop_proportion(age_g(i),4) = pop_proportion(age_g(i),4) - sum(whonum == 4);
    pop_proportion(age_g(i),5) = pop_proportion(age_g(i),5) - sum(whonum == 5);
    pop_proportion(age_g(i),6) = pop_proportion(age_g(i),6) - sum(whonum == 6);
    pop_proportion(age_g(i),7) = pop_proportion(age_g(i),7) - sum(whonum == 7);
end

for i = 1:length(age_l)
    whonum = randi([3,7],sex_age_diff(age_l(i)),1);
    pop_proportion(age_l(i),3) = pop_proportion(age_l(i),3) + sum(whonum == 3);
    pop_proportion(age_l(i),4) = pop_proportion(age_l(i),4) + sum(whonum == 4);
    pop_proportion(age_l(i),5) = pop_proportion(age_l(i),5) + sum(whonum == 5);
    pop_proportion(age_l(i),6) = pop_proportion(age_l(i),6) + sum(whonum == 6);
    pop_proportion(age_l(i),7) = pop_proportion(age_l(i),7) + sum(whonum == 7);
end

pop_proportion(:,8) = sum(pop_proportion(:,1:2),2);
pop_proportion(:,9) = sum(pop_proportion(:,3:7),2);

age_popsize_diff  = sum(pop_proportion(:,9)) - pop_size;
if age_popsize_diff < 0
    age_d_id = randi(length(sex_age_diff),abs(age_popsize_diff),1);
    for i = 1:abs(age_popsize_diff)
        whonum = randi([3,7],1);
        pop_proportion(age_d_id(i),whonum) = pop_proportion(age_d_id(i),whonum) + 1;
    end
end

if age_popsize_diff > 0
    age_g_id = find(pop_proportion(:,9)>0);
    rand_index = randperm(length(age_g_id));
    draw_rand_index = rand_index(1:abs(age_popsize_diff));
    age_d_id = age_g_id(draw_rand_index);
    for i = 1:abs(age_popsize_diff)
        age_sample = [ones(pop_proportion(age_d_id(i),3),1)+2; 
                ones(pop_proportion(age_d_id(i),4),1)+3;
                ones(pop_proportion(age_d_id(i),5),1)+4;
                ones(pop_proportion(age_d_id(i),6),1)+5;
                ones(pop_proportion(age_d_id(i),7),1)+6;
                ];

        rand_index = randperm(pop_proportion(age_d_id(i),9));
        draw_rand_index = rand_index(1);
        whonum = age_sample(draw_rand_index);
        pop_proportion(age_d_id(i),whonum) = pop_proportion(age_d_id(i),whonum) - 1;
    end
end


pop_proportion(:,8) = sum(pop_proportion(:,1:2),2);
pop_proportion(:,9) = sum(pop_proportion(:,3:7),2);


sex_age_diff = pop_proportion(:,8)-pop_proportion(:,9);
sex_g  = find(sex_age_diff>0);
sex_l  = find(sex_age_diff<0);

for i = 1:length(sex_g)
    sex_sample = [zeros(pop_proportion(sex_g(i),1),1); 
            ones(pop_proportion(sex_g(i),2),1)
            ];
    rand_index = randperm(pop_proportion(sex_g(i),8));
    draw_rand_index = rand_index(1:sex_age_diff(sex_g(i)));
    whonum = sex_sample(draw_rand_index);
    pop_proportion(sex_g(i),1) = pop_proportion(sex_g(i),1) - sum(whonum == 0);
    pop_proportion(sex_g(i),2) = pop_proportion(sex_g(i),2) - sum(whonum == 1);
end

for i = 1:length(sex_l)
    whonum = randi([1,2],-sex_age_diff(sex_l(i)),1);
    pop_proportion(sex_l(i),1) = pop_proportion(sex_l(i),1) + sum(whonum == 1);
    pop_proportion(sex_l(i),2) = pop_proportion(sex_l(i),2) + sum(whonum == 2);
end
pop_proportion(:,8) = sum(pop_proportion(:,1:2),2);
pop_proportion(:,9) = sum(pop_proportion(:,3:7),2);





count_start = 1;

for j= 1:length(sex_age_diff)
    if pop_proportion(j,8)~=0

        sex_j = [zeros(pop_proportion(j,1),1); ones(pop_proportion(j,2),1)];
        age_j = [ones(pop_proportion(j,3),1); 
            ones(pop_proportion(j,4),1)+1;
            ones(pop_proportion(j,5),1)+2;
            ones(pop_proportion(j,6),1)+3;
            ones(pop_proportion(j,7),1)+4;
            ];
        sex_shuffle = randperm(size(sex_j,1));
        age_shuffle = randperm(size(age_j,1));
        count_end = count_start + pop_proportion(j,8) - 1;
        people(count_start:count_end,3) = sex_j(sex_shuffle);
        people(count_start:count_end,2) = age_j(age_shuffle);
        n_tpu_all = 1;
        n_tpu = floor(log10(pop_tpu_stat(j,1)))+1;

        if n_tpu == 3
            people(count_start:count_end,11) = pop_tpu_stat(j,1);
            tpu_index = find(TPU2D(:,1)==pop_tpu_stat(j,1));
            if TPU2D(tpu_index,3)==0
                people(count_start:count_end,7) = TPU2D(tpu_index,2);
            else
                reg_choice = randi(2,pop_proportion(j,8),1);
                reg_candidate = TPU2D(tpu_index,2:3);
                people(count_start:count_end,7) = reg_candidate(reg_choice);
            end
                
        else
            tpu_all = str2num(num2str(pop_tpu_stat(j,1))')';
            tpu_candidate = [];
            if tpu_all(end)==1
                n_tpu_all = (length(tpu_all)-2)/2;
                for n1 = 1:n_tpu_all
                    tpu_n1 = tpu_all(1)*100 + tpu_all(n1*2)*10 + tpu_all(n1*2+1);
                    tpu_candidate = [tpu_candidate;tpu_n1];
                end
            else
                n_tpu_all = (length(tpu_all)-3);
                for n2 = 1:n_tpu_all
                    tpu_n2 = tpu_all(1)*100 + tpu_all(2)*10 + tpu_all(n2+2);
                    tpu_candidate = [tpu_candidate;tpu_n2];
                end
            end
            tpu_randi = sort(randi(n_tpu_all,pop_proportion(j,8),1));
            people(count_start:count_end,11) = tpu_candidate(tpu_randi);
            tpu_candidate_chose = unique(tpu_candidate(tpu_randi));
            K_slice = count_start-1;
            K_end = count_start;
            for k = 1:length(tpu_candidate_chose)
                K_end = K_end + sum(people(count_start:count_end,11)==tpu_candidate_chose(k)) - 1;
                K_slice = [K_slice K_end];
                K_end = K_end + 1;
            end
            for k = 1:length(tpu_candidate_chose)
                tpu_index = find(TPU2D(:,1)==tpu_candidate_chose(k));
                if TPU2D(tpu_index,3)==0
                    people(((K_slice(k)+1):K_slice(k+1)),7) = TPU2D(tpu_index,2);
                else
                    reg_choice = randi(2,(K_slice(k+1)-K_slice(k)),1);
                    reg_candidate = TPU2D(tpu_index,2:3);
                    people(((K_slice(k)+1):K_slice(k+1)),7) = reg_candidate(reg_choice);
                end
            end

        end
        count_start = count_end + 1;
   end
end

%% Add TPU_edge_index
people = sortrows(people,11);
for i = 1:size(TPU_edge,1)
    people_TPU_index = find(people(:,11)==TPU_edge(i,1));
    if people_TPU_index
        people_TPU_index_count = length(people_TPU_index);
        people_TPU_index_start = min(people_TPU_index);
        people_TPU_index_end = max(people_TPU_index);
        index_list = randsample(TPU_edge(i,3):TPU_edge(i,4),people_TPU_index_count,true);
        people(people_TPU_index_start:people_TPU_index_end,12) = index_list;
    end
end

%% Add Coordinate
people = sortrows(people,12);

theta = 0:pi/6:2*pi;
sin_theta = sin(theta);
cos_theta = cos(theta);
r= 0.001;

coordinate_points_set = [];
coordinate_set = [];
coordinate_center_set = [];
for i = 1:length(TPU_shp)
    people_TPU_shp_index = find(people(:,12)==i);
    if people_TPU_shp_index
        people_TPU_shp_index_count = length(people_TPU_shp_index);
        people_TPU_shp_index_start = min(people_TPU_shp_index);
        people_TPU_shp_index_end = max(people_TPU_shp_index);
        n = 0;
        while (n < people_TPU_shp_index_count)
             range_xy = diff(TPU_shp(i).BoundingBox);
             rand_delta = rand(1,2);
             l_x = TPU_shp(i).BoundingBox(1,1);
             d_y = TPU_shp(i).BoundingBox(1,2);
             generated_coordinate = [l_x + rand_delta(1)*range_xy(1) d_y + rand_delta(2)*range_xy(2)];
             xq = generated_coordinate(1);
             yq = generated_coordinate(2);
             xv = TPU_shp(i).X;
             yv = TPU_shp(i).Y;
             if inpolygon(xq,yq,xv,yv)
                 n = n + 1;
                 coordinate_set = [coordinate_set;generated_coordinate];
                 coordinate_center_set = [coordinate_center_set;xq,yq];
                 coordinate_points.Geometry = 'Polygon';
                 coordinate_points.BoundingBox = TPU_shp(i).BoundingBox;
                 coordinate_points.X = xq + sin_theta*r;
                 coordinate_points.Y = yq + cos_theta*r;
%                  coordinate_points.ID = TPU_shp(i).ID;
                 coordinate_points_set = [coordinate_points_set;coordinate_points];
             end
        end
    end
end
people(:,13) = coordinate_center_set(:,1);
people(:,14) = coordinate_center_set(:,2);



%% Region dict   S_reg_dict
% All_region_dict
region_dic = [];
region_count = 1;

%S_reg_dict
s_region_dic = [];
s_region_count = 1;

for reg_id = 1:18
    reg_id_people_index = find(people(:,7)==reg_id);
    reg_count(reg_id) = length(reg_id_people_index);
    if reg_id_people_index
        region_dic(1,region_count) = reg_id;
        region_dic(2,region_count) = length(reg_id_people_index);
        region_dic(3:2 + region_dic(2,region_count),region_count) = people(people(:,7)==reg_id,1);
        region_count = region_count + 1;
        s_reg_id_people_index = find(people(people(:,7)==reg_id,4)==1);
        if s_reg_id_people_index
            s_region_dic(1,s_region_count) = reg_id;
            s_region_dic(2,s_region_count) = length(s_reg_id_people_index);
            people_s_reg_id = people(people(:,7)==reg_id,:);
            s_region_dic(3:2 + s_region_dic(2,s_region_count),s_region_count) = people_s_reg_id(s_reg_id_people_index,1);
            s_region_count = s_region_count + 1;
        else
            continue;
        end
    else
        continue
    end
    
end

%% Add Coordinate
S = shaperead('HKDistrict18.shp');
T_ID = [8;7;9;17;14;1;2;3;12;13;4;18;6;5;10;11;15;16];
for i = 1:18
    S(i).ID = T_ID(i);
end
Table_S = struct2table(S);
sorted_Table_S = sortrows(Table_S,'ID');
S = table2struct(sorted_Table_S);


%% Import cases
import_cases = csvread('Input_cases.csv');    % 1st col region  2nd col infected_start_day
import_cases_inf_start = import_cases(:,2);

% valid_reg = 0:17;

%% Run
figure
mapshow(S, 'FaceColor', 'white');
for day = 1:n_days
    tic
    i_new_index =[];
    i_new_uid = [];
    i_index = get_i_id(people);    %index
    
    i2r_index = check_i2r(people,i_index);
    i_keep_index = setdiff(i_index,i2r_index);
    
    dS_s(day) = - binornd(S_n(day), beta * I_n(day) / pop_size, 1, 1);
    
    if abs(dS_s(day)) > S_n(day)
        dS_s(day) = -S_n(day);
    end 
    
    contact_fail = 0;
    
    if -dS_s(day) > 0        
        infect_origin = randi(I_n(day),-dS_s(day),1);
        infect_dis_med = rand(-dS_s(day),1);
        
        for S_count = 1:(-dS_s(day))           
            origin_index = i_index(infect_origin(S_count));
            valid_reg = s_region_dic(1,s_region_dic(2,:) > 0); %When sus_people is very less    ==valid_reg = check_valid_dis(S_id,people)==             
            dis_origin = dcd_dis(:,people(origin_index,7));           
            max_dis = max(dis_origin);            
            [candidate_reg, candidate_reg_len]  = get_candidate_reg(max_dis,infect_dis_med(S_count),dis_origin,dis_mu);            
            reg_choosed_id = candidate_reg(randi(candidate_reg_len)); %choosed id_may not exist
            
            if ~isempty( find(valid_reg == reg_choosed_id,1))
                reg_index = find(s_region_dic(1,:) == reg_choosed_id);
                total_index = 2 + s_region_dic(2,reg_index);
                candidate_reg_uid = s_region_dic(3:total_index,reg_index);
                candidate_reg_id_len = length(candidate_reg_uid);
                infect_uid = candidate_reg_uid(randi(candidate_reg_id_len));
                
                uid_index_in_s = find(candidate_reg_uid == infect_uid) + 2;
                s_region_dic(uid_index_in_s:(total_index-1),reg_index) = s_region_dic((uid_index_in_s + 1):total_index,reg_index);
                s_region_dic(total_index,reg_index) = 0;
                s_region_dic(2,reg_index) = s_region_dic(2,reg_index) - 1;
                if s_region_dic(2,reg_index)==0
                    s_region_dic(:,reg_index) = [];
                end
                
                i_new_uid = [i_new_uid infect_uid];
                i_new_index = [i_new_index find(people(:,1) == infect_uid)];
            else
                contact_fail = contact_fail + 1;
            end

            
        end
    end
    dS_s(day) = dS_s(day) + contact_fail;
    dI_s(day) = -dS_s(day) - length(i2r_index);
    dR_s(day) = length(i2r_index);
    Inf_time = truncated_poisson(mu,1,(-dS_s(day)));

    people(i_new_index,4) = 0;
    people = infect_happen(people,i_new_index,Inf_time,day);

    people(i2r_index,5) = 0;
    people(i2r_index,6) = 1;
    people(i_keep_index,8) = people(i_keep_index,8)-1;
    
    S_n = [S_n count_s(people)];
    I_n = [I_n count_i(people)];
    R_n = [R_n count_r(people)];
    i_now_infected_index = get_i_id(people);
    new_recover_coordinate = coordinate_points_set(i2r_index);
    pause(0.5)
    hold on
    mapshow(new_recover_coordinate,'FaceColor', 'green'); 
    
    new_infected_coordinate = coordinate_points_set(i_now_infected_index);

    pause(0.5)
    hold on
    mapshow(new_infected_coordinate,'FaceColor', 'red'); 
    toc
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
    person_index = people(:,1) == person_id;
    whos_state = people(person_index,:); 
end


function i2r_id = check_i2r(x,i_id)
    dur_time = x(i_id,8);
    i2r_id = i_id(dur_time == 0);
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

%Input max_dis_of_origin, dis_media, all_dis_origin, mu 
function [candidate_reg, candidate_reg_len]  = get_candidate_reg(max_dis,infect_dis_med,dis_origin,mu)
    p_dis = relative_p(max_dis, mu);
    x_id = find(p_dis >= infect_dis_med);
    contact_dis = x_id(1)-1;
    candidate_reg = find(dis_origin == contact_dis);
    candidate_reg_len = length(candidate_reg);

end


function reg_s_id = get_reg_s_id(reg,x)
    s_id = get_s_id(x);
    reg_s_id = s_id(x(s_id,7)==reg);
end


function valid_reg = check_valid_dis(s_id,people)
    valid_reg = unique(people(s_id,7));
end

function p_n = relative_p(num_dis,cdf_mu)
    dis = 1:(num_dis+1);
    p = expcdf(dis,cdf_mu);
    diff_p = [p(1) diff(p)];
    diff_p_n = diff_p + diff_p * (1 - max(p));
    p_n = cumsum(diff_p_n);
    p_n(num_dis+1)=1;
end

