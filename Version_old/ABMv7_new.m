% clear all;
% clc
% close all;
rng(1)

% Population Generation 
%% Global varibables
global people 

%% Default parameters 
% Population parameters
pop_size = 15000;               %Hong Kong population_size : 7396000   https://www.censtatd.gov.hk/tc/web_table.html?id=216#
exp_initial = 2;
pre_initial = 1;
asym_initial = 2;              %round(0.1*pop_size)    Initial infected people
sym_initial = 2;

% Simulation parameters
start_day = datetime('today');
end_day = 0;
n_days = 10;

%NPI_impact factor

npi_factor = 1;


% Basic disease transmission parameters
R0 = 2.5;           %basic reproduction number
mu = 1/2.9;         %2.9 a rate of moving from I to R
r = 0.55;          %ratio of transmission rate for asymptomatic over symptomatic cases
epsilon = 1/2.9;   %a rate of moving from E to P
gamma = 1/2.3;     %a rate of moving from P to I
p = 30.8/100;      %Asymptomatic rate
beta = R0*mu/(mu*r/gamma+p*r+1-p);   %a rate of moving from S to E  0.6644

N_mean_contact = 20;
Var_contact = 0.6;

dist_mu = 1;

% Efficacy of protection measures
iso_factor = 0;
quar_factor = 0;

%% Initialization Peoplemeta

uid = randperm(pop_size)';                      % 1st
age = zeros(pop_size,1);                        % 2nd    https://www.censtatd.gov.hk/tc/web_table.html?id=216#  
%<15 1 15-24 2 25-44 3 45-64 4 65+ 5
sex = zeros(pop_size,1);                  % 3rd male 0 female 1
I_S = zeros(pop_size,1);             % S = Susceptible  4th
I_E = zeros(pop_size,1);             % I = Infectious   5th
I_P = zeros(pop_size,1);             % I = Infectious   6th
I_IA = zeros(pop_size,1);            % I = Infectious   7th
I_IS = zeros(pop_size,1);            % I = Infectious   8th
I_R = zeros(pop_size,1);             % R = Recovered    9th
district = zeros(pop_size,1);        % 18 districts     10th
% https://hkplace.fandom.com/wiki/%E5%8D%81%E5%85%AB%E5%8D%80?file=Map_of_Hong_Kong_District_zh-hant.png
%Central and Western 1 Wan Chai 2 Eastern 3 Southern 4 Yau Tsim Mong 5 Sham Shui Po 6 
%Kowloon City 7 Wong Tai Sin 8 Kwun Tong 9 Kwai Tsing 10 Tsuen Wan 11
%Tuen Mun 12 Yuen Long 13 North 14 Tai Po 15 Sha Tin 16 
%Sai Kung 17 Islands 18

d2d_dis = csvread('numeric_DCD_dis.csv'); 

dur_exp2pre = zeros(pop_size,1);          % Duration for people with infectious to recover  11th
dur_pre2inf = zeros(pop_size,1);          % Duration for people with infectious to recover  12th
dur_inf2rec = zeros(pop_size,1);          % Duration for people with infectious to recover  13th
exp_time = NaN(pop_size,1);               % exp_time = exposed time   14th
pre_time = NaN(pop_size,1);               % exp_time = exposed time   15th
inf_time = NaN(pop_size,1);               % inf_time = infected time  16th
rec_time = NaN(pop_size,1);               % rec_time = recovery time  17th

TPU = zeros(pop_size,1);   % TPU  18th
TPU_edge_index = zeros(pop_size,1);  % TPU_edge_index  19th
Coordinate_X = zeros(pop_size,1);   %Coordinate_X 20th
Coordinate_Y = zeros(pop_size,1);   %Coordinate_Y 21th

people = [uid age sex I_S I_E I_P I_IA I_IS I_R district dur_exp2pre dur_pre2inf dur_inf2rec exp_time pre_time inf_time rec_time TPU TPU_edge_index Coordinate_X Coordinate_Y];

%Defined by users:
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

s_color = [55,126,184]/256;
exp_color = [152,78,163]/256;
pre_color = [255,127,0]/256;
i_color = [228,26,28]/256;
r_color = [77,175,74]/256;


district_count = zeros(1,18);


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
pop_tpu_stat(120,1) = str2num(char(text(121,1)));

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
            people(count_start:count_end,18) = pop_tpu_stat(j,1);
            tpu_index = find(TPU2D(:,1)==pop_tpu_stat(j,1));
            if TPU2D(tpu_index,3)==0
                people(count_start:count_end,10) = TPU2D(tpu_index,2);
            else
                district_choice = randi(2,pop_proportion(j,8),1);
                district_candidate = TPU2D(tpu_index,2:3);
                people(count_start:count_end,10) = district_candidate(district_choice);
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
            people(count_start:count_end,18) = tpu_candidate(tpu_randi);
            tpu_candidate_chose = unique(tpu_candidate(tpu_randi));
            K_slice = count_start-1;
            K_end = count_start;
            for k = 1:length(tpu_candidate_chose)
                K_end = K_end + sum(people(count_start:count_end,18)==tpu_candidate_chose(k)) - 1;
                K_slice = [K_slice K_end];
                K_end = K_end + 1;
            end
            for k = 1:length(tpu_candidate_chose)
                tpu_index = find(TPU2D(:,1)==tpu_candidate_chose(k));
                if TPU2D(tpu_index,3)==0
                    people(((K_slice(k)+1):K_slice(k+1)),10) = TPU2D(tpu_index,2);
                else
                    district_choice = randi(2,(K_slice(k+1)-K_slice(k)),1);
                    district_candidate = TPU2D(tpu_index,2:3);
                    people(((K_slice(k)+1):K_slice(k+1)),10) = district_candidate(district_choice);
                end
            end

        end
        count_start = count_end + 1;
   end
end

%% Add TPU_edge_index
people = sortrows(people,18);
for i = 1:size(TPU_edge,1)
    people_TPU_index = find(people(:,18)==TPU_edge(i,1));
    if people_TPU_index
        people_TPU_index_count = length(people_TPU_index);
        people_TPU_index_start = min(people_TPU_index);
        people_TPU_index_end = max(people_TPU_index);
        index_list = randsample(TPU_edge(i,3):TPU_edge(i,4),people_TPU_index_count,true);
        people(people_TPU_index_start:people_TPU_index_end,19) = index_list;
    end
end

%% Add Coordinate
people = sortrows(people,19);

theta = 0:pi/6:2*pi;
sin_theta = sin(theta);
cos_theta = cos(theta);
r= 0.001;

coordinate_points_set = [];
coordinate_set = [];
coordinate_center_set = [];
for i = 1:length(TPU_shp)
    people_TPU_shp_index = find(people(:,19)==i);
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
people(:,20) = coordinate_center_set(:,1);
people(:,21) = coordinate_center_set(:,2);



%% Region dict   S_district_dict
% All_distriction_dict
distriction_dic = [];
distriction_count = 1;

%S_district_dict
s_distriction_dic = [];
s_distriction_count = 1;

for district_id = 1:18
    district_id_people_index = find(people(:,10)==district_id);
    district_count(district_id) = length(district_id_people_index);
    if district_id_people_index
        distriction_dic(1,distriction_count) = district_id;
        distriction_dic(2,distriction_count) = length(district_id_people_index);
        distriction_dic(3:2 + distriction_dic(2,distriction_count),distriction_count) = people(people(:,10)==district_id,1);
        distriction_count = distriction_count + 1;
        s_district_id_people_index = find(people(people(:,10)==district_id,4)==1);
        if s_district_id_people_index
            s_distriction_dic(1,s_distriction_count) = district_id;
            s_distriction_dic(2,s_distriction_count) = length(s_district_id_people_index);
            people_s_district_id = people(people(:,10)==district_id,:);
            s_distriction_dic(3:2 + s_distriction_dic(2,s_distriction_count),s_distriction_count) = people_s_district_id(s_district_id_people_index,1);
            s_distriction_count = s_distriction_count + 1;
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
import_cases = csvread('Input_cases.csv');    % 1st col distriction  2nd col infected_start_day
import_cases_inf_start = import_cases(:,2);

% valid_district = 0:17;

%% Run
figure
mapshow(S, 'FaceColor', 'white');
% mapshow(coordinate_points_set(S_ini_id), 'FaceColor', s_color);
mapshow(coordinate_points_set(Exp_ini_id), 'FaceColor', exp_color);
mapshow(coordinate_points_set(Pre_ini_id), 'FaceColor', pre_color);
mapshow(coordinate_points_set(Asym_ini_id), 'FaceColor', i_color);
mapshow(coordinate_points_set(Sym_ini_id), 'FaceColor', i_color);
mapshow(coordinate_points_set(R_ini_id), 'FaceColor', r_color);
for day = 1:n_days
    tic
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
    
    N_contact = round(lognrnd(log(N_mean_contact),Var_contact,Pre_n(day) +Asym_n(day)+Sym_n(day),1));
    infectious_count = 1;
    dist_u = rand(Pre_n(day) +Asym_n(day)+Sym_n(day),1);
    for pre_count = 1:length(pre_index)
        origin_dis = people(pre_index(pre_count),10);
        dist_origin = d2d_dis(:,origin_dis);           
        max_dis = max(dist_origin);            
        [candidate_dis, candidate_dis_len]  = get_candidate_district(max_dis,dist_u(infectious_count),dist_origin,dist_mu);            
        dis_choosed = candidate_dis(randi(candidate_dis_len)); %choosed id_may not exist
        dis_people = get_dis_people(people,dis_choosed);
        if ~isempty(dis_people)
            dis_s = find(dis_people(:,4)==1);
            if dis_s
                dis_s_n = numel(dis_s);
                dis_n_contact = min(dis_s_n,N_contact(infectious_count));
                dis_exp_n = binornd(dis_n_contact, beta * r);
                if dis_exp_n > 0
                    dis_exp_indx = dis_s(randperm(dis_s_n,dis_exp_n));
                    exp_uid = dis_people(dis_exp_indx,1);
                    people = exp_happen(people,exp_uid,dis_exp_n,day,epsilon,gamma,mu);
                end                
            end
        end
        infectious_count = infectious_count + 1;
    end
    
    for asym_count = 1:length(asym_index)
        origin_dis = people(asym_index(asym_count),10);
        dist_origin = d2d_dis(:,origin_dis);           
        max_dis = max(dist_origin);            
        [candidate_dis, candidate_dis_len]  = get_candidate_district(max_dis,dist_u(infectious_count),dist_origin,dist_mu);            
        dis_choosed = candidate_dis(randi(candidate_dis_len)); %choosed id_may not exist
        dis_people = get_dis_people(people,dis_choosed);
        if ~isempty(dis_people)
            dis_s = find(dis_people(:,4)==1);
            if dis_s
                dis_s_n = numel(dis_s);
                dis_n_contact = min(dis_s_n,N_contact(infectious_count));
                dis_exp_n = binornd(dis_n_contact, beta * r);
                if dis_exp_n > 0
                    dis_exp_indx = dis_s(randperm(dis_s_n,dis_exp_n));
                    exp_uid = dis_people(dis_exp_indx,1);
                    people = exp_happen(people,exp_uid,dis_exp_n,day,epsilon,gamma,mu);
                end                
            end
        end
        infectious_count = infectious_count + 1;
    end
    
    for sym_count = 1:length(sym_index)
        origin_dis = people(sym_index(sym_count),10);
        dist_origin = d2d_dis(:,origin_dis);           
        max_dis = max(dist_origin);            
        [candidate_dis, candidate_dis_len]  = get_candidate_district(max_dis,dist_u(infectious_count),dist_origin,dist_mu);            
        dis_choosed = candidate_dis(randi(candidate_dis_len)); %choosed id_may not exist
        dis_people = get_dis_people(people,dis_choosed);
        if ~isempty(dis_people)
            dis_s = find(dis_people(:,4)==1);
            if dis_s
                dis_s_n = numel(dis_s);
                dis_n_contact = min(dis_s_n,N_contact(infectious_count));
                dis_exp_n = binornd(dis_n_contact, beta);
                if dis_exp_n > 0
                    dis_exp_indx = dis_s(randperm(dis_s_n,dis_exp_n));
                    exp_uid = dis_people(dis_exp_indx,1);
                    people = exp_happen(people,exp_uid,dis_exp_n,day,epsilon,gamma,mu);
                end                
            end
        end
        infectious_count = infectious_count + 1;
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
    
    s_now_index = get_s_id(people);
    s2e_index = setdiff(s_index, s_now_index);
    s2e_coordinate = coordinate_points_set(s2e_index);
    mapshow(s2e_coordinate,'FaceColor', exp_color); 
    
    e2p_coordinate = coordinate_points_set(e2p_index);
    mapshow(e2p_coordinate,'FaceColor', pre_color); 
    
    p2i_coordinate = coordinate_points_set(p2i_index);
    mapshow(p2i_coordinate,'FaceColor', i_color); 
    
    asym2r_coordinate = coordinate_points_set(asym2r_index);
    mapshow(asym2r_coordinate,'FaceColor', r_color); 
    
    sym2r_coordinate = coordinate_points_set(sym2r_index);
    mapshow(sym2r_coordinate,'FaceColor', r_color); 
    pause(0.5)
    toc
end


%% Definition of function

% Count Number
function s_number = count_s(x)
    s_number = sum(x(:,4)); 
end
function dis_people = get_dis_people(x,dis)
    dis_people = x(x(:,10)==dis,:);
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

%Input max_dis_of_origin, dis_media, all_dis_origin, mu 
% function [candidate_district, candidate_district_len]  = get_candidate_district(max_dis,infect_dis_med,dis_origin,mu)
%     p_dis = relative_p(max_dis, mu);
%     x_id = find(p_dis >= infect_dis_med);
%     contact_dis = x_id(1)-1;
%     candidate_district = find(dis_origin == contact_dis);
%     candidate_district_len = length(candidate_district);
% 
% end

function [candidate_district, candidate_district_len]  = get_candidate_district(max_dis,district_u,dist_origin,mu)
    p_dist = relative_p(max_dis, mu);
    x_id = find(p_dist >= district_u);
    contact_dist = x_id(1)-1;
    candidate_district = find(dist_origin == contact_dist);
    candidate_district_len = length(candidate_district);
end


function district_s_id = get_district_s_id(district,x)
    s_id = get_s_id(x);
    district_s_id = s_id(x(s_id,7)==district);
end


function valid_district = check_valid_dis(s_id,people)
    valid_district = unique(people(s_id,7));
end

function p_n = relative_p(num_dis,cdf_mu)
    dis = 1:(num_dis+1);
    p = expcdf(dis,cdf_mu);
    diff_p = [p(1) diff(p)];
    diff_p_n = diff_p + diff_p * (1 - max(p));
    p_n = cumsum(diff_p_n);
    p_n(num_dis+1)=1;
end

