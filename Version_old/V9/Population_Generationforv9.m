clear all
clc
close all

% Population Generation 
%% Global varibables
global people 

%% Default parameters 
% Population parameters
pop_size = 1500;               %Hong Kong population_size : 7396000   https://www.censtatd.gov.hk/tc/web_table.html?id=216#

dis_mu = 0.5;
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
Coordinate_Y = zeros(pop_size,1);   %Coordinate_Y


people = [uid age sex I_S I_E I_P I_IA I_IS I_R district dur_exp2pre dur_pre2inf dur_inf2rec exp_time pre_time inf_time rec_time TPU TPU_edge_index Coordinate_X Coordinate_Y];


d2d_dist = csvread('numeric_DCD_dis.csv'); 

%% Read statistic data
tic
[data,text]  = xlsread('TPU_DCD_STATISTIC.xlsx');
TPU2D  = xlsread('TPU_D.xlsx');
TPU_shp = shaperead('Boundaries_of_TPU_for_2016_Population_By_Census_of_Hong_Kong.shp');
TPU_edge=[];
[TPU_shp, TPU_edge] = TPU_processing(TPU_shp,TPU_edge,TPU2D);

toc
%% Generate pop distribution

pop_tpu_stat = data;   
pop_tpu_stat(1:215,2:9) = pop_tpu_stat(1:215,2:9)./pop_tpu_stat(215,2);
pop_tpu_stat(120,1) = str2num(char(text(121,1)));
[people,pop_tpu_stat,pop_proportion] = pop_process(people,pop_tpu_stat,pop_size,TPU2D);

%% Add TPU_edge_index

people = add_tpu_edge_index(people,TPU_edge);

%% Add Coordinate
tic
coordinate_points_set =[];
[people,coordinate_points_set] = add_coordinate(people,TPU_shp,coordinate_points_set);
load('coordinate_points_set_record.mat');
coordinate_points_set_record = [coordinate_points_set_record;coordinate_points_set];
save("coordinate_points_set_record.mat",'coordinate_points_set_record')
toc

%% District Number count

N_dis = zeros(18,1);
N_dis = count_N_dis(people,N_dis);

%% Contact probability
d2d_p = relative_p(d2d_dist,dis_mu,N_dis);
p_contact = get_contact_probability(d2d_p);
cum_p_contact = cumsum(p_contact);


%% Add Coordinate

S = shaperead('HKDistrict18.shp');
T_ID = [8;7;9;17;14;1;2;3;12;13;4;18;6;5;10;11;15;16];
for i = 1:18
    S(i).ID = T_ID(i);
end
Table_S = struct2table(S);
sorted_Table_S = sortrows(Table_S,'ID');
S = table2struct(sorted_Table_S);

dis_index_pool = get_the_index_pool(people);



save("Population_Generation.mat")

function district_visiting = get_visit_district(district_u,dis_p,origin_dis)
    p_dist = dis_p(:,origin_dis);
    x_id = find(p_dist() >= district_u);
    district_visiting = x_id(1);
end

function d2d_p = relative_p(d2d_dist,cdf_mu,N_dis)
    N_matrix = repmat(N_dis,1,18);
    d2d_p = exp(-cdf_mu*d2d_dist).*N_matrix;
end

function [TPU_shp, TPU_edge] = TPU_processing(TPU_shp,TPU_edge,TPU2D)
    Table_TPU = struct2table(TPU_shp);
    sorted_Table_TPU = sortrows(Table_TPU,'TPU');
    TPU_shp = table2struct(sorted_Table_TPU);
    TPU_edge_id_all = cat(1,TPU_shp.TPU);
    for i = 1:size(TPU2D,1)
        TPU_edge_i_index = find(TPU_edge_id_all==TPU2D(i));
        TPU_edge_count_i = length(TPU_edge_i_index);
        TPU_edge_count_i_start = min(TPU_edge_i_index);
        TPU_edge_count_i_end = max(TPU_edge_i_index);
        TPU_edge_i = [TPU2D(i),TPU_edge_count_i,TPU_edge_count_i_start,TPU_edge_count_i_end];
        TPU_edge = [TPU_edge;TPU_edge_i];
    end
end

function [people,pop_tpu_stat,pop_proportion] = pop_process(people,pop_tpu_stat,pop_size,TPU2D)
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
end

function people = add_tpu_edge_index(people,TPU_edge)
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
end

function [people,coordinate_points_set] = add_coordinate(people,TPU_shp,coordinate_points_set)
    people = sortrows(people,19);

    theta = 0:pi/6:2*pi;
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    r= 0.001;

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
                     coordinate_points.TPU_shp_index = i;
                     coordinate_points.Center = [xq,yq];
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
end

function N_dis = count_N_dis(people,N_dis)
    for i=1:18
        dis_index = find(people(:,10)==i);
        N_dis(i) = length(dis_index);
    end
end

function p_contact = get_contact_probability(d2d_p)
    p_contact_un = d2d_p;
    v=sum(p_contact_un);
    D=diag(v);
    p_contact=p_contact_un*(D^-1);
end

function dis_index_pool = get_the_index_pool(people)
    for i =1:18
        dis_index_pool(i).index = find(people(:,10)==i);
        dis_index_pool(i).uid = people(dis_index_pool(i).index,1);
    end
end

function d2d_p_normal = normal_symm(d2d_p)
    A=d2d_p;
    B=triu(A);
    sum_value = sum(sum(B));
    B=B/sum_value;
    d2d_p_normal=tril(B',-1)+B;
end




