clear all
clc
close all

% Population Generation 
%% Global varibables
global People 

%% Default parameters 
% Population parameters
popSize = 1500;               %Hong Kong population_size : 7396000   https://www.censtatd.gov.hk/tc/web_table.html?id=216#
d_mu = 0.5;

%% Initialization Peoplemeta

uid = randperm(popSize)';     uidCol = 1;        % UID
age = zeros(popSize,1);       ageCol = 2;        % https://www.censtatd.gov.hk/tc/web_table.html?id=216#  
%<15 1 15-24 2 25-44 3 45-64 4 65+ 5
sex = zeros(popSize,1);       sexCol = 3;        % male 0 female 1
indS = zeros(popSize,1);      indSCol = 4;        % S = Susceptible  
indE = zeros(popSize,1);      indECol = 5;        % E = Exposed   
indP = zeros(popSize,1);      indPCol = 6;        % P = Presymptomatic  
indAsym = zeros(popSize,1);   indAsymCol = 7;       % IA = Asym_Infectious   
indSym = zeros(popSize,1);    indSymCol = 8;       % IS = Sym_Infectious  
indRe = zeros(popSize,1);     indReCol = 9;        % R = Recovered/Removed    
indQ = zeros(popSize,1);      indQCol = 10;       % Q = Quarantine    
indC = zeros(popSize,1);      indCCol = 11;       % C = Comfirmation  
Dist = zeros(popSize,1);      DistCol = 12;      % 18 Districts      
% https://hkplace.fandom.com/wiki/%E5%8D%81%E5%85%AB%E5%8D%80?file=Map_of_Hong_Kong_District_zh-hant.png
%Central and Western 0 Wan Chai 2 Eastern 3 Southern 4 Yau Tsim Mong 5 Sham Shui Po 6 
%Kowloon City 7 Wong Tai Sin 8 Kwun Tong 9 Kwai Tsing 10 Tsuen Wan 11
%Tuen Mun 12 Yuen Long 13 North 14 Tai Po 15 Sha Tin 16 
%Sai Kung 17 Islands 18

durE2P = zeros(popSize,1);     durE2PCol = 13;       % Duration from Exposed to Presymptomatic 
durP2AsymOrSym = zeros(popSize,1);   durP2AsymOrSymCol = 14;     % Duration from Presymptomatic to Asym/Sym
durAsymOrSym2Re = zeros(popSize,1);  durAsymOrSym2ReCol = 15;    % Duration from Asym/Sym to Recovered
timeOfContactTracing = zeros(popSize,1);  timeOfContactTracingCol = 16;
testingDelay = zeros(popSize,1);     testingDelayCol = 17;
timeOfE = NaN(popSize,1);            timeOfECol = 18;   
timeOfP = NaN(popSize,1);            timeOfPCol = 19;   
timeOfAsymOrSym = NaN(popSize,1);    timeOfAsymOrSymCol = 20;          
timeOfRe = NaN(popSize,1);           timeOfReCol = 21;   
timeOfQ = NaN(popSize,1);            timeOfQCol = 22;
timeOfC = NaN(popSize,1);            timeOfCCol = 23;

TPU = zeros(popSize,1);              TPUCol = 24; 
idxOfTPUEdge = zeros(popSize,1);     idxOfTPUEdgeCol = 25;
xOfCoord = zeros(popSize,1);         xOfCoordCol = 26;
yOfCoord = zeros(popSize,1);         yOfCoordCol = 27;


People = [uid age sex indS indE indP indAsym indSym indRe indQ indC Dist durE2P durP2AsymOrSym durAsymOrSym2Re timeOfContactTracing testingDelay timeOfE timeOfP timeOfAsymOrSym timeOfRe timeOfQ timeOfC TPU idxOfTPUEdge xOfCoord yOfCoord];


dOfDist2Dist = csvread('numeric_DCD_dis.csv'); 

%% Read statistic data
tic
[data,text]  = xlsread('TPU_DCD_STATISTIC.xlsx');
TPU2D  = xlsread('TPU_D.xlsx');
TPUShp = shaperead('Boundaries_of_TPU_for_2016_Population_By_Census_of_Hong_Kong.shp');
TPUEdge=[];
[TPUShp, TPUEdge] = fprocessTPU(TPUShp,TPUEdge,TPU2D);
toc

%% Generate pop distribution

statsOfPopInTPU = data;   
statsOfPopInTPU(1:215,2:9) = statsOfPopInTPU(1:215,2:9)./statsOfPopInTPU(215,2);
statsOfPopInTPU(120,1) = str2num(char(text(121,1)));
[People,statsOfPopInTPU,popSizeInTPU] = fGeneratePop(People,statsOfPopInTPU,popSize,TPU2D,ageCol,sexCol,DistCol,TPUCol);

%% Add TPU_edge_index

People = fAddIdxOfTPUEdge(People,TPUEdge,TPUCol,idxOfTPUEdgeCol);

%% Add Coordinate
tic
coordinate_points_set =[];
[People,coordinate_points_set] = add_coordinate(People,TPUShp,coordinate_points_set,idxOfTPUEdgeCol,xOfCoordCol,yOfCoordCol);
load('coordinate_points_set_record.mat');
coordinate_points_set_record = [coordinate_points_set_record;coordinate_points_set];
save("coordinate_points_set_record.mat",'coordinate_points_set_record')
toc

%% District Number count

popSizeOfDists = zeros(18,1);
popSizeOfDists = fGetPopSizeOfDist(People,popSizeOfDists,DistCol);

%% Contact probability
postPOfContactRate = fGetPostPOfContactRate(dOfDist2Dist,d_mu,popSizeOfDists);
cumPostPOfContactRate = cumsum(postPOfContactRate);


%% Add Coordinate

S = shaperead('HKDistrict18.shp');
idOfTPU = [8;7;9;17;14;1;2;3;12;13;4;18;6;5;10;11;15;16];
for i = 1:18
    S(i).ID = idOfTPU(i);
end
tableS = struct2table(S);
sortedTableS = sortrows(tableS,'ID');
S = table2struct(sortedTableS);

idxAnduidOfDists = fGetidxAnduidOfDists(People,uidCol,DistCol);

save("popGeneration.mat")


function postPOfContactRate = fGetPostPOfContactRate(dOfDist2Dist,d_mu,popSizeOfDists)
    popSizeOfDistsMatrix = repmat(popSizeOfDists,1,18);
    unPostPOfContactRate = exp(-d_mu * dOfDist2Dist).*popSizeOfDistsMatrix;
    v = sum(unPostPOfContactRate);
    D = diag(v);
    postPOfContactRate = unPostPOfContactRate*(D^-1);
end

function [TPUShp, TPUEdge] = fprocessTPU(TPUShp,TPUEdge,TPU2D)
    tableTPU = struct2table(TPUShp);
    sortedTableTPU = sortrows(tableTPU,'TPU');
    TPUShp = table2struct(sortedTableTPU);
    allIdOfTPUEdge = cat(1,TPUShp.TPU);
    for i = 1:size(TPU2D,1)
        ithIdxOfTPUEdge = find(allIdOfTPUEdge==TPU2D(i));
        TPU_edge_count_i = length(ithIdxOfTPUEdge);
        TPU_edge_count_i_start = min(ithIdxOfTPUEdge);
        TPU_edge_count_i_end = max(ithIdxOfTPUEdge);
        TPU_edge_i = [TPU2D(i),TPU_edge_count_i,TPU_edge_count_i_start,TPU_edge_count_i_end];
        TPUEdge = [TPUEdge;TPU_edge_i];
    end
end

function [People,statsOfPopInTPU,popSizeInTPU] = fGeneratePop(People,statsOfPopInTPU,popSize,TPU2D,ageCol,sexCol,DistCol,TPUCol)
    popSizeInTPU = round(statsOfPopInTPU(1:214,3:9) * popSize);
    popSizeInTPU(:,8) = sum(popSizeInTPU(:,1:2),2);
    popSizeInTPU(:,9) = sum(popSizeInTPU(:,3:7),2);


    sex_age_diff = popSizeInTPU(:,8)-popSizeInTPU(:,9);
    age_g = find(sex_age_diff<0);
    age_l = find(sex_age_diff>0);

    for i = 1:length(age_g)
        age_sample = [ones(popSizeInTPU(age_g(i),3),1)+2; 
                ones(popSizeInTPU(age_g(i),4),1)+3;
                ones(popSizeInTPU(age_g(i),5),1)+4;
                ones(popSizeInTPU(age_g(i),6),1)+5;
                ones(popSizeInTPU(age_g(i),7),1)+6;
                ];

        rand_index = randperm(popSizeInTPU(age_g(i),9));
        draw_rand_index = rand_index(1:-sex_age_diff(age_g(i)));
        whonum = age_sample(draw_rand_index);
        popSizeInTPU(age_g(i),3) = popSizeInTPU(age_g(i),3) - sum(whonum == 3);
        popSizeInTPU(age_g(i),4) = popSizeInTPU(age_g(i),4) - sum(whonum == 4);
        popSizeInTPU(age_g(i),5) = popSizeInTPU(age_g(i),5) - sum(whonum == 5);
        popSizeInTPU(age_g(i),6) = popSizeInTPU(age_g(i),6) - sum(whonum == 6);
        popSizeInTPU(age_g(i),7) = popSizeInTPU(age_g(i),7) - sum(whonum == 7);
    end

    for i = 1:length(age_l)
        whonum = randi([3,7],sex_age_diff(age_l(i)),1);
        popSizeInTPU(age_l(i),3) = popSizeInTPU(age_l(i),3) + sum(whonum == 3);
        popSizeInTPU(age_l(i),4) = popSizeInTPU(age_l(i),4) + sum(whonum == 4);
        popSizeInTPU(age_l(i),5) = popSizeInTPU(age_l(i),5) + sum(whonum == 5);
        popSizeInTPU(age_l(i),6) = popSizeInTPU(age_l(i),6) + sum(whonum == 6);
        popSizeInTPU(age_l(i),7) = popSizeInTPU(age_l(i),7) + sum(whonum == 7);
    end

    popSizeInTPU(:,8) = sum(popSizeInTPU(:,1:2),2);
    popSizeInTPU(:,9) = sum(popSizeInTPU(:,3:7),2);

    age_popsize_diff  = sum(popSizeInTPU(:,9)) - popSize;
    if age_popsize_diff < 0
        age_d_id = randi(length(sex_age_diff),abs(age_popsize_diff),1);
        for i = 1:abs(age_popsize_diff)
            whonum = randi([3,7],1);
            popSizeInTPU(age_d_id(i),whonum) = popSizeInTPU(age_d_id(i),whonum) + 1;
        end
    end

    if age_popsize_diff > 0
        age_g_id = find(popSizeInTPU(:,9)>0);
        rand_index = randperm(length(age_g_id));
        draw_rand_index = rand_index(1:abs(age_popsize_diff));
        age_d_id = age_g_id(draw_rand_index);
        for i = 1:abs(age_popsize_diff)
            age_sample = [ones(popSizeInTPU(age_d_id(i),3),1)+2; 
                    ones(popSizeInTPU(age_d_id(i),4),1)+3;
                    ones(popSizeInTPU(age_d_id(i),5),1)+4;
                    ones(popSizeInTPU(age_d_id(i),6),1)+5;
                    ones(popSizeInTPU(age_d_id(i),7),1)+6;
                    ];

            rand_index = randperm(popSizeInTPU(age_d_id(i),9));
            draw_rand_index = rand_index(1);
            whonum = age_sample(draw_rand_index);
            popSizeInTPU(age_d_id(i),whonum) = popSizeInTPU(age_d_id(i),whonum) - 1;
        end
    end


    popSizeInTPU(:,8) = sum(popSizeInTPU(:,1:2),2);
    popSizeInTPU(:,9) = sum(popSizeInTPU(:,3:7),2);


    sex_age_diff = popSizeInTPU(:,8)-popSizeInTPU(:,9);
    sex_g  = find(sex_age_diff>0);
    sex_l  = find(sex_age_diff<0);

    for i = 1:length(sex_g)
        sex_sample = [zeros(popSizeInTPU(sex_g(i),1),1); 
                ones(popSizeInTPU(sex_g(i),2),1)
                ];
        rand_index = randperm(popSizeInTPU(sex_g(i),8));
        draw_rand_index = rand_index(1:sex_age_diff(sex_g(i)));
        whonum = sex_sample(draw_rand_index);
        popSizeInTPU(sex_g(i),1) = popSizeInTPU(sex_g(i),1) - sum(whonum == 0);
        popSizeInTPU(sex_g(i),2) = popSizeInTPU(sex_g(i),2) - sum(whonum == 1);
    end

    for i = 1:length(sex_l)
        whonum = randi([1,2],-sex_age_diff(sex_l(i)),1);
        popSizeInTPU(sex_l(i),1) = popSizeInTPU(sex_l(i),1) + sum(whonum == 1);
        popSizeInTPU(sex_l(i),2) = popSizeInTPU(sex_l(i),2) + sum(whonum == 2);
    end
    popSizeInTPU(:,8) = sum(popSizeInTPU(:,1:2),2);
    popSizeInTPU(:,9) = sum(popSizeInTPU(:,3:7),2);





    count_start = 1;

    for j= 1:length(sex_age_diff)
        if popSizeInTPU(j,8)~=0

            sex_j = [zeros(popSizeInTPU(j,1),1); ones(popSizeInTPU(j,2),1)];
            age_j = [ones(popSizeInTPU(j,3),1); 
                ones(popSizeInTPU(j,4),1)+1;
                ones(popSizeInTPU(j,5),1)+2;
                ones(popSizeInTPU(j,6),1)+3;
                ones(popSizeInTPU(j,7),1)+4;
                ];
            sex_shuffle = randperm(size(sex_j,1));
            age_shuffle = randperm(size(age_j,1));
            count_end = count_start + popSizeInTPU(j,8) - 1;
            People(count_start:count_end,sexCol) = sex_j(sex_shuffle);
            People(count_start:count_end,ageCol) = age_j(age_shuffle);
            n_tpu_all = 1;
            n_tpu = floor(log10(statsOfPopInTPU(j,1)))+1;

            if n_tpu == 3
                People(count_start:count_end,TPUCol) = statsOfPopInTPU(j,1);
                tpu_index = find(TPU2D(:,1)==statsOfPopInTPU(j,1));
                if TPU2D(tpu_index,3)==0
                    People(count_start:count_end,DistCol) = TPU2D(tpu_index,2);
                else
                    district_choice = randi(2,popSizeInTPU(j,8),1);
                    district_candidate = TPU2D(tpu_index,2:3);
                    People(count_start:count_end,DistCol) = district_candidate(district_choice);
                end

            else
                tpu_all = str2num(num2str(statsOfPopInTPU(j,1))')';
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
                tpu_randi = sort(randi(n_tpu_all,popSizeInTPU(j,8),1));
                People(count_start:count_end,TPUCol) = tpu_candidate(tpu_randi);
                tpu_candidate_chose = unique(tpu_candidate(tpu_randi));
                K_slice = count_start-1;
                K_end = count_start;
                for k = 1:length(tpu_candidate_chose)
                    K_end = K_end + sum(People(count_start:count_end,TPUCol)==tpu_candidate_chose(k)) - 1;
                    K_slice = [K_slice K_end];
                    K_end = K_end + 1;
                end
                for k = 1:length(tpu_candidate_chose)
                    tpu_index = find(TPU2D(:,1)==tpu_candidate_chose(k));
                    if TPU2D(tpu_index,3)==0
                        People(((K_slice(k)+1):K_slice(k+1)),DistCol) = TPU2D(tpu_index,2);
                    else
                        district_choice = randi(2,(K_slice(k+1)-K_slice(k)),1);
                        district_candidate = TPU2D(tpu_index,2:3);
                        People(((K_slice(k)+1):K_slice(k+1)),DistCol) = district_candidate(district_choice);
                    end
                end

            end
            count_start = count_end + 1;
       end
    end
end

function People = fAddIdxOfTPUEdge(People,TPUEdge,TPUCol,idxOfTPUEdgeCol)
    People = sortrows(People,TPUCol);
    for i = 1:size(TPUEdge,1)
        people_TPU_index = find(People(:,TPUCol)==TPUEdge(i,1));
        if people_TPU_index
            people_TPU_index_count = length(people_TPU_index);
            people_TPU_index_start = min(people_TPU_index);
            people_TPU_index_end = max(people_TPU_index);
            index_list = randsample(TPUEdge(i,3):TPUEdge(i,4),people_TPU_index_count,true);
            People(people_TPU_index_start:people_TPU_index_end,idxOfTPUEdgeCol) = index_list;
        end
    end
end

function [People,coordinate_points_set] = add_coordinate(People,TPUShp,coordinate_points_set,idxOfTPUEdgeCol,xOfCoordCol,yOfCoordCol)
    People = sortrows(People,idxOfTPUEdgeCol);

    theta = 0:pi/6:2*pi;
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    r= 0.001;

    coordinate_set = [];
    coordinate_center_set = [];
    for i = 1:length(TPUShp)
        people_TPU_shp_index = find(People(:,idxOfTPUEdgeCol)==i);
        if people_TPU_shp_index
            people_TPU_shp_index_count = length(people_TPU_shp_index);
            people_TPU_shp_index_start = min(people_TPU_shp_index);
            people_TPU_shp_index_end = max(people_TPU_shp_index);
            n = 0;
            while (n < people_TPU_shp_index_count)
                 range_xy = diff(TPUShp(i).BoundingBox);
                 rand_delta = rand(1,2);
                 l_x = TPUShp(i).BoundingBox(1,1);
                 d_y = TPUShp(i).BoundingBox(1,2);
                 generated_coordinate = [l_x + rand_delta(1)*range_xy(1) d_y + rand_delta(2)*range_xy(2)];
                 xq = generated_coordinate(1);
                 yq = generated_coordinate(2);
                 xv = TPUShp(i).X;
                 yv = TPUShp(i).Y;
                 if inpolygon(xq,yq,xv,yv)
                     n = n + 1;
                     coordinate_set = [coordinate_set;generated_coordinate];
                     coordinate_center_set = [coordinate_center_set;xq,yq];
                     coordinate_points.Geometry = 'Polygon';
                     coordinate_points.BoundingBox = TPUShp(i).BoundingBox;
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
    People(:,xOfCoordCol) = coordinate_center_set(:,1);
    People(:,yOfCoordCol) = coordinate_center_set(:,2);
end

function popSizeOfDists = fGetPopSizeOfDist(People,popSizeOfDists,DistCol)
    for i=1:18
        dis_index = find(People(:,DistCol)==i);
        popSizeOfDists(i) = length(dis_index);
    end
end

function idxAnduidOfDists = fGetidxAnduidOfDists(People,uidCol,DistCol)
    for i =1:18
        idxAnduidOfDists(i).index = find(People(:,DistCol)==i);
        idxAnduidOfDists(i).uid = People(idxAnduidOfDists(i).index,uidCol);
    end
end

% function d2d_p_normal = normal_symm(d2d_p)
%     A=d2d_p;
%     B=triu(A);
%     sum_value = sum(sum(B));
%     B=B/sum_value;
%     d2d_p_normal=tril(B',-1)+B;
% end




