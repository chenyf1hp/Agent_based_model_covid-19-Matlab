%% Add Coordinate
clc 
clear all
close all
pop_size= 1500;
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

TPU_index = 3;
rand_index = randperm(TPU_edge(TPU_index,2),1);
index_list = TPU_edge(TPU_index,3):TPU_edge(TPU_index,4);
TPU_edge_index = index_list(rand_index);


% fields = fieldnames(TPU_shp);
% fields(2:4)=[];
% fields(3)=[];
% edge = rmfield(TPU_shp,fields);


% sum_generation = sum(pop_proportion);
% if sum_generation ~= pop_size
%     randsample = randi([1,length(pop_tpu_stat(:,5))],abs(pop_size-sum_generation),1);
%     if sum_generation < pop_size
%         pop_proportion(randsample) = pop_proportion(randsample) + 1;
%     else
%         pop_proportion(randsample) = pop_proportion(randsample) - 1;
%     end
% end

% for pop_dis_id = 1:length(pop_tpu_stat(:,5))
%     if pop_proportion(pop_dis_id) == 0
%         continue;
%     else
%         up2now = sum(pop_proportion(1:pop_dis_id));
%         start_index = up2now - pop_proportion(pop_dis_id) + 1;
%         end_index = start_index + pop_proportion(pop_dis_id) - 1;
%         people(start_index:end_index,2) = pop_tpu_stat(pop_dis_id,3);
%         people(start_index:end_index,3) = pop_tpu_stat(pop_dis_id,2)+1;
%         people(start_index:end_index,7) = pop_tpu_stat(pop_dis_id,1)+1;
%     end
% end
% T_ID = [8;7;9;17;14;1;2;3;12;13;4;18;6;5;10;11;15;16];
% for i = 1:18
%     S(i).ID = T_ID(i);
% end
% Table_S = struct2table(S);
% sorted_Table_S = sortrows(Table_S,'ID');
% S = table2struct(sorted_Table_S);
% fields = fieldnames(S);
% fields(2:4)=[];
% fields(3)=[];
% edge = rmfield(S,fields);
% 
% 
% theta = 0:pi/6:2*pi;
% sin_theta = sin(theta);
% cos_theta = cos(theta);
% r= 0.001;
% 
% 
% coordinate_points_set = [];
% coordinate_set = [];
% for i = 1:18
%     if reg_count(i) > 0
%          n = 0;
%          while(n < reg_count(i))
%              range_xy = diff(edge(i).BoundingBox);
%              rand_delta = rand(1,2);
%              l_x = edge(i).BoundingBox(1,1);
%              d_y = edge(i).BoundingBox(1,2);
%              generated_coordinate = [l_x + rand_delta(1)*range_xy(1) d_y + rand_delta(2)*range_xy(2)];
%              xq = generated_coordinate(1);
%              yq = generated_coordinate(2);
%              xv = edge(i).X;
%              yv = edge(i).Y;
%              if inpolygon(xq,yq,xv,yv)
%                  n = n + 1;
%                  coordinate_set = [coordinate_set;generated_coordinate];
%                  coordinate_points.Geometry = 'Polygon';
%                  coordinate_points.BoundingBox = edge(i).BoundingBox;
%                  coordinate_points.X = xq + sin_theta*r;
%                  coordinate_points.Y = yq + cos_theta*r;
%                  coordinate_points.ID = edge(i).ID;
%                  coordinate_points_set = [coordinate_points_set;coordinate_points];
%              end
%          end
%     end
% 
% end
% generated_people = edge;