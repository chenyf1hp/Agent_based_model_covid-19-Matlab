clear all
clc
close all
S = shaperead('Boundaries_of_TPU_for_2016_Population_By_Census_of_Hong_Kong.shp');

figure
mapshow(S, 'FaceColor', 'white');

