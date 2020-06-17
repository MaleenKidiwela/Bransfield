% Creates the srEvent for the Orca Tomography Shooting (inefficient looping)
% Shots at 15 m depth for tomography; 5 m depth for MCS

clear srEvent

t = readtable('ShotFiles_ZoeAtSea/orca_tomo_shotfile_final.txt');

for i = 1:size(t,1)
  srEvent.id(i,1) = table2array(t(i,1));
  srEvent.type(i,1) = 1;
  srEvent.latitude(i,1) = table2array(t(i,4));
  srEvent.longitude(i,1) = table2array(t(i,5));
  srEvent.elevation(i,1) = -0.015;
end

save srEventOrcaTomo_v1 srEvent