% Create a srStation for deployed seismometer locations at Orca volcano

clear srStation

t = readtable('orcaDeploy.csv');

for i = 1:size(t,1)
  srStation.name(i,1) =  table2cell(t(i,1));
  srStation.longitude(i,1) = table2array(t(i,2));
  srStation.latitude(i,1) = table2array(t(i,3));
  srStation.elevation(i,1) = -table2array(t(i,4))/1000;
end

save srStationOrcaDeploy_v1 srStation