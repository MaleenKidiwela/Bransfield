%Script- tttables_simple.m
%modified for BRAVOSEIS Experiment. 

rTable = [0:.1:20]';%km
zTable = [.001,0.1:0.1:6]'; %kmm
%Vinterp
vp = v/1000; % Average Sounding Velocity
%15m shot depth was used for the BRAVOSEIS experiment
shotDepth=0.015; 


%travel times
for i=1:length(rTable)
    for j=1:length(zTable)
        
        ttTable(i,j) = sqrt( (rTable(i)).^2 +...
            (zTable(j)-shotDepth)^2) / vp ;
       
    end
end

%derivatives

for i=1:length(rTable)
    for j=1:length(zTable)
        
        dtdrTable(i, j) = [ (rTable(i)) ./...
            (vp*sqrt( (rTable(i)).^2 +...
            (zTable(j)-shotDepth).^2))];
        
        dtdzTable(i, j) = [(zTable(j)-0.009) ./...
            (vp*sqrt( (rTable(i)).^2+...
            (zTable(j)-shotDepth).^2 ))];           
    
    end
end


save TTtables_simple.mat ttTable dtdrTable dtdzTable rTable zTable vp


