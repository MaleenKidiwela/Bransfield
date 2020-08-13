%tlPick Generator

%generates a .dat files for each station in stationIn using the finalized
%srStationa and srEvent file

load('srStation_orca.mat'); %load srStation
load('srEvent_orca.mat') %load srEvent
dirName = '/Users/earthnote/Desktop/Brainsfield/Stingray/Picks_syn';%Set up directory
channelSpecific = false;
orderEventid = true;
eventID = srEvent.id;
stationIn ={'BRA13','BRA14','BRA15','BRA16','BRA18','BRA19','BRA20','BRA21','BRA22','BRA23','BRA24','BRA25','BRA26','BRA27'};
%specify which stations to make files for

for j = 1:length(stationIn)
    
    for k = 1:length(srEvent.x)
       
       dx = srStation.x(j)-srEvent.x(k);
       dy = srStation.y(j)-srEvent.y(k);
       dz = srStation.z(j)-srEvent.z(k);
       d = sqrt((dx^2)+(dy^2)+(dz^2));
       timePred(k,j) = d/1.456;
        
    end

end

for i =1:14
    
    tlPick.time = timePred(:,i);
    tlPick.eventid = srEvent.id;
    logic = isnan(tlPick.time);
    tlPick.time(logic)= [];
    tlPick.eventid(logic) =[];
    tlPick.station(1:length(tlPick.eventid)) = stationIn(i);
    tlPick.station=tlPick.station';
    tlPick.channel(1:length(tlPick.eventid))= 3;
    tlPick.channel=tlPick.channel';
    tlPick.unc(1:length(tlPick.eventid)) = 0.005;
    tlPick.unc=tlPick.unc';
    tlPick.filtLim0(1:length(tlPick.eventid)) = 5;
    tlPick.filtLim0=tlPick.filtLim0';
    tlPick.filtLim1(1:length(tlPick.eventid)) = 0;
    tlPick.filtLim1=tlPick.filtLim1';
    tlPick.filtOrder(1:length(tlPick.eventid)) = 4;
    tlPick.filtOrder=tlPick.filtOrder';
    tlPick.phase(1:length(tlPick.eventid))= {'Pw'}';
    tlPick.phase=tlPick.phase';
    tlPick.user(1:length(tlPick.eventid))= {'earthnote'};
    tlPick.user=tlPick.user';
    
    status = save_tlPick(dirName, tlPick, channelSpecific, orderEventid);
    clear tlPick
end