%% error ellipse BRAVOSEIS

srGeometry  = load_srGeometry('srGeometryOrca_v1.mat');
srStation   = load_srStation('srStationOrcaDeploy_v1.mat',srGeometry);
srEvent     = load_srEvent('srEventOrcaTomo_v1.mat',srGeometry);
tlPickDir = '/Users/earthnote/Desktop/Bransfield/Stingray/Picks'
channelIn =[]
phaseOut = {}
rlim = [0 Inf]

stat ={'BRA13','BRA14','BRA15','BRA16','BRA18','BRA19','BRA20','BRA21','BRA22','BRA23','BRA24','BRA25','BRA26','BRA27'}
v = 1456;

for ii = 1:13
stationIn =stat(ii); % to isolate what you want to relocate
PhaseIn ={'Pw'};  
clear p
create_bravoseis_grd_xyz; %creates SrElevation
tlArrival = tlPick2tlArrival(srEvent, srStation, tlPickDir, stationIn, PhaseIn, phaseOut, rlim, channelIn);

TTtables_simple;
run_obsloc_bravoseis;
a(ii).x1= s.a;
b(ii).y1 = s.b;


figure(1)
subplot(5,3,ii)
plot(a(ii).x1,b(ii).y1,'k')
title(stat(ii))
axis('equal')
hold on

end
