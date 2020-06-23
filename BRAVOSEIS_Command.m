%Commands for running Tomography

srGeometry  = load_srGeometry('srGeometryOrca_v1.mat');
srStation   = load_srStation('srStationOrcaDeploy_v1.mat',srGeometry);
srEvent     = load_srEvent('srEventOrcaTomo_v1.mat',srGeometry);

tlPickDir = '/Users/earthnote/Desktop/Bransfield/Stingray/Picks'
v= 1456;
stationIn ={'BRA13','BRA14','BRA15','BRA16','BRA18','BRA19','BRA20','BRA21','BRA22','BRA23','BRA24','BRA25','BRA26','BRA27'}
%stationIn ={'BRA25'} % to isolate what you want to relocate
PhaseIn ={'Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw'}
%PhaseIn ={'Pw'}
channelIn =[]
phaseOut = {}
rlim = [0 Inf]

clear p
create_bravoseis_grd_xyz; %creates SrElevation
tlArrival = tlPick2tlArrival(srEvent, srStation, tlPickDir, stationIn, PhaseIn, phaseOut, rlim, channelIn);

TTtables_simple;
run_obsloc_bravoseis;

%% Hitogram of misfits

e= s.time-s.timePred;
e= e(:);
e=e(~isnan(e));
figure
histfit(e)
std(e)

a=sum(~isnan(s.time),2);
figure
histogram(a)


