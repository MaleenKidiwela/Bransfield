%Commands for running Tomography

srGeometry  = load_srGeometry('srGeometryOrca_v1.mat');
srStation   = load_srStation('srStationOrcaDeploy_v1.mat',srGeometry);
srEvent     = load_srEvent('srEventOrcaTomo_v1.mat',srGeometry);

tlPickDir = '/Users/earthnote/Desktop/Brainsfield/Stingray/Picks'
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

%% error ellipse

srGeometry  = load_srGeometry('srGeometryOrca_v1.mat');
srStation   = load_srStation('srStationOrcaDeploy_v1.mat',srGeometry);
srEvent     = load_srEvent('srEventOrcaTomo_v1.mat',srGeometry);
tlPickDir = '/Users/earthnote/Desktop/Brainsfield/Stingray/Picks'
channelIn =[]
phaseOut = {}
rlim = [0 Inf]

stat ={'BRA13','BRA14','BRA15','BRA16','BRA18','BRA19','BRA20','BRA21','BRA22','BRA23','BRA24','BRA25','BRA26'}
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

%% plotting velocity variability

srGeometry  = load_srGeometry('srGeometryOrca_v1.mat');
srStation   = load_srStation('srStationOrcaDeploy_v1.mat',srGeometry);
srEvent     = load_srEvent('srEventOrcaTomo_v1.mat',srGeometry);

tlPickDir = '/Users/earthnote/Desktop/Brainsfield/Stingray/Picks';
stationIn ={'BRA13','BRA14','BRA15','BRA16','BRA18','BRA19','BRA20','BRA21','BRA22','BRA23','BRA24','BRA25','BRA26','BRA27'};
PhaseIn ={'Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw','Pw'};
channelIn =[];
phaseOut = {};
rlim = [0 Inf];

for vvv= 1445:1465
    v = vvv;
    create_bravoseis_grd_xyz; %creates SrElevation
    tlArrival = tlPick2tlArrival(srEvent, srStation, tlPickDir, stationIn, PhaseIn, phaseOut, rlim, channelIn);
    TTtables_simple;
    run_obsloc_bravoseis;
    rms(vvv).a=[s.rmsStationNorm];
    clearvars -except vvv rms srGeometry srStation srEvent tlPickDir stationIn PhaseIn channelIn phaseOut rlim tlArrival srElevation
    cd '/Users/earthnote/Desktop/Brainsfield'
    vvv
end
%%

surf([1445:1465],[13,14,15,16,18,19,20,21,22,23,24,25,26],[rms(1445:1465).a])
xlabel('velocity')
ylabel('Station')
zlabel('RMS-Norm')

%%
a = [1445:1465];
%rms = [rms(1445:1465).a];
for i=1:14
    b = rms(i,:);
    plot([1445:1465],b)
    hold on
    clear b idx
end

for i=1:14
    b = rms(i,:);
    idx = islocalmin(b);
    plot(a(idx),b(idx),'*r')
    text(a(idx), b(idx), num2str(a(idx)),'FontSize',9)
    hold on
    clear b idx
end

legend('BRA13','BRA14','BRA15','BRA16','BRA18','BRA19','BRA20','BRA21','BRA22','BRA23','BRA24','BRA25','BRA26')

%%
plot([1445:1465],rms(1,:));
xlabel('velocity (m/s)')
ylabel('rmsNorm')
title('Normalized RMS variability with velocity')
