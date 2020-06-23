%% plotting velocity variability - use obslocNP in run_obsloc_bravoseis > s = obslocNP(p);

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

rms = [rms(1445:1465).a];

a = [1445:1465];
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
    minVelocity(i)= a(idx);
    clear b idx
end

legend('BRA13','BRA14','BRA15','BRA16','BRA18','BRA19','BRA20','BRA21','BRA22','BRA23','BRA24','BRA25','BRA26')

%% plotting laterial min velocity variability

log = logical(zeros(16,1))
log(1) = 1
log(6) = 1
srStation.longitude(log) = []
srStation.latitude (log) = []

figure
plot(srStation.longitude,srStation.latitude,'r*')
hold on
plot(srEvent.longitude,srEvent.latitude,'b.')
for i = 1:14
    text(srStation.longitude(i), srStation.latitude(i), num2str(minVelocity(i)),'FontSize',14)
    hold on
end

figure
[X,Y] = meshgrid(-58.7:0.005:-58.1,-62.6:0.005:-62.3);
in = griddata(srStation.longitude,srStation.latitude,[minVelocity'],X,Y)
X = X(:);
Y = Y(:);
imagesc(X,Y,in)
hold on
plot(srEvent.longitude,srEvent.latitude,'k.')
plot(srStation.longitude,srStation.latitude,'r*')
%% plotting 3d velocity plot
figure
surf([1445:1465],[13,14,15,16,18,19,20,21,22,23,24,25,26],[rms(1445:1465).a])
xlabel('velocity')
ylabel('Station')
zlabel('RMS-Norm')
%% Normalized RMS variability with velocity
plot([1445:1465],rms(1,:));
xlabel('velocity (m/s)')
ylabel('rmsNorm')
title('Normalized RMS variability with velocity')