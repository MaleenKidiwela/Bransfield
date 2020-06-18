%% plotting velocity variability - use obslocNP in run_obsloc_bravoseis > s = obslocN(p);

srGeometry  = load_srGeometry('srGeometryOrca_v1.mat');
srStation   = load_srStation('srStationOrcaDeploy_v1.mat',srGeometry);
srEvent     = load_srEvent('srEventOrcaTomo_v1.mat',srGeometry);

tlPickDir = '/Users/earthnote/Desktop/Bransfield/Stingray/Picks';
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
    cd '/Users/earthnote/Desktop/Bransfield'
    vvv
end


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

legend('BRA13','BRA14','BRA15','BRA16','BRA18','BRA19','BRA20','BRA21','BRA22','BRA23','BRA24','BRA25','BRA26','BRA27')

%%
surf([1445:1465],[13,14,15,16,18,19,20,21,22,23,24,25,26],[rms(1445:1465).a])
xlabel('velocity')
ylabel('Station')
zlabel('RMS-Norm')

plot([1445:1465],rms(1,:));
xlabel('velocity (m/s)')
ylabel('rmsNorm')
title('Normalized RMS variability with velocity')
