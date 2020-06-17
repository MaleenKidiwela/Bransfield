%% Examine Outliers
%load Obsloc_Sstructure_12km_wo_changein_z.mat
clear
clc
load('Obsloc_structure_stations_xy_events_xy');
misfit = s.time - s.timePred;
length_col = length(s.time(:,1));

time1 = clock;
disp(['Time: ' int2str(time1(1,4)) ' hr ' int2str(time1(1,5)) ' min ' int2str(time1(1,6)) ' sec']);
% obs_num_ent = [1:12,14,16:40,42:68];
% obs_num_ent = [19:40,42:68];
% for q=1:length(obs_num_ent);
%     obs_num_entry = obs_num_ent(q);
obs_num_entry = input(['Type the OBS you wish to look at   ']);  %Input OBS you wish to examine.

if obs_num_entry == 13 || obs_num_entry==15 || obs_num_entry==41 || obs_num_entry==45 || obs_num_entry==54
    disp(['There is no data for this OBS.'])
else if obs_num_entry == 14
        obs_num = obs_num_entry-1;
    else if obs_num_entry >15 && obs_num_entry < 41;
            obs_num = obs_num_entry-2;
        else if obs_num_entry > 41 && obs_num_entry < 45;
                obs_num = obs_num_entry-3;
            else if obs_num_entry > 45 && obs_num_entry < 54;
                    obs_num = obs_num_entry-4;
                else if obs_num_entry > 54;
                        obs_num = obs_num_entry-5;
                    else
                        obs_num = obs_num_entry;
                    end
                end
            end
        end
    end
end

% misfitentry = 0.02;
% misfitid = find(abs(misfit(:,obs_num))>misfitentry);
% event = s.srEvent.id(misfitid);

disp(['Enter (1) to Look at the outlier data points for the OBS selected']);  %Input OBS you wish to examine.
disp(['(2) to Look at a specific Event.   ']);
disp(['(3) to Look at a specific Events within a certain range.   ']);
eva = input(['']);
if eva == 1; % To Examine at outlining shots
    misfitentry = input(['What cutoff do you want for the outling data (ex. 0.02)   ']);
    misfitid = find(abs(misfit(:,obs_num))>misfitentry);
    event = s.srEvent.id(misfitid);
end
if eva ==2; % To Examine specific shots
    event_num = input(['Type the Event Number you wish to look at   ']);  %Input event id you wish to examine.
    misfitid = find(event_num == s.srEvent.id);
    event = s.srEvent.id(misfitid);
end
if eva ==3; % To Examine  shots within a certain range of the obs
    range1 = input(['Type the distance (km) around the OBS for the shots you wish to look at   ']);  %Input event id you wish to examine.
    for i = obs_num;  % For OBS in question, find shots within givin range
        dx(:,i)  = s.srStation.x(i) - s.srEvent.x;
        dy(:,i)  = s.srStation.y(i) - s.srEvent.y;
        dxy(:,i) = sqrt(dy(:,i).^2 + dx(:,i).^2); % km
        Ishot = find(dxy(:,i)<=range1);
        event = s.srEvent.id(Ishot)
    end
end

svel = 1.486;    % km/s average sounding vel for obs depth range (Using average SV of 1486 m/s from XBT data)
chan = 4; % Channel 4 is the hydrophone channel
% depl is used for the A,B,C extra sites (experimental OBSs).  Skip these
% for now.
%%%%%%%% COME BACK TO THIS %%%%%%%%
depl = '';
lat_1deg = 110.574; %km

% AR picking parameters
arord = 5;
twin = 0.8;
twin2 = 2.5;  % window for plotting data

% Filter parameters
fmini =  20;    % Hz
fmax = 80;      % Hz
dboct = 12;

%% Load OBS locations and depths
obs = read_obs_loc;
lat_deg = mean(obs.lat);   %The average latitude for the OBS's was 47.9772.
%%%%%% Make sure to set the deployment value for the experimental sites.
%%%%%%%%% We need to read this from some archival location!!! %%%%

%% Tell whether Scripps or Whoi OBS
for p = 1:length(obs.name)
    if obs.name(p) <= 10;
        instrument(p) = 1;
    elseif ((obs.name(p) >= 19) && (obs.name(p) <= 28));
        instrument(p) = 1; %'WHOI   ';
    elseif obs.name(p) == 58;
        instrument(p) = 1; %'WHOI   ';
    elseif obs.name(p) == 61;
        instrument(p) = 1; %'WHOI   ';
    else
        instrument(p) = 2; %'SCRIPPS';
    end
end


%% load the Shot lines
shot(1,64) = struct( 'number', [], 'datetime', [], 'lat',[], 'lon',[], ...
    'shiplat',[], 'shiplon',[], 'depth',[], 'linename',[]);
for i = [1:35,37:46];
    if i < 10;
        obsnum = ['0', int2str(i)];
    end
    if i >= 10;
        obsnum = int2str(i);
    end
    if i == 41;
        obsnum = '02A';
    end
    if i == 42;
        obsnum = '02R';
    end
    if i == 43;
        obsnum = '03R';
    end
    if i == 44;
        obsnum = '05A';
    end
    if i == 45;
        obsnum = '10A';
    end
    if i == 46;
        obsnum = '23R';
    end
    filename = ['/research/data/ETOMO/Data/obsip_shotlogs/MGL0910_', obsnum, '.obsip'];
    shot(i) = read_obsip_shotlog(filename);
end

waterPick(1,length(obs.name)) = struct( 'obsname', [], 'shot', [], ...
    'artime',[], 'sdv',[]);

%% Main loop over each OBS  - do parallel loop here
for i = obs_num;   %Manually typing in OBS's numbers.   %%%%length(obs.name);% For each OBS
    if  obs_num_entry == 65;        % OBS S15A is obs.name 65
        depl = 'A';
        segyname = make_filename(depl,obs.name(15),chan);  % Made the name of the SEGY file for this OBS
    else if obs_num_entry == 66;   % OBS S15B is obs.name 66
            depl = 'B';
            segyname = make_filename(depl,obs.name(15),chan);
        else if obs_num_entry == 67;   % OBS S15C is obs.name 67
                depl = 'C';
                segyname = make_filename(depl,obs.name(15),chan);
            else if obs_num_entry == 68;   % OBS S46A is obs.name 68
                    depl = 'A';
                    segyname = make_filename(depl,obs.name(46),chan);
                else
                    segyname = make_filename(depl,obs.name(obs_num_entry),chan);  % Made the name of the SEGY file for this OBS
                end, end, end, end
    
    % shotid = []; depdist =[]; tw=[];
    % For each OBS, find ALL shots within 8km range
    %     for j = 1:length(shot);    % For each shot line
    %
    %         % Distance of each shot in this line to the OBS.
    %         dlat = obs.lat(i) - shot(j).lat;
    %         dlon = obs.lon(i) - shot(j).lon;
    %         dx = dlat.*lat_1deg;
    %         dy = dlon.*lat_1deg.*cosd(lat_deg);
    %         d = (dy.^2 + dx.^2).^.5;
    
    %         % Find thdose shots that are within 8 km
    %         Ishot = find(d<=8);
    %         if isempty(Ishot), continue, end   % Skip to next itteration over j
    %         shotid = [shotid; shot(j).number(Ishot) ];
    %
    %         % Predict travel time through water column
    %         ll = sqrt( d(Ishot).^2 + (obs.depth(i)/1000).^2 );
    %         depdist = [depdist; ll ]; %convert OBS depth to km
    %     end
    
    %shotid = s.srEvent.id(find(~isnan(misfit(:,obs_num))));

    for iii = 1:length(event);
        eventid(iii,1) = find(s.srEvent.id==event(iii));
    end
    shotid = event;        
    disp(['OBS# ' int2str(obs_num_entry) ' (Channel ' int2str(chan) ') : Number of traces = ' int2str(length(shotid)) ', ' segyname ]);
    
    if instrument(i) == 1;  % Display whether WHOI or Scripps OBS
        disp('WHOI OBS');
    else
        disp('SCRIPPS OBS');
    end
    
    % Generate a window starting 1 sec before  and ending
    % one second after the smallest and largest predicted time for this
    % line
    tw = 4.7;
    tlim = [min(tw)-twin/2 max(tw)+twin/2];
    
    
    %% Load and filter the seismic data for this trace
    % read_segy.m is slow so loading all the appropriate shots for from all
    % lines for this OBS first
    
    waterPick(i).obsname = obs.name(i);
    ARtime =[];  Sdv = []; Shot = [];
    % Assign some default values for read_segy that we are not using
    tlim=[0 0];vred=0;; rlim=[0 0]; rpass=0; decfac=1;
    [seis_unfilt,hdr,reel]=read_segy(segyname,shotid,tlim,vred,rlim,rpass,decfac,chan);    %***********
    
    tcorr=[hdr(1).delay]/1e3;
    nsamp=[hdr(1).nsamp];
    sampint=[hdr(1).sampint]/1e6;   % is the same for each trace
    %        seis = seis_unfilt; % EEH When fitlering is commented out
    
    % Filter the data using parts of load_data and the
    % m-file filter_data from upicker
    %----Add a filter ramp up buffer of 2 s to seis_unfilt
    [m,n]=size(seis_unfilt);
    nfilterbuf = min(round(2/sampint)+1,m);
    seis_filterbuf1 = -seis_unfilt(nfilterbuf:-1:2,:);
    seis_filterbuf2 = -seis_unfilt(end-1:-1:end-nfilterbuf+1,:);
    for ii=1:n
        seis_filterbuf1(:,ii) = seis_filterbuf1(:,ii)+2*seis_unfilt(1,ii);
        seis_filterbuf2(:,ii) = seis_filterbuf2(:,ii)+2*seis_unfilt(end,ii);
    end
    seis_unfilt=[seis_filterbuf1; seis_unfilt; seis_filterbuf2];
    
    %----Filter if necessary
    if fmini>0 | fmax>0
        filter_data
        seis_unfilt=seis_unfilt(nfilterbuf:end-nfilterbuf+1,:);
    else
        seis=seis_unfilt(nfilterbuf:end-nfilterbuf+1,:);
        seis_unfilt=seis_unfilt(nfilterbuf:end-nfilterbuf+1,:);
        disp(' '); drawnow
    end
    % End of filtering data
    
    %% Run AR picker for each shot
    
    for k = 1:length(shotid);
        %outlier_ind(k) = find(shotid==event(k,1));
        %tw(k) = s.timePred(find(event(k,1)==s.srEvent.id),obs_num);
        tw(k) = s.time(eventid(k),obs_num);
        
        % Taken from pick_ar.m from upicker   - William Wilcock
        % Calls on the function ar_onset.m
        %AR picker code lifted directly from Doug Toomey and Rob Dunn to
        %make a pick and estimate an uncertainty.  I do not fully understand how this
        %works for two reasons
        % 1.  Since the variable AICS is 1 sample shorter than the input portion of the
        %     seismo gram, there seems to be a 1 sample uncertainty in the pick.  I have
        %     resolved this uncertainty by choosing a lag which seems to correspond
        %     better with visual picks.
        % 2.  The error estimation is a complete mystery to me.
        
        
        %----Make the AR pick
        istart = round((tw(k) - hdr(1).delay/1000 - twin/2)/sampint);
        istart2 = round((tw(k) - hdr(1).delay/1000 - twin2/2)/sampint);
        seisar = seis(istart+1:istart+twin/sampint,k);
        seisar_unfilt = seis_unfilt(istart+1:istart+twin/sampint,k);
        kst=1;
        ken=length(seisar);
        klo=kst+2*arord+1;
        khi=ken-2*arord-1;
        aics=ar_onset(seisar,kst,ken,klo,khi,arord);
        if instrument(1,obs.name) == 1 ; % Whoi Instrument
            [v,index]=min(aics);
        else if i == 60;
                [v,index]=min(aics);
            else
                [v,index]=max(diff(aics));  % Scripps Instrument
            end
        end
        index=index-1+klo;  %EEH
        %index=index+klo;
        %lagar=index-round(tpre/sampint)-1; %EEEH ??
        %----Estimate a rough std deviation - use drt subphrase
        l=aics-v<15;
        ll=find(l);
        llow=min(ll);
        lhig=max(ll);
        sdv(k)= 0.008;  % SET to 8 msec %(lhig-llow+1)*sampint/6.;
        % End of the section taken from pick_ar.m
        
        %% Plot result and evaulate
        reply =0;
        
        while reply ~= 121 & reply ~= 110,
            artime(k) =  s.time(eventid(k),obs_num);  % Water Wave previously picked;
            artimePred(k) = s.timePred(eventid(k),obs_num);  % Water Wave predicted;
            artime1(k) =  hdr(k).delay/1000 + (istart+index)*sampint;
            
            subplot(311), hold off
            plot(istart*sampint +[kst:ken]*sampint,seisar,'LineWidth',1.4)
            hold on
            plot(istart*sampint +[kst:ken]*sampint,seisar_unfilt-mean(seisar_unfilt),'--')
            yy = ylim;
            %plot([tw(k) tw(k)], yy,'g','LineWidth',1.2) % predicted time
            plot([artime(k) artime(k)],yy,'r','LineWidth',1.4)
            plot([artime(k)-sdv(k) artime(k)-sdv(k)], yy,'r--')
            plot([artime(k)+sdv(k) artime(k)+sdv(k)], yy,'r--')
            plot([artime1(k) artime1(k)],yy,'g','LineWidth',1.4)
            plot([artimePred(k) artimePred(k)],yy,'m','LineWidth',1.4)
            xax = xlim;
                        if (instrument(i) == 1);
                            ylim([(-0.5e4) (0.5e4)]);   % WHOI OBS
                        elseif (i ~= 65) && (i ~= 67);
                            ylim([(-1.5e6) (1.5e6)]);   % SCRIPPS OBS
                        end
                        if i == 35 || 60;
                            ylim([(-0.5e4) (0.5e4)]);   % SCRIPPS OBS
                        end
                        if i == 41,
                            ylim([(-1.5e4) (1.5e4)]);   % SCRIPPS OBS
                        end
                        if i == 61
                            ylim([(-1.5e6) (1.5e6)]);   % SCRIPPS OBS
                        end
                        if i == 66;
                            ylim([(-0.5e4) (0.5e4)]);
                        end
            xlabel('time, sec');
            ylabel('Amplitude');
            title([' OBS# ' int2str(obs.name(obs_num_entry)) ', Channel ' int2str(chan) ', Shot Number ' int2str(s.srEvent.id(eventid(k))) ', Shot ' int2str(k) ' of ' int2str(length(shotid)) ' Picked Time -red, Recommended Pick - green, Predicted time - mauve' ],'FontSize',11)
            
            % Plot AICS criterion
            subplot(312),hold off
            plot(istart*sampint +[1+klo:length(aics)+klo]*sampint,aics,'.-')
            hold on
            plot(istart*sampint + index*sampint, min(aics) ,'r*','MarkerSize',12)
            xlim(xax)
            ylim([(min(aics)-500) (max(aics)+500)])
            
            % Plot fitlered and unfiltered seismic trace
            subplot(313),hold off
            t = hdr(1).delay/1000 + (1:hdr(1).nsamp)*sampint;
            plot(t(istart2+1:istart2+twin2/sampint),seis(istart2+1:istart2+twin2/sampint,k))
            hold on
            yy = ylim;
            
            plot([artime(k) artime(k)],yy,'r')
            plot([artime(k)-sdv(k) artime(k)-sdv(k)], yy,'r--')
            plot([artime(k)+sdv(k) artime(k)+sdv(k)], yy,'r--')
            xlabel('time, sec');
            ylabel('Amplitude');
            
            t = hdr(1).delay/1000 + (1:hdr(1).nsamp)*sampint;
            plot(t(istart2+1:istart2+twin2/sampint),seis_unfilt(istart2+1:istart2+twin2/sampint,k) - ...
                mean(seis_unfilt(istart2+1:istart2+twin2/sampint,k)),'--')
            hold on
            yy = ylim;
            
            plot([artime(k) artime(k)],yy,'r')
            plot([artime(k)-sdv(k) artime(k)-sdv(k)], yy,'r--')
            plot([artime(k)+sdv(k) artime(k)+sdv(k)], yy,'r--')
            title(['Filtered from ' int2str(fmini) ' to ' int2str(fmax) ' and unfiltered'])
            hold off
            drawnow
            %, pause(0.1)
            %disp('Hit return to go to the next shot'), pause
            
            % The User indicates whether or not to keep this pick
            reply =0;
            while reply ~=110 &&  reply ~=121 && reply ~=1,  %110 = 'n'; 121 = 'y'' 1 = left mouse click
                disp('Do you want to keep this pick? Type y/n in the graph window.  You can also change the pick by clicking on the graph with the left mouse button.  [y]: ');
                [gx,gy,reply] = ginput(1);
                if isempty(reply) % Yes
                    reply = 121;
                elseif reply == 110 % No
                    artime(k) = nan; sdv(k) = nan;
                elseif reply == 1 % Change pick location
                    index = (gx - istart*sampint)/sampint;
                end, end, end,
        if artime(k) ~= artime1(k) && artime(k) ~= nan;
            artime(k) = artime1(k);
        end
        waterPick(obs_num_entry).shot(k) = shotid(k);
        waterPick(obs_num_entry).artime(k) = artime(k);
        waterPick(obs_num_entry).sdv(k) = sdv(k);
        %eval(['save outlier_change_waterPick_OBS_' int2str(obs_num_entry) ' waterPick'])  % Saving the results up to now to an intermediate file.  eval(['save waterPick' date ' waterPick'])
        
    end
end
%end



time2 = clock;
time = time2 - time1;
disp(['Amount of time to finish OBS# ' int2str(obs.name(obs_num_entry)) ': ' int2str(time(1,4)) ' hr ' int2str(time(1,5)) ' min ' int2str(time(1,6)) ' sec']);

%
% load Obsloc_Sstructure_12km_changein_z
% n = [4,5,11,12,13,15,17,30,35,36,41,49,50,56,57,64];
% obs_num = 68;
% %s = nan(length_col,68);
% for j = 1:obs_num  %length(s.xStation);
%     if j ~= n;
%         file_name = ['/Volumes/research/users/awells3/OBS Location/Outlier Data/outlier_change_waterPick_OBS_', int2str(j), '.mat'];
%         load(file_name);
%         if length(waterPick(j).shot)>1
%             for k = 1:length(waterPick(j).shot);
%                 ind(k) = find(waterPick(j).shot(k) == s.srEvent.id);
%             end
%         else
%         ind = find(waterPick(j).shot == s.srEvent.id);
%         end
%
%
%         obs_num_ent = [1:12,14,16:40,42:68];
% for q=1:length(obs_num_ent);
%     obs_num_entry = obs_num_ent(q);
%

%obs_num_entry = input(['Type the OBS you wish to look at   ']);  %Input OBS you wish to examine.
% if obs_num_entry == 13 || obs_num_entry==15 || obs_num_entry==41
%         disp(['There is no data for this OBS.'])
% else if obs_num_entry == 14
%         obs_num = obs_num_entry-1;
% else if obs_num_entry >15 && obs_num_entry < 41;
%         obs_num = obs_num_entry-2;
% else if obs_num_entry > 41;
%         obs_num = obs_num_entry-3;
%     else
%         obs_num = obs_num_entry;
%     end
%     end
%     end
%         s.time(ind,j) = waterPick(j).artime;
%
%
%
%     end
%     clear waterPick ind
% end
%
%
%
% save Obsloc_Sstructure_12km_changein_z_minus_outliers s



