% File to run the residual plots for the ETOMO experiment.  

% Inputs
%   s         - obsloc output structure (from file
%               run_obsloc_etomo_updated.m
%   ista      - Indicies of stations in s structure to plot 
%                (not station names or indicies in srStation input to obsloc)
%   doMap     - Logical 
%               0 - Make one plot of residuals versus event index
%                  color coded by station
%               1 - Make one map per station of color coded residuals
%               2 - Make one map for all stations of color coded residuals
%   doNorm    - Logical 0/1 to plot residual or residual normalized to
%               picking error
%
% Could do enhancements such as offsetting of overlapping events and 
% labeling of events 

clc;
load Obsloc_Sstructure_12km_changein_z 
%Obsloc_Sstructure_12km_wo_changein_z 


%% Input variables.
time1 = clock;
disp(['Time (' int2str(time1(1,4)) ' hr ' int2str(time1(1,5)) ' min ' int2str(time1(1,6)) ' sec)']);
   
disp(['Type the OBS station numbers you would like to run']);
ista = input(['(Example 1:10) ']);  %Input OBS numbers to be run.
for i = 1:length(ista);
    if ista(i) < 13
        ista(i) = ista(i);
    end
    if ista(i) == 13
        disp(['No OBS Information'])
    end
    if ista(i) == 41
        disp(['No OBS Information'])
    end
    if ista(i) == 15
        disp(['No OBS Information'])
    end
    if ista(i) > 13 && ista(i) <15
        ista(i) = ista(i)-1;
    else if ista(i) > 15 && ista(i) <41
            ista(i) = ista(i)-2;
        else if  ista(i) > 41
                ista(i) = ista(i)-3;
            end
        end
    end
end

disp(['Indicate whether plotting']);
disp(['(0) residuals versus event index color coded by station']);
disp(['(1) one map per station of color coded residuals']);
disp(['(2) one map for all stations of color coded residuals']);
doMap = input(['']);  
disp(['Indicate whether plotting']);
disp(['(0) residual']);
disp(['(1) residual normalized to picking error']);
doNorm = input(['']);
disp(['To locate all residual errors per specific lines, enter (1)'])
disp(['If not (0)   ']);
lin_num = input(['']);
if lin_num == 1
    line_number = input(['Which line number would you like to look at?  ']);
else
    line_number = 0;
end
plot_obsloc_event_residual(s,ista,doMap,doNorm,lin_num,line_number);
time2 = clock;
time = time2 - time1;
disp(['Amount of time to finish OBS# ' cell2mat(s.srStation.name(ista(1))) ':' cell2mat(s.srStation.name(ista(length(ista)))) ' (' int2str(time(1,4)) ' hr '  int2str(time(1,5)) ' min ' int2str(time(1,6)) ' sec)']);

    