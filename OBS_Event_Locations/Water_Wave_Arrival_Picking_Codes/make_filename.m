function filename=make_filename(depl,rcvr,chan,filetype,suffix)
% function filename=make_filename(depl,rcvr,chan,filetype,suffix)
%                                             
% Makes a filename for segy or pickfiles
% depl        - Deployment string [OBS test deployments (A,B or C)]
% rcvr        - Receiver number
% chan        - Channel number
% filetype    - 'segy' or 'pick' (default is 'segy')
% suffix      - Pick identifier (only relevant for pick files)

if nargin<4;
  filetype='segy';
  suffix=' ';
end
if nargin<5
  suffix=' ';
end
suffix=deblank(suffix);

wherestuffis
if exist('whoi_mat','var')==0; load whoi_obs.mat; end
row = whoi_mat(:,1)==rcvr;    %#ok<NODEF>
if ~any(row)
    insttype = 'soi';
else
    insttype = 'whoi';
end

if filetype(1)~='p'      
    if strcmp(insttype,'soi')==1
        if strcmp(depl,'A')==1 || strcmp(depl,'B')==1
            if chan==1; chanstr='TrilZ';
            elseif chan==2; chanstr='TrilX';
            elseif chan==3; chanstr='TrilY';
            elseif chan==4; chanstr='DPG'; 
            end
        else
            if chan==1; chanstr='L28Z';
            elseif chan==2; chanstr='L28X';
            elseif chan==3; chanstr='L28Y';
            elseif chan==4; chanstr='HYD'; 
            end  
        end
        if strcmp(depl,'A')==0 && strcmp(depl,'B')==0 && strcmp(depl,'C')==0
            depl = '';
        end
        filename=[dirsegy '/SIO_segy/OBS' int2str(rcvr) depl '_' chanstr '_etomo.segy'];
    else
        dNum = whoi_mat(row,2);
        if rcvr < 10; instNum = ['0' int2str(rcvr)]; else instNum = int2str(rcvr); end        
        if dNum < 10; dNumStr = ['0' int2str(dNum)]; else dNumStr = int2str(dNum); end
        filename=[dirsegy '/WHOI_segy/obs' instNum '_etomo_D' dNumStr '.segy'];
    end  
else
  if length(suffix)>6
    filename=suffix;
  else
    if rcvr < 10; instNum = ['0' int2str(rcvr)]; else instNum = int2str(rcvr); end
    filename=[dirpicks '/etomo_obs' instNum depl '_ch' int2str(chan) '_picks_' ...
              suffix '.mat'];  
  end
end
