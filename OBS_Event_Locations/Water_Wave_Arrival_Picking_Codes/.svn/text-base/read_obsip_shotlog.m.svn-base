 function shot = read_obsip_shotlog(filename)

% Function to read the obsip_shotlog file and return the data below header 
% in the structure shot with the fields.  NOTE that this log does not
% record the gun depth.
%
% INPUTS - the name of the obsip file to read e.g 'MGL0910_09.obsip'
%
% OUTPUTS
% SHOT structure with the fields
%   NUMBER      - the shot number
%   DATETIME    - the date time string formatted YEAR-MONTH-DAY HR:MIN:SEC
%   LAT         - the latitude of the shot
%   LON         - the longitude of the shot
%   SHIPLAT     - the latitude of the ship
%   SHIPLON     - the longitude of the ship
%   DEPTH       - the water depth of the ship at the shot point
%   LINENAME    - cruise line name MGL0910_LINE#


%file = input('Input file name : ','s');

fid = fopen(filename,'r');

header1 = fgetl(fid);        %#obsipshotfile v1.0 
header2 = fgetl(fid);        %#shotnumber date time sourceLat sourceLon shipLat shipLon waterDepth sciTag

c = fread(fid,'uchar');
c = setstr(c');
if length(filename) == 58;
c = reshape(c,101,length(c)/101)';
end
if length(filename) == 57;
c = reshape(c,100,length(c)/100)';
end

shot.number     = str2num(c(:,1:7));
shot.datetime   = c(:,8:34);
shot.lat        = str2num(c(:,36:45));
shot.lon        = str2num(c(:,47:58));
shot.shiplat    = str2num(c(:,60:68));
shot.shiplon    = str2num(c(:,70:80));
shot.depth      = str2num(c(:,82:88));
shot.linename   = c(:,89:99);

fclose(fid);

