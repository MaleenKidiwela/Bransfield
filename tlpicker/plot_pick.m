function [status,handles] = plot_pick(time, x0Trace, range, redVel, static, pickWidth, lineType,use);
%Plots one set of pick ticks or a pick line on a plot_recordSection.m record section
%
%Usage
%  [status,h] = plot_pick(time,x0Trace,redVel,static, pickWidth,lineType,use)
%
%Inputs
%  time      - Vector of pick times (1 per trace)
%  x0Trace   - Vector of pick X values (1 per trace)
%  range     - Vector of trace ranges to be used with reduction velocity
%              (empty gives default of none)
%  redVel    - Reduction velocity (0, NaN or empty for none)
%  static    - Vector of static corrections (subtracted from picks)
%              (empty gives default of none)
%  pickWidth - Width of pick ticks (Default 1)
%  lineType  - Code to specify the line type (Default is 'b-')
%              If the first character is '-' picks are plotted a line
%              joining the pick times and not as ticks
%              ticks.
%              If the last character is an interger it gives the line width
%              which is otherwise set to 0.5
%              The rest of the linetype is interpreted as the line type,
%              color and symbol by plot
%  use       - Optional input that flags picks that are used.  Picks with
%              use = false are plotted as in faint colors
%
%Outputs
%  status    - Status of exectution
%              0 - Okay
%  handles   - Vector ot tick handles
 

status = 0;
handles = [];

% Process inputs
if nargin<3
  range = [];
end
if nargin<4
  redVel = 0;
end
if nargin<5
  static = [];
end
if nargin<6
  pickWidth = 1;
end
if nargin<7
  lineType = 'b-';
end
lineType = deblank_fb(lineType);
if isempty(lineType)
  lineType = 'b-';
end
if lineType(1:1)=='-'
  oneLine = true;
  lineType = lineType(2:end);
  if isempty(lineType)
    lineType = 'b-';
  end
else
  oneLine = false;
end
if ~isnan(str2double(lineType(end)));
  lineWidth = str2double(lineType(end));
  lineType = lineType(1:end-1);
  if isempty(lineType)
    lineType = 'b-';
  end
else
  lineWidth = 1; %%%% Changed from 0.5, Nov. 23, 2016 (G. Arnoux)
end  
if nargin<8
  use = true(size(time));
end

% Apply reduction velocity
if ~isempty(redVel)
  if ~isnan(redVel)
    if redVel
      time = time - range./redVel;
    end
  end
end

% Apply a static correction
if ~isempty(static)
  if ~isnan(static)
    if any(static)
      time = time - static;
    end
  end
end

% Plotting
if oneLine
  [x0Trace, i] = sort(x0Trace);
  time = time(i);
  handles = plot(x0Trace,time,lineType,'linewidth',lineWidth);
else
  time = [time(:)'; time(:)'];
  x0Trace = [x0Trace(:)' - pickWidth/2; x0Trace(:)' + pickWidth/2];
  handles = plot(x0Trace,time, lineType,'linewidth',lineWidth);
end

if ~all(use)
  for i=1:length(handles)
    if ~use(i)
      set(handles(i),'color',[1 0.6 0.6])
    end
  end
end  
  
