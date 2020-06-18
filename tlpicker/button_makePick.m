% Script to make picks
%
% 7/11/12 - Fixed 'Q' option so initial plot parameters are restored
% 1/9/20  - Augmented try/catch segments so that invalid values do not replace original values 

% Flag a ginput problem for Mac keyboard entries (x,y are wrong!)
% Doesn't appear to be a problem in Matlab R2015a or later versions
ginputBug = false;

% Turn off menu
menu_off
drawnow

% Need to have loaded SEGY
if ~doneWhat.segy
  beep
  warning('Failure - No SEGY loaded (need to load SEGY, plot a record section and load picks)')

% and to have plotted it
elseif ~doneWhat.plot
  beep
  warning('Failure - No record section plotted (need to plot a record section and load picks)')
  
else
  % If picks have not been loaded they will be loaded
  if ~doneWhat.loadPick
    button_loadPick
  end
  
  % If picks are successfully loaded do the picking
  if doneWhat.loadPick
    currentMenuKeep = currentMenu;
    
    % Print picking options
    
    disp(' ')
    disp('  *** PICKING OPTIONS ***')
    disp('  Left Mouse   (Click)       - Make a pick (nearest trace to click)');
    disp('  Middle Mouse (Shift-click) - Assign uncertainty of menu picking parameter unc to an existing pick');
    disp('  Right Mouse  (Ctrl-Click)  - Delete a pick');
    disp('  a - Return the plot to the initial plot area (i.e., X limits and time limits when MAKE PICKS was selected)');
    disp('  b - Rubberband a new plot area');
    disp('  c - Change the menu plotting parameter clip and replot the record section');
    disp('  d - Set axes to default values (i.e. set the X-axis limits and Time Limits menu plotting parameters to 0)');
    disp('  e - Change the menu picking parameter unc');
    disp('  f - Rubberband X limits to assign menu uncertainty to a set of existing picks');
    disp('  g - Rubberband X limits to delete a set of existing picks');
    disp('  h - Assign uncertainty of menu picking parameter unc to an existing pick*');
    disp('  i - Return all the menu parameters to their initial values');
    disp('  j - Change the menu plotting parameter Label Inc');
    disp('  k - Change the menu picking parameter Comment');
    disp('  l - Change the menu plotting parameter Label Option');
    disp('  m - Mouse in picks by clicking points on a line with picks interpolated linearly inbetween');
    disp('  n - Recenter plot about current mouse click*');
    disp('  o - Zoom out about cursor click by a factor of x 2*');
    disp('  p - Replot the record section - useful if one changes parameters directly on the menu');
    disp('  q - Quit picking and return all the menu parameters to their starting values');
    disp('  r - Reverse the X-axes');
    disp('  s - Change the menu plotting parameter X scale');
    disp('  t - Change the menu plotting parameter X Scale Option');
    disp('  u - Change the menu picking parameter User');
    disp('  v - Change the menu plotting parameter Reduction Velocity');
    disp('  w - Change the menu plotting parameter Wiggle');
    disp('  x - Quit picking and and leave all the menu parameters at their current values');
    disp('  y - Not Used');
    disp('  z - Zoom in about cursor click by a factor of x 0.5*');
    disp('  - - Increase the picking error of the selected pick by approximately sqrt(2)*');
    disp('  + - Decrease the picking error of the selected pick by approximately 1/sqrt(2)*');
    disp('  / - Toggle the "use" status of an existing pick*')
    disp('  [ - Rubberband X limits to assign "use" status to true')
    disp('  ] - Rubberband X limits to assign "use" status to false')
    disp('  1 - Move area to the left a half screen');
    disp('  2 - Move area to the right a half screen');
    disp('  9 - Move area down a half screen');
    disp('  0 - Move area up a half screen');
    disp('  *Options requires a 2nd mouse click if ginputBug = true')


    % Graphical picking options here in a while loop with switch/case/otherwise;
    moreToDo = true;
    while moreToDo

      % Mouse input
      figure(hFig)
      button = [];
      while isempty(button)
        [xm,ym,button] = ginput(1);
        if isempty(button)
          beep
          warning('Warning : Mouse entry with carriage return is not allowed')
        elseif debugPicking
          ginputRecord = [ginputRecord; [xm ym button]];
        end
      end
      if ginputBug
        if button==abs('-') || button==abs('+') || button==abs('o') || ...
           button==abs('z') || button==abs('n') || button==abs('h') || ...
           button==abs('/')  %% Must apply only to certain options
          disp('This option requires a mouse entry because of the ginput keyboard x-y are unreliable')
          xm = [];
          while isempty(xm)
            [xm,ym,test] = ginput(1);
            if isempty(xm)
              beep
              warning('Mouse entry with carriage return is not allowed')
            elseif test~=1 && test~=2 && test~=3
              beep
              warning('Keyboard Entry is not allowed')   
            elseif debugPicking
              ginputRecord = [ginputRecord; [xm ym button]];
            end
          end

        end
      end

      % Index of picked trace
      [~, i] = min(abs(x0Trace-xm));

      % Remove static and reduction velocity from pick
      ymCor = ym + traceMetaData.static(i);
      if ~isnan(lastPlot.redVel)
        if lastPlot.redVel
          ymCor = ymCor + abs(traceMetaData.range(i))./(lastPlot.redVel);
        end
      end

      % Old pick plotted y
      ymOld = active.time + traceMetaData.static;
      if ~isnan(lastPlot.redVel)
        if lastPlot.redVel
          ymOld = ymOld - abs(traceMetaData.range)./(lastPlot.redVel);
        end
      end

      switch button

        % Make a pick
        case 1
          if isgraphics(handlesActive.time(i))
            delete(handlesActive.time(i))
            delete(handlesActive.unc1(i))
            delete(handlesActive.unc2(i))
          end
          handlesActive.time(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ym ym],'r-','linewidth',2);
          handlesActive.unc1(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ym ym]-currentMenu.active.unc,'r:');
          handlesActive.unc2(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ym ym]+currentMenu.active.unc,'r:');
          active.channel(i) = traceMetaData.channel(i);
          active.time(i) = ymCor;
          active.unc(i) = currentMenu.active.unc;
          active.phase(i) = {currentMenu.active.phase};
          if doneWhat.filter
            active.filtLim0(i) = lastFilter.lim0;
            active.filtLim1(i) = lastFilter.lim1;
            active.filtOrder(i) = lastFilter.order;
            active.filtZeroPhase(i) = lastFilter.zeroPhase;
          else
            active.filtLim0(i) = 0;
            active.filtLim1(i) = 0;
            active.filtOrder(i) = 0;
            active.filtZeroPhase(i) = false;
          end          
          active.scale(i) = absScale(i);
          active.user(i) = {currentMenu.active.user};
          c = clock;
          active.lddate(i) = date2secnds(c(1),julday(c(2),c(3),c(1)),c(4),c(5),c(6));
          active.use(i) = true;
          active.comment(i) = {currentMenu.active.comment};
          active.picked(i) = true;
          active.updated(i) = true;

        % Assign an error
        case {2, abs('h')}
          if active.picked(i)
            active.unc(i) = currentMenu.active.unc;
            delete(handlesActive.unc1(i))
            delete(handlesActive.unc2(i))
            handlesActive.unc1(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ymOld(i) ymOld(i)]-currentMenu.active.unc,'r:');
            handlesActive.unc2(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ymOld(i) ymOld(i)]+currentMenu.active.unc,'r:');
            active.updated(i) = true;
          else
            warning('Cannot assign a pick error to an unpicked trace')
          end

        % Delete a pick
        case 3
          if active.picked(i)
            ym = NaN;
            if isgraphics(handlesActive.time(i))
              delete(handlesActive.time(i))
              delete(handlesActive.unc1(i))
              delete(handlesActive.unc2(i))
%               handlesActive.time(i) = NaN;
%               handlesActive.unc1(i) = NaN;
%               handlesActive.unc2(i) = NaN;
            end
            active.time(i) = NaN;
            active.picked(i) = false;
            active.updated(i) = true;
          end

        % Return to initial plot area
        case abs('a')
          currentMenu.plot.xlim0 = currentMenuKeep.plot.xlim0;
          currentMenu.plot.xlim1 = currentMenuKeep.plot.xlim1;
          currentMenu.plot.tlim0 = currentMenuKeep.plot.tlim0;
          currentMenu.plot.tlim1 = currentMenuKeep.plot.tlim1;
          set(handlesMenu.plot_xlim0,'string',num2str(currentMenuKeep.plot.xlim0));
          set(handlesMenu.plot_xlim1,'string',num2str(currentMenuKeep.plot.xlim1));
          set(handlesMenu.plot_tlim0,'string',num2str(currentMenuKeep.plot.tlim0));
          set(handlesMenu.plot_tlim1,'string',num2str(currentMenuKeep.plot.tlim1));
          upToDate.plot = false;
          button_plot;

        % Rubber band a new plot area  
        case abs('b')
          disp('Rubberband in a new plot area')
          [xlim0, xlim1, tlim0, tlim1] = rubberband(hFig);
          if xlim0==xlim1 || tlim0==tlim1
            warning('Rubberband area  must have different minimum & maximum values')
          else
            currentMenu.plot.xlim0 = xlim0;
            currentMenu.plot.xlim1 = xlim1;
            currentMenu.plot.tlim0 = tlim0;
            currentMenu.plot.tlim1 = tlim1;
            set(handlesMenu.plot_xlim0,'string',num2str(currentMenu.plot.xlim0));
            set(handlesMenu.plot_xlim1,'string',num2str(currentMenu.plot.xlim1));
            set(handlesMenu.plot_tlim0,'string',num2str(currentMenu.plot.tlim0));
            set(handlesMenu.plot_tlim1,'string',num2str(currentMenu.plot.tlim1));
            upToDate.plot = false;
            button_plot;
          end        

        % Update clip parameter  
        case abs('c')
          disp('Enter a new value of the clip parameter and hit return')
          [~,~,string] = ginput;
          tempNum = currentMenu.plot.clip;
          try
            eval(['currentMenu.plot.clip = ' char(string') ';']);
            set(handlesMenu.plot_clip,'string',num2str(currentMenu.plot.clip));
            upToDate.plot = false;
            button_plot
          catch
            warning(['Clip parameter "' char(string') '" is invalid'])
            eval(['currentMenu.plot.clip = ' num2str(tempNum) ';']);
            set(handlesMenu.plot_clip,'string',num2str(currentMenu.plot.clip));
          end

        % Set axes to default values  
        case abs('d')
            currentMenu.plot.xlim0 = 0;
            currentMenu.plot.xlim1 = 0;
            currentMenu.plot.tlim0 = 0;
            currentMenu.plot.tlim1 = 0;
            set(handlesMenu.plot_xlim0,'string',num2str(currentMenu.plot.xlim0));
            set(handlesMenu.plot_xlim1,'string',num2str(currentMenu.plot.xlim1));
            set(handlesMenu.plot_tlim0,'string',num2str(currentMenu.plot.tlim0));
            set(handlesMenu.plot_tlim1,'string',num2str(currentMenu.plot.tlim1));
            upToDate.plot = false;
            button_plot;

        % Update pick uncertainty (error)    
        case abs('e')
          disp('Enter a new value of the picking uncertainty (error)')
          [~,~,string] = ginput;
          tempNum = currentMenu.active.unc;
          try
            eval(['currentMenu.active.unc = ' char(string') ';']);
            set(handlesMenu.active_unc,'string',num2str(currentMenu.active.unc));
          catch
            warning(['Picking uncertainty (error) "' char(string') '" is invalid'])
            eval(['currentMenu.active.unc = ' num2str(tempNum) ';']);
            set(handlesMenu.active_unc,'string',num2str(currentMenu.active.unc));
          end

        %Rubberband setting of pick uncertainty
        case abs('f')
          disp('Rubberband in X-limits to set pick uncertainty')
          [xlim0, xlim1] = rubberband(hFig);
          index = find(x0Trace>=xlim0 & x0Trace<=xlim1);
          for i = index
            if active.picked(i)
              active.unc(i) = currentMenu.active.unc;
              delete(handlesActive.unc1(i))
              delete(handlesActive.unc2(i))
              handlesActive.unc1(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ymOld(i) ymOld(i)]-currentMenu.active.unc,'r:');
              handlesActive.unc2(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ymOld(i) ymOld(i)]+currentMenu.active.unc,'r:');
              active.updated(i) = true;
            end
          end

        % Delete picks by rubber banding
        case abs('g')
          disp('Rubberband in X-limits to delete picks')
          [xlim0, xlim1] = rubberband(hFig);
          index = find(x0Trace>=xlim0 & x0Trace<=xlim1);
          for i = index
            if active.picked(i)
              ym = NaN;
              if isgraphics(handlesActive.time(i))
                delete(handlesActive.time(i))
                delete(handlesActive.unc1(i))
                delete(handlesActive.unc2(i))
%                 handlesActive.time(i) = NaN;
%                 handlesActive.unc1(i) = NaN;
%                 handlesActive.unc2(i) = NaN;
              end
              active.time(i) = NaN;
              active.picked(i) = false;
              active.updated(i) = true;
            end
          end

        % Return to initial menu parameters
        case abs('i')
          currentMenu = currentMenuKeep;
          upToDate.plot = false;
          button_plot

        % Adjust the label increment
        case abs('j')
          disp('Enter a new value of the Label Increment and hit return')
          [~,~,string] = ginput;
          tempNum = currentMenu.plot.labelIncrement;
          try
            eval(['currentMenu.plot.labelIncrement = ' char(string') ';']);
            set(handlesMenu.plot_labelIncrement,'string',num2str(currentMenu.plot.labelIncrement));
            upToDate.plot = false;
            button_plot
          catch
            warning(['Label Increment "' char(string') '" is invalid'])
            eval(['currentMenu.plot.labelIncrement = ' num2str(tempNum) ';']);
            set(handlesMenu.plot_labelIncrement,'string',num2str(currentMenu.plot.labelIncrement));
          end

        % Change the comment
        case abs('k')
          disp('Enter a new comment and hit return')
          [~,~,string] = ginput;
          currentMenu.active.comment = char(string');
          set(handlesMenu.active_comment,'string',currentMenu.active.comment);
          upToDate.plot = false;
          button_plot
 
        % Change the label option  
        case abs('l')
          disp('Enter a new value of the Label Option parameter and hit return')
          [~,~,string] = ginput;
          tempString = currentMenu.plot.labelOption;
          try
            currentMenu.plot.labelOption = char(string');
            set(handlesMenu.plot_labelOption,'string',currentMenu.plot.labelOption);
            upToDate.plot = false;
            button_plot
          catch
            warning(['Label Option parameter "' char(string') '" is invalid'])
            currentMenu.plot.labelOption = tempString;
            set(handlesMenu.plot_labelOption,'string',currentMenu.plot.labelOption);
          end

        case abs('m')
          disp('Mouse in pick lines (interpolation between multiple plots) and hit return')
          [xmVec, ymVec] = ginput;
          [xmVec,i] = sort(xmVec);
          ymVec = ymVec(i);
          if length(xmVec)<2 
            warning('Need to enter at least two points to interpolate picks')
          elseif ~all(diff(xmVec))
            warning('Points must have distict X-values for pick interpolation')
          else
            index = find(x0Trace>=xmVec(1) & x0Trace<=xmVec(end));
            for i = index
              if isgraphics(handlesActive.time(i))
                delete(handlesActive.time(i))
                delete(handlesActive.unc1(i))
                delete(handlesActive.unc2(i))
              end
              ym = interp1(xmVec, ymVec, x0Trace(i), 'spline');
              ymCor = ym + traceMetaData.static(i);
              if ~isnan(lastPlot.redVel)
                if lastPlot.redVel
                  ymCor = ymCor + abs(traceMetaData.range(i))./(lastPlot.redVel);
                end
              end
              handlesActive.time(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ym ym],'r-','linewidth',2);
              handlesActive.unc1(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ym ym]-currentMenu.active.unc,'r:');
              handlesActive.unc2(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ym ym]+currentMenu.active.unc,'r:');
              active.channel(i) = traceMetaData.channel(i);
              active.time(i) = ymCor;
              active.unc(i) = currentMenu.active.unc;
              active.phase(i) = {currentMenu.active.phase};
              if doneWhat.filter
                active.filtLim0(i) = lastFilter.lim0;
                active.filtLim1(i) = lastFilter.lim1;
                active.filtOrder(i) = lastFilter.order;
                active.filtZeroPhase(i) = lastFilter.zeroPhase;
              else
                active.filtLim0(i) = 0;
                active.filtLim1(i) = 0;
                active.filtOrder(i) = 0;
                active.filtZeroPhase(i) = false;
              end          
              active.scale(i) = absScale(i);
              active.user(i) = {currentMenu.active.user};
              c = clock;
              active.lddate(i) = date2secnds(c(1),julday(c(2),c(3),c(1)),c(4),c(5),c(6));
              active.use(i) = true;
              active.comment(i) = {currentMenu.active.comment};
              active.picked(i) = true;
              active.updated(i) = true;
            end
          end

        % Recenter plot     
        case abs('n')
          temp = xlim;
          temp = xm + [-0.5 0.5] * diff(temp);
          currentMenu.plot.xlim0 = temp(1);
          currentMenu.plot.xlim1 = temp(2);
          temp = ylim;
          temp = ym + [-0.5 0.5] * diff(temp);
          currentMenu.plot.tlim0 = temp(1);
          currentMenu.plot.tlim1 = temp(2);
          set(handlesMenu.plot_xlim0,'string',num2str(currentMenu.plot.xlim0));
          set(handlesMenu.plot_xlim1,'string',num2str(currentMenu.plot.xlim1));
          set(handlesMenu.plot_tlim0,'string',num2str(currentMenu.plot.tlim0));
          set(handlesMenu.plot_tlim1,'string',num2str(currentMenu.plot.tlim1));
          upToDate.plot = false;
          button_plot;

        % Zoom out x2     
        case abs('o')
          temp = xlim;
          temp = xm + [-1 1] * diff(temp);
          currentMenu.plot.xlim0 = temp(1);
          currentMenu.plot.xlim1 = temp(2);
          temp = ylim;
          temp = ym + [-1 1] * diff(temp);
          currentMenu.plot.tlim0 = temp(1);
          currentMenu.plot.tlim1 = temp(2);
          set(handlesMenu.plot_xlim0,'string',num2str(currentMenu.plot.xlim0));
          set(handlesMenu.plot_xlim1,'string',num2str(currentMenu.plot.xlim1));
          set(handlesMenu.plot_tlim0,'string',num2str(currentMenu.plot.tlim0));
          set(handlesMenu.plot_tlim1,'string',num2str(currentMenu.plot.tlim1));
          upToDate.plot = false;
          button_plot;

        % Plot Record Section again  
        case abs('p')
          upToDate.plot = false;
          button_plot

        % Quit picking  
        case abs('q')
          moreToDo = false;
          currentMenu = currentMenuKeep;
          set(handlesMenu.plot_xOption,'string',currentMenu.plot.xOption);
          set(handlesMenu.plot_xlim0,'string',num2str(currentMenu.plot.xlim0));
          set(handlesMenu.plot_xlim1,'string',num2str(currentMenu.plot.xlim1));
          set(handlesMenu.plot_xScaleOption,'string',currentMenu.plot.xScaleOption);
          set(handlesMenu.plot_xScale,'string',num2str(currentMenu.plot.xScale));
          set(handlesMenu.plot_clip,'string',num2str(currentMenu.plot.clip));
          set(handlesMenu.plot_wiggleOption,'string',num2str(currentMenu.plot.wiggleOption));
          set(handlesMenu.plot_tlim0,'string',num2str(currentMenu.plot.tlim0));
          set(handlesMenu.plot_tlim1,'string',num2str(currentMenu.plot.tlim1));
          set(handlesMenu.plot_defaultTitle,'Value',currentMenu.plot.defaultTitle);
          set(handlesMenu.plot_demedian,'Value',currentMenu.plot.demedian);
          set(handlesMenu.plot_title,'string',currentMenu.plot.title);
          set(handlesMenu.plot_labelIncrement,'string',num2str(currentMenu.plot.labelIncrement));
          set(handlesMenu.plot_labelOption,'string',currentMenu.plot.labelOption);
          doneWhat.plot = false;
          button_plot

        % Revese x axes  
        case abs('r')
          temp = currentMenu.plot.xlim0;
          currentMenu.plot.xlim0 = currentMenu.plot.xlim1;
          currentMenu.plot.xlim1 = temp;
          set(handlesMenu.plot_xlim0,'string',num2str(currentMenu.plot.xlim0));
          set(handlesMenu.plot_xlim1,'string',num2str(currentMenu.plot.xlim1));
          upToDate.plot = false;
          button_plot;

        % Change scale value
        case abs('s')
          disp('Enter a new value of the X scaling factor and hit return')
          [~,~,string] = ginput;
          tempNum = currentMenu.plot.xScale;
          try
            eval(['currentMenu.plot.xScale = ' char(string') ';']);
            set(handlesMenu.plot_xScale,'string',num2str(currentMenu.plot.xScale));
            upToDate.plot = false;
            button_plot
          catch
            warning(['X scaling factor "' char(string') '" is invalid'])
            eval(['currentMenu.plot.xScale = ' num2str(tempNum) ';']);
            set(handlesMenu.plot_xScale,'string',num2str(currentMenu.plot.xScale));
          end

        % Change scale option parameter
        case abs('t')
          disp('Enter a new value of the X scaling option (type) and hit return')
          [~,~,string] = ginput;
          tempString = currentMenu.plot.xScaleOption;
          try
            currentMenu.plot.xScaleOption =  char(string');
            set(handlesMenu.plot_xScaleOption,'string',currentMenu.plot.xScaleOption);
            upToDate.plot = false;
            button_plot
          catch
            warning(['X scaling Option "' char(string') '" is invalid'])
            currentMenu.plot.xScaleOption =  tempString;
            set(handlesMenu.plot_xScaleOption,'string',currentMenu.plot.xScaleOption);
          end

        % Update User
        case abs('u')
          disp('Enter a new user and hit return')
          [~,~,string] = ginput;
          currentMenu.active.user = char(string');
          set(handlesMenu.active_user,'string',currentMenu.active.user);
          upToDate.plot = false;
          button_plot

         % Set reduction velocity
        case abs('v')
          disp('Enter a new value of the reduction velocity and hit return')
          [~,~,string] = ginput;
          tempNum = currentMenu.plot.redVel;
          try
            eval(['currentMenu.plot.redVel = ' char(string') ';']);
            set(handlesMenu.plot_redVel,'string',num2str(currentMenu.plot.redVel));
            upToDate.plot = false;
            button_plot
          catch
            warning(['Recuction velocity parameter "' char(string') '" is invalid'])
            eval(['currentMenu.plot.redVel = ' char(string') ';']);
            set(handlesMenu.plot_redVel,'string',num2str(currentMenu.plot.redVel));
          end

        % Wiggle Plot Option
        case abs('w')
          disp('Enter a new value of the wiggle option and hit return')
          [~,~,string] = ginput;
          tempNum = currentMenu.plot.wiggleOption;
          try
            eval(['currentMenu.plot.wiggleOption = ' char(string') ';']);
            set(handlesMenu.plot_wiggleOption,'string',num2str(currentMenu.plot.wiggleOption));
            upToDate.plot = false;
            button_plot
          catch
            warning(['Wiggle Option parameter "' char(string') '" is invalid'])
            eval(['currentMenu.plot.wiggleOption = ' num2str(tempNum) ';']);
            set(handlesMenu.plot_wiggleOption,'string',num2str(currentMenu.plot.wiggleOption));
          end

        % Quit picking and keep menu parameters
        case abs('x')
          moreToDo = false;

        % Zoom in x 1/2     
        case abs('z')
          temp = xlim;
          temp = xm + [-0.25 0.25] * diff(temp);
          currentMenu.plot.xlim0 = temp(1);
          currentMenu.plot.xlim1 = temp(2);
          temp = ylim;
          temp = ym + [-0.25 0.25] * diff(temp);
          currentMenu.plot.tlim0 = temp(1);
          currentMenu.plot.tlim1 = temp(2);
          set(handlesMenu.plot_xlim0,'string',num2str(currentMenu.plot.xlim0));
          set(handlesMenu.plot_xlim1,'string',num2str(currentMenu.plot.xlim1));
          set(handlesMenu.plot_tlim0,'string',num2str(currentMenu.plot.tlim0));
          set(handlesMenu.plot_tlim1,'string',num2str(currentMenu.plot.tlim1));
          upToDate.plot = false;
          button_plot;

        % Increase pick error
        case abs('-')
         if active.picked(i)
            active.unc(i) = sd_round(active.unc(i)/sqrt(2),2);
            delete(handlesActive.unc1(i))
            delete(handlesActive.unc2(i))
            handlesActive.unc1(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ymOld(i) ymOld(i)]-active.unc(i),'r:');
            handlesActive.unc2(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ymOld(i) ymOld(i)]+active.unc(i),'r:');
            active.updated(i) = true;
         else
            warning('Cannot modify a pick error for an unpicked trace')
         end

        % Decrease pick error    
        case abs('+')
         if active.picked(i)
            active.unc(i) = sd_round(active.unc(i)*sqrt(2),2);
            delete(handlesActive.unc1(i))
            delete(handlesActive.unc2(i))
            handlesActive.unc1(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ymOld(i) ymOld(i)]-active.unc(i),'r:');
            handlesActive.unc2(i) = plot(x0Trace(i)+[-pickWidth/2 pickWidth/2],[ymOld(i) ymOld(i)]+active.unc(i),'r:');
            active.updated(i) = true;
         else
            warning('Cannot modify a pick error for an unpicked trace')
         end
         
        % Change "use" flag
        case abs('/')
          if active.picked(i)
            if ~active.use(i)
              active.use(i) = true;
              set(handlesActive.time(i),'color','r')
              set(handlesActive.unc1(i),'color','r')
              set(handlesActive.unc2(i),'color','r')
              active.updated(i) = true;
            else
              active.use(i) = false;
              set(handlesActive.time(i),'color',[1 0.6 0.6])
              set(handlesActive.unc1(i),'color',[1 0.6 0.6])
              set(handlesActive.unc2(i),'color',[1 0.6 0.6])
              active.updated(i) = true;
            end
          else
            warning('Cannot change "use" flag status for an unpicked trace')
          end
         
        % Rubberband use = true     
        case abs('[')
          disp('Rubberband in X-limits to set "use" flag to true')
          [xlim0, xlim1] = rubberband(hFig);
          index = find(x0Trace>=xlim0 & x0Trace<=xlim1);
          for i = index
            if active.picked(i)
              active.use(i) = true;
              set(handlesActive.time(i),'color','r')
              set(handlesActive.unc1(i),'color','r')
              set(handlesActive.unc2(i),'color','r')
              active.updated(i) = true;
            end
          end

        % Rubberband use = false     
        case abs(']')
          disp('Rubberband in X-limits to set "use" flag to false')
          [xlim0, xlim1] = rubberband(hFig);
          index = find(x0Trace>=xlim0 & x0Trace<=xlim1);
          for i = index
            if active.picked(i)
              active.use(i) = false;
              set(handlesActive.time(i),'color',[1 0.6 0.6])
              set(handlesActive.unc1(i),'color',[1 0.6 0.6])
              set(handlesActive.unc2(i),'color',[1 0.6 0.6])
              active.updated(i) = true;
            end
          end

        % Move plot to the left a half screen
        case abs('1')
          temp = xlim;
          if ~strcmp(get(gca,'xdir'),'reverse')
            temp = mean(temp) + [-1 0] * diff(temp);
          else
            temp = mean(temp) + [0 1] * diff(temp);
          end
          currentMenu.plot.xlim0 = temp(1);
          currentMenu.plot.xlim1 = temp(2);
          set(handlesMenu.plot_xlim0,'string',num2str(currentMenu.plot.xlim0));
          set(handlesMenu.plot_xlim1,'string',num2str(currentMenu.plot.xlim1));
          upToDate.plot = false;
          button_plot;

        % Move plot to the right a half screen
        case abs('2')
          temp = xlim;
          if ~strcmp(get(gca,'xdir'),'reverse')
            temp = mean(temp) + [0 1] * diff(temp);
          else
            temp = mean(temp) + [-1 0] * diff(temp);
          end
          currentMenu.plot.xlim0 = temp(1);
          currentMenu.plot.xlim1 = temp(2);
          set(handlesMenu.plot_xlim0,'string',num2str(currentMenu.plot.xlim0));
          set(handlesMenu.plot_xlim1,'string',num2str(currentMenu.plot.xlim1));
          upToDate.plot = false;
          button_plot;


        % Move plot down a half screen
        case abs('9')
          temp = ylim;
          if ~strcmp(get(gca,'ydir'),'reverse')
            temp = mean(temp) + [-1 0] * diff(temp);
          else
            temp = mean(temp) + [0 1] * diff(temp);
          end
          currentMenu.plot.tlim0 = temp(1);
          currentMenu.plot.tlim1 = temp(2);
          set(handlesMenu.plot_tlim0,'string',num2str(currentMenu.plot.tlim0));
          set(handlesMenu.plot_tlim1,'string',num2str(currentMenu.plot.tlim1));
          upToDate.plot = false;
          button_plot;

        % Move plot up a half screen
        case abs('0')
          temp = ylim;
          if ~strcmp(get(gca,'ydir'),'reverse')
            temp = mean(temp) + [0 1] * diff(temp);
          else
            temp = mean(temp) + [-1 0] * diff(temp);
          end
          currentMenu.plot.tlim0 = temp(1);
          currentMenu.plot.tlim1 = temp(2);
          set(handlesMenu.plot_tlim0,'string',num2str(currentMenu.plot.tlim0));
          set(handlesMenu.plot_tlim1,'string',num2str(currentMenu.plot.tlim1));
          upToDate.plot = false;
          button_plot;

      end

      if ~upToDate.plot
        button_plot
      end

    end

    if any(active.updated)
      doneWhat.madePick = true;
    end
    
  end
  
end

% Turn on menu
menu_on
drawnow

      
