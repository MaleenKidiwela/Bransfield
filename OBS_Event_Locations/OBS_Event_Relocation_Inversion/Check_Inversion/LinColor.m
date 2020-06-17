function [H,wech2] = LinColor(Vector, Vlims, HOTT)
% LinColor is a function which takes a Nx1 vector of data and
% outputs a matrix of Nx3 numbers corresponding to red,green, & blue color
% handles for a colorbar.
% INPUTS:  Vector - is a Nx1 vector of data by which you want color code
%          Vlims  - is a 2 element start/end point vector for the color scheme
%                   ([Vmin Vmax] is default)
%                   
%          HOTT - is a string corresponding to the type of colormap you
%                 want.
%                 choose HOTT color scheme: jet, hsv, hot, cool, spring, 
%                 summer, autumn, winter, gray, bone, copper, pink, or lines.
%                 Personally, I like 'spring'.
%
%
% OUTPUT:  H - is a Nx3 vector corresponding to the color (each column is
%              a number relating to the amount of red,green, or blue)
%
%
% To utilize: feed in vector to get H. Then loop over the vector and plot 
%             giving a new H(i,:) color each time.
% 
% Example....
%
% z=100*rand(100,1);
% H=LinColor(z,[min(z) max(z)],'jet');
% for k=1:length(z)
%   plot(z(k),1,'o','MarkerFaceColor',H(k,:),'MarkerEdgeColor',H(k,:));hold on;
% end
% caxis([min(z) max(z)]);
% colorbar
%
% Fixed input checking, W. Wilcock, March 5

if nargin<3
  HOTT='jet';
end

if nargin<3
  Vmin=min(Vector);  % the earliest day you want to start plotting
  Vmax=max(Vector);  % the day after the latest day you want to plot
else
  Vmin=Vlims(1);
  Vmax=Vlims(2);
end

% cmap=colormap;

wech=10000; % an arbitrary big number to keep a smooth colormap gradient
% set h (a 3 by wech matrix) to the chosen color scheme:
eval(sprintf('h=colormap(%s(wech));',HOTT)); % insert "HOTT" colorscheme
this=linspace(Vmin, Vmax,length(h)); % set up x variable for interpolation...
for i=1:length(Vector) % loop over locations to be plotted
  H(i,1)=interp1(this,h(:,1),Vector(i),'linear'); % interpolate to get red color from day
  H(i,2)=interp1(this,h(:,2),Vector(i),'linear'); % interpolate to get green color from day
  H(i,3)=interp1(this,h(:,3),Vector(i),'linear'); % interpolate to get blue color from day
end
wech2=colormap(HOTT);
% colormap(cmap);
return