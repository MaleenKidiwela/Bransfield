% Turn all the menu items on
% 1/9/20 Need set background color to something else and then back to current value to get it to appear 
temp = get(handlesMenu.figure1,'children');
for i = 1:length(temp)
  set(temp(i),'enable','on')
  tempNum = get(temp(i),'backgroundcolor');
  set(temp(i),'backgroundcolor','red');
  set(temp(i),'backgroundcolor',tempNum);
  drawnow
end
