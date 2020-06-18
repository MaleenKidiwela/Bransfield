function string=cell2string(c)
% Creates a spaced delimited string from a character cell structure

if isempty(c)
  string = '';
  
else
  string = cell2mat(c(1));
  
  for i = 2:length(c)
    string = [string ' ' cell2mat(c(i))];
  end
end
