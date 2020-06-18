function c=string2cell(string,delimChar);
% Breaks up a string into a cell structure based on space or other demarking characters

if nargin<2
  delimChar = ' ';
end

if isempty(string)
  c = {};
  
else
  igap = find(string ~= abs(delimChar));

  if isempty(igap) 
    c = {};
    
  else
    i0 = [1 find(diff(igap)>1)+1];
    i1 = [find(diff(igap)>1) length(igap)];
    for j = 1:length(i0)
      c(j) = {string(igap(i0(j)):igap(i1(j)))};
    end
  end
end
