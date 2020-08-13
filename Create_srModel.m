%creating srModel
load('1DVeloProf.mat')

%we first create a cordinate system
srModel.ghead(1) = -2;
srModel.ghead(2) = -0.5;
srModel.ghead(3) = 700; %x ranges from -8 to 6 and will have 50 nodes per km
srModel.ghead(4) = 450; %y ranges from -4 to 5 and will have 50 nodes per km
srModel.ghead(5) = 100; %z ranges from 0-10 km with 100 nodes.
srModel.ghead(6) = 0.02;
srModel.ghead(7) = 0.02;
srModel.ghead(8) = 0.1;
VProf=VProf(2:101,2);
srModel.P.u = ones(700,450,100);
srModel.P.anis_phi = [];
srModel.P.anis_fraction = [];
srModel.P.anis_theta = [];
srModel.S =[];
srModel.anis_sym = [];

for i = 1:srModel.ghead(3)
    for j = 1:srModel.ghead(4)
        
        srModel.P.u(i,j,:) = [VProf];
        
    end
end

