% This file generates the Data for Figure 4 and 5 in the work
%work 'Distance-based resourc quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruﬂ
% Code by: Lucas Tendick
% Last Update: 12.09.2022

m = 2;
dim = [2,3,5];
for cnt = 1:size(dim,2)
clear Max 
clear Max_eta

d = dim(cnt);
Max = MUB_Measurements(d,m);
Max_eta(:,:,:,1) = Max(:,:,:,1);

eta_cnt = 0;
for eta = 0.1:0.1:1
    eta_cnt = eta_cnt+1;
    
   for a = 1:d
       Max_eta(:,:,a,2) = eta*Max(:,:,a,2)+(1-eta)*trace(Max(:,:,a,2))*eye(d)/d;
   end
   [IncompDualOpt(d,eta_cnt),~,~,~,p(:,d,eta_cnt)] = IncompDiamondNormDualOptDistribution(Max_eta);

end


end