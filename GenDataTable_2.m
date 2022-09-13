% This file is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% It generates the MUB bounds from Table 2, i.e., the upper and lower
% bounds in Eq. (43). Furthermore the corresponding incompatibility
% is calculated numerically and the (better) lower bound from Eq. (42)
% is also calculated



% Requires: CVX (http://cvxr.com/cvx/), A solver like SDPT3, Sedumi or
% Mosek, etc. 

% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruﬂ
% Code by: Lucas Tendick
% Last Update: 01.09.2022

% Note: Denpending on your available RAM, it will not be possible
% to go to high dimensions. In particular here the problem could 
% most likely occure for d =7 for the settings m=7 and m=8. 
% In general the runtime of this program is relatively long if performed
% for all values in Table 2.
% Try also using a different solver to improve the runtime. 
% However in that case we do not even have to run the SDP since the lower bound
% in Eq. (42) is analytically proven to be tight (to coincide with the
% actual incompatibility in that case.
% We can however, just run the primal and the dual of the SDP
% 'by hand' to convince us, that the values coincide. 'By hand'
% means here, that we can run the primal and dual SDP separetly
% and enforce a particular feasible point, that we have shown 
% analytically to be feasible.
% We incooperate this here in a simple style by computing the distance 
% to Max_eta which are the depolarized measurements with eta = (d*T-m)/(d*m-m); 


dim = [2,3,5,7];

for cnt = 1:size(dim,2)
clear Max
clear Max_eta
d = dim(cnt);    
    
for m = 2:d+1

Max = MUBMeasurements(d,m);
[~,~,o,~] = size(Max);
p = ones(1,m)/m; % input distribution
    
    
    Incomp_up(m,d) = (d-1)*(m-1)/((d+1)*m); 
    Incomp_low(m,d) = 1-(1/m)*(1+(m-1)/sqrt(d));
    
    % Preparation for the lower bound from Eq. (42)
    SingleParty = genSinglePartyArray(d,m);
    Ndet = size(SingleParty,3);
    
    det_strategy = zeros(d,d,Ndet);
    E = zeros(1,Ndet);
    
for lam = 1:Ndet
    for x = 1:m
        for a = 1:o
            det_strategy(:,:,lam) = det_strategy(:,:,lam)+SingleParty(a,x,lam)*Max(:,:,a,x);
        end
    end
    E(lam) = lambda_max(det_strategy(:,:,lam));
    
end
T(m,d) = max(E);
Incomp_low_improved(m,d) = 1-T(m,d)/m;
[Incomp(m,d),~,~,~] = IncompatibilityDiamondNorm(Max,p);

eta = (d*T(m,d)-m)/(d*m-m);
for x = 1:m
    for a = 1:o
        Max_eta(:,:,a,x) = eta*Max(:,:,a,x)+(1-eta)*trace(Max(:,:,a,x))*eye(d)/d;
    end
end

[Incomp_eta(m,d)] = AssemblageDiamondDistance(Max,Max_eta,ones(1,m)/m);
end

end





