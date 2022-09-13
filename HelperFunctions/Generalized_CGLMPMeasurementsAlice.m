function Max = Generalized_CGLMPMeasurementsAlice(d,m)

% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% The function Generalized_CGLMPMeasurementsAlice generates 
% generalized CGLMP Measurements for Alice according to the paper
% https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.97.170409


% Inputs: 
% 1) Dimension d
% 2) Number of measurements m

%Outputs:
% 1) Max = M_a|x, the gerneralized CGLMP measurements for Alice

% Requires: nothing

% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruﬂ
% Code by: Lucas Tendick
% Last Update: 01.09.2022



Id = eye(d);

% Set up some basis vectors

for l = 1:d
basis_vec(:,l) =  Id(:,l);
end

Max = zeros(d,d,d,m);
Max_vec = zeros(d,d,m); % vector representation

for x = 1:m
    alp(x) = (x-1/2)/m;
 for a = 1:d    
    
     for k = 0:d-1
            Max_vec(:,a,x) = Max_vec(:,a,x)+exp(2*pi*1i/d*k*((a-1)-alp(x)))*basis_vec(:,k+1);
     end
     Max_vec(:,a,x) = (1/sqrt(d))*Max_vec(:,a,x);
     Max(:,:,a,x) = Max_vec(:,a,x)*Max_vec(:,a,x)';
 end   
end



end
