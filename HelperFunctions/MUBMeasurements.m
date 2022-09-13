function Max = MUBMeasurements(p,m)

% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% The function MUBMeasurements generates a set of m MUB
% in a prime dimension p

% Inputs: 

% 1) number of settings m
% 2) prime dimension p

%Outputs:
% 1) Max = M_a|x, which is an assemblage of MUB measurements


% Requires: QETLAB (http://www.qetlab.com) for the function (GenPauli)


% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruﬂ
% Code by: Lucas Tendick
% Last Update: 01.09.2022




d = p;  % d has to be a prime number

% number of settings for measurements, can be in the 
% interval [1,d+1]
  
o = d; % number of outcomes is d (rank-1 projeactive measurements)         

% GenPauli gives U = X^j*Z^k in d dimensions with integers
% j and k which are in the interval [0,d-1]

X = GenPauli(1,0,d);
Z = GenPauli(0,1,d); 

[U(:,:,1),~] = eig(X);
[U(:,:,2),~] = eig(Z);

for n = 1:d-1
        [U(:,:,2+n),~] = eig(X*(Z^n));
end

Max = zeros(d,d,o,m);
for x = 1:m
    for a = 1:o
        Max(:,:,a,x) = U(:,a,x)*U(:,a,x)';
    end
end



end
