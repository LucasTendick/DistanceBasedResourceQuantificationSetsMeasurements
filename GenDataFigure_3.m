% This file generates the Data for Figure 3 in the work
%work 'Distance-based resourc quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruﬂ
% Code by: Lucas Tendick
% Last Update: 01.09.2022

% Relies on the functions in
% 'DistanceBasedResourceQuantificationSetsMeasurements'
% Uses implicitly Helper Functions, HelperFunctionsQetlab,
% HelperFunctionsSteeringReview
% Uses explicitly the function 'MaxEntangled' from Qetlab
% to generate the maximally entangled state
% Uses explicitly the function 'genAssemblage' from 
% SteeringReview


m = 2;
for d = 2:9
rho = MaxEntangled(d)*MaxEntangled(d)';

X = d*ones(1,m);
Y = d*ones(1,m);

Max = Generalized_CGLMPMeasurementsAlice(d,m);
Mby = Generalized_CGLMPMeasurementsBob(d,m);
sigma = genAssemblage(rho,Max);

p = 1/m*ones(1,m);
mu = 1/m^2*ones(m,m);

for x = 1:length(X)
    for y = 1:length(Y)
        for a = 1:d
            for b = 1:d
                q(a,b,x,y) = trace(Tensor(Max(:,:,a,x),Mby(:,:,b,y))*rho);
            end
        end
    end
end

[Incomp(m,d),~] = IncompatibilityDiamondNorm(Max,p);
[STD(m,d), ~] = SteeringTraceDistance(sigma,p);
[NL(m,d),~] = NonlocalityTraceDistance(q,X,Y,mu);


end