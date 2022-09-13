function  DoubleParty = DoublePartyStrategy(X,Y)

% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% The function DoublePartyStrategy generates the local deterministic
% strategies (probability distributions) that Alice and Bob have
% in a Bell experiment.

% Inputs: 
% 1)a vector X, whose length determines the number of measurement settings
% each entry of the vector determines how many outcomes are given in that
% setting. For instance X = [2,2,2] denotes that there are 3 measurements
% with 2 outcomes each. The same holds for Y (Bob). Note that in this
% version, the number of outcomes in each settings has to be equivalent

% Outputs:
% 1) DoubleParty, a 5 dim array containing the local deterministic
% strategies

% Requires: (https://github.com/paulskrzypczyk/steeringreview) 
% for the function (SinglePartyArray). 

% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruﬂ
% Code by: Lucas Tendick
% Last Update: 01.09.2022


NdetA = prod(X); % number of strategies for Alice
NdetB = prod(Y); % number of strategies for Bob

Ndet = NdetA*NdetB; % number of total strategies

% Generate the local strategies for Alice, respectively Bob

SinglePartyA = genSinglePartyArray(X(1,1),length(X));
SinglePartyB = genSinglePartyArray(Y(1,1),length(Y));


% Initialise their combined strategies

DoubleParty = zeros(X(1,1),Y(1,1),length(X),length(Y),Ndet);

cnt = 0;

for lamA = 1:NdetA
    for lamB = 1:NdetB
        cnt = cnt+1;
        for x = 1:length(X)
            for a = 1:X(1,1)
                for y = 1:length(Y)
                    for b = 1:Y(1,1)
                        DoubleParty(a,b,x,y,cnt) = SinglePartyA(a,x,lamA)*SinglePartyB(b,y,lamB);
                    end
                end
            end
        end
    end
end


end