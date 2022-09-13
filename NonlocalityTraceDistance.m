function [NL,q_loc] = NonlocalityTraceDistance(p,X,Y,mu)

% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% The quantifier itself is from the paper: Phys. Rev. A 97, 022111 (2018)

% The function NonlocalityTraceDistance quantifies the 
% nonlocality according to Eq. (26) respectively SDP Eq. (59) (or Eq.
% (C5).

% Inputs: 
% 1) a propability distribution p = p(a,b|x,y) which is a 4 dimensional
% array containing the outcome probabilities of Alice (a) and Bob (b)
% given that they chose setting x, respectively y. 

%2) Vectors X,Y that have length of the number of measurements
% ma (Alice), respectively mb (Bob). Each entry of a vector describes
% the nuber of outcomes in that setting. X = [2,2,2] means that
% Alice has 3 measurements with 2 outcomes in each setting.
% Note: This versio requires the same number of outcomes in each 
% setting for Alice and Bob

%3) A matrix mu that contains the input probabilities mu = mu(x,y)
% which incodes with which probability Alice and Bob chose 
% a given setting

%Outputs:
% 1) NL, the value of the incompatibility of p
% 2) q_loc, the closest local probability distribution


% Requires: CVX (http://cvxr.com/cvx/), A solver like SDPT3, Sedumi or
% Mosek, etc., % DoublePartyStrategy

% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruß
% Code by: Lucas Tendick
% Last Update: 01.09.2022

% Function can be used within CVX (as a convex constraint)

% generate the local deterministic probability distributions/strategies
% of Alice and Bob
DoubleParty = DoublePartyStrategy(X,Y); 

ma = size(X,2); % Number of Settings of Alice
mb = size(Y,2); % Number of Settings of Bob

Ndet = size(DoubleParty,5); 
% Number of local deterministic strategies

cvx_begin % start cvx for the SDP.

%cvx_precision best
% note that the precision can sometimes lead to 
% errors depeding on your solver and your local settings
% the standard precision is typically good enough

variable q_lam(Ndet) 
% variable for the distribution of the local hidden
% variable

expression q_loc(X(1,1),Y(1,1),length(X),length(Y));
% place holder for the closest local probability distribution

% generate the closest local probability distribution

for x = 1:ma
    oa = X(1,x);
   for y = 1:mb
       ob = Y(1,y);
    for a = 1:oa
       for b = 1:ob
          for lam = 1: Ndet
              
              q_loc(a,b,x,y) = q_loc(a,b,x,y)+q_lam(lam)*DoubleParty(a,b,x,y,lam);
              
          end
       end
    end
   end
end

% Enable these constraints for a 'consistent' quantifier

%for x = 1:ma
%   sum(q_loc(:,:,x,:),1) ==  sum(p(:,:,x,:),1);
%end

%for y = 1:mb
%   sum(q_loc(:,:,:,y),2) ==  sum(p(:,:,:,y),2);
%end

% Constraints


for lam = 1:Ndet
   q_lam(lam) >= 0; 
end

sum(q_lam) == 1;

Obj = 0; % initial value of the objective function
  
% generate the final objective function

for x = 1:ma 
    oa = X(1,x);
   for a = 1:oa
       for y = 1:mb
           ob = Y(1,y);
           for b = 1:ob
                Obj = Obj + mu(x,y)*norm(q_loc(a,b,x,y)-p(a,b,x,y));  %create obj_function
           end
       end
   end
end


minimise Obj


cvx_end

NL = 1/2*Obj;



end