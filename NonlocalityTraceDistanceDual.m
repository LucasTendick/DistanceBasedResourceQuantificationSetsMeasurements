function [NLDual,C_abxy,L] = NonlocalityTraceDistanceDual(p,X,Y,mu)

% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% The quantifier itself is from the paper: Phys. Rev. A 97, 022111 (2018)

% The function NonlocalityTraceDistanceDual quantifies the 
% nonlocality according to Eq. (26) respectively SDP Eq. (56) (or Eq.
% (C3).

% Inputs: 
% 1) A propability distribution p = p(a,b|x,y) which is a 4 dimensional
% array containing the outcome probabilities of Alice (a) and Bob (b)
% given that they chose setting x (Alice), respectively y (Bob). 

%2) Vectors X,Y that have length of the number of measurements
% ma (Alice), respectively mb (Bob). Each entry of a vector describes
% the nuber of outcomes in that setting. X = [2,2,2] means that
% Alice has 3 measurements with 2 outcomes in each setting.
% Note: This version requires the same number of outcomes in each 
% setting for Alice and Bob

%3) A matrix mu that contains the input probabilities mu = mu(x,y)
% which incodes with which probability Alice and Bob chose 
% a given setting

%Outputs:
% 1) NLDual, the value of the incompatibility of p
% 2) C_abxy, the coefficients of the optimal Bell inequality
% 3) L, the local bound of the Bell inequality


% Requires: CVX (http://cvxr.com/cvx/), A solver like SDPT3, Sedumi or
% Mosek, etc., % DoublePartyStrategy

% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruﬂ
% Code by: Lucas Tendick
% Last Update: 01.09.2022

% Function can be used within CVX (as a convex constraint)

DoubleParty = DoublePartyStrategy(X,Y);

ma = size(X,2); % Number of Settings of Alice
mb = size(Y,2); % Number of Settings of Bob
oa = X(1,1);
ob = Y(1,1);

Ndet = size(DoubleParty,5); 
% Number of local deterministic strategies


cvx_begin 

variable C_abxy(oa,ob,ma,mb)
% variables for the Bell inequallity
expressions local_strategy(Ndet,1)
% place holder the the local_deterministic strategies
variable L(1)
% variable for the local bound

Obj = 0; % inital value of the objective function

for x = 1:ma
    for y = 1:mb
        for a = 1:oa
            for b = 1:ob
                 C_abxy(a,b,x,y) <= 1; 
                 C_abxy(a,b,x,y) >= 0;
                 % The Bell coefficients should be between 0 and 1
                Obj = Obj+mu(x,y)*C_abxy(a,b,x,y)*p(a,b,x,y);
            end
        end
    end
end

% Constraints for the local deterministic strategies

for lam = 1:Ndet
   for x = 1:ma
       for y = 1:mb
           for a = 1:oa
               for b = 1:ob
                   local_strategy(lam,1) = local_strategy(lam,1)+mu(x,y)*C_abxy(a,b,x,y)*DoubleParty(a,b,x,y,lam);
               end
           end
       end
   end
   local_strategy(lam,1) <= L;
end

Obj = Obj-L;

maximize real(Obj)

cvx_end

NLDual = real(Obj);

end

