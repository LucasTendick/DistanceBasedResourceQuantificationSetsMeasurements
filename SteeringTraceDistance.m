function [STD, tau] = SteeringTraceDistance(sigma,p)

% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)
% The quantifier itself is from the paper: Phys. Rev. A 97, 022338 (2018)

% The function SteeringTraceDistance quantifies the 
% steerability according to Eq. (24) respectively SDP Eq. (60) (or Eq.
% (C6).

% Inputs: 
% 1) a 4-D array, containing the steering conditional states. The first two
% dimensions contain the assemblage elements (the Matrix representation),
% while the remaining two dimensions are (a,x), such that sigma(:,:,a,x) = \sigma_a|x.
% The outcome 'a' is the 3rd and the setting 'x' the 4th dimension.

% 2) A probability vector denoting the probabilities for choosing the
% settings x with probabbility p(x). For example p = ones(1,m)/m is
% a uniform distribution for m settings

%Outputs:
% 1) STD, the value of the steerability of sigma
% 2) tau, the closest LHS assemblage 


% Requires: CVX (http://cvxr.com/cvx/), A solver like SDPT3, Sedumi or
% Mosek, etc., Steeringreview-master (https://arxiv.org/abs/1604.00501) (https://github.com/paulskrzypczyk/steeringreview) 
% for the function (genSinglePartyArray). 

% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruﬂ
% Code by: Lucas Tendick
% Last Update: 01.09.2022

% Note: The number of outcomes is assumed to be equal for each setting.
% Function can be used within CVX (as a convex constraint)

[d,~,o,m] = size(sigma); %Get the dimension d, the number of outcomes o and 
                         %the number of settings m

Ndet = o^m; % caluclate the number of deterministic strategies
SingleParty = genSinglePartyArray(o,m); 
% get all the deterministic assignements

cvx_begin % start cvx for the SDP.

cvx_precision best
% note that the precision can sometimes lead to 
% errors depeding on your solver and your local settings
% the standard precision is typically good enough

variable siglam(d,d,Ndet) hermitian semidefinite
% Define variables for the local hiddens states

variable tau(d,d,o,m) hermitian semidefinite
% Define variables for the closest LHS assemblage

% Constraints

tau == squeeze(sum(repmat(siglam,[1,1,1,o,m])...
        .*permute(repmat(SingleParty,[1,1,1,d,d]),[4,5,3,1,2]),3));
    
   % tau_a|x = sum_lambda D(a|x,lambda) \sigmlam_lambda
  % See als JMPOVMs. 
  
    
trace(sum(siglam,3)) == 1;  % Normalization

% Turn the following constraint on if you want 
% a 'consistent' quantifier
%for x = 1:m
%   sum(tau(:,:,:,x),3) == sum(sigma(:,:,:,x),3);
%end

Obj = 0; % Set up the objective function
  

for x = 1:m
   for a = 1:o
       Obj = Obj + 1/2*p(x)*TraceNorm(sigma(:,:,a,x)-tau(:,:,a,x));  %create objective function
   end
end

minimize real(Obj)

STD = real(Obj);

cvx_end


end

