function [STDDual,Gax,L] = SteeringTraceDistanceDual(sigma,p)

% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)
% The quantifier itself is from the paper: Phys. Rev. A 97, 022338 (2018)

% The function SteeringTraceDistanceDual quantifies the 
% steerability according to Eq. (24) respectively SDP Eq. (66) (or Eq.
% (C12).

% Inputs: 
% 1) A 4-D array, containing the steering conditional states. The first two
% dimensions contain the assemblage elements (the Matrix representation),
% while the remaining two dimensions are (a,x), such that sigma(:,:,a,x) = sigma_a|x.
% The outcome 'a' is the 3rd and the setting 'x' the 4th dimension.

% 2) a probability vector denoting the probabilities for choosing the
% settings x with probabbility p(x). For example p = ones(1,m)/m is
% a uniform distribution for m settings

%Outputs:
% 1) STD_Dual, the value of the steerability of sigma
% 2) Gax, steering inequality 'coefficients' (4 Dimensional array, first two
% dimensions are the matrices, third dimension is the setting 'a', fourth dimension is the setting 'x')
% 3) The classical/local bound L


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

%cvx_precision best
% note that the precision can sometimes lead to 
% errors depeding on your solver and your local settings
% the standard precision is typically good enough

variable Gax(d,d,o,m) hermitian semidefinite
% Define variables for the steering inequality 

variable L(1)
% Define variable for the classical/local bound of the steering
% inequality

expressions local_strategy(d,d,Ndet)
% Place holder for the deterministic strategies/the constraints for 
% all lambda

Obj = 0; % set the initial value of the objective

% Set up the objective and the norm constraint for the
% Gax

% Constraints

for x = 1:m
    for a = 1:o
        Obj = Obj+p(x)*trace(Gax(:,:,a,x)*sigma(:,:,a,x));
        norm(Gax(:,:,a,x)) <= 1;
        % norm is the spectral norm (largest singular value)
    end
end

Obj = real(Obj-L); % final form of the Objective

% Set up the constraint that LHS assemblages can't violate
% the steering inequality

for lam = 1:Ndet
    for x = 1:m
        for a = 1:o
            local_strategy(:,:,lam) = local_strategy(:,:,lam)+p(x)*SingleParty(a,x,lam)*Gax(:,:,a,x);
        end
    end
    L*eye(d)-local_strategy(:,:,lam) == hermitian_semidefinite(d);
end

maximize Obj

cvx_end

STDDual = realObj;

end

