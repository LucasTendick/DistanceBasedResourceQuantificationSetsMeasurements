function [IncompDual,Cax,rho_x,L] = IncompatibilityDiamondNormDual(Max,p)

% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% The function IncompatibilityDiamondNorm quantifies the 
% incompatibility according to Eq. (21) respectively SDP Eq. (70) (or Eq.
% (E2).

% Inputs: 
% 1) a 4-D array, containing the POVM elements. The first two
% dimensions contain the POVM elements (the Matrix representation),
% while the remaining two dimensions are (a,x), such that Max(:,:,a,x) = M_a|x.
% The outcome 'a' is the 3rd and the setting 'x' the 4th dimension.

% 2) a probability vector denoting the probabilities for choosing the
% settings x with probabbility p(x). For example p = ones(1,m)/m is
% a uniform distribution for m settings

%Outputs:
% 1) IncompDual, the value of the incompatibility of M_a|x
% 2) Cax, 4 dimensional 'coefficients' (matricies) of the hyperplane 
% separating M_a|x for the set of jointly measureable assemblages
% 3) rho_x, states that constraint the SDP such that Incomp <= 1;
% 4) L, the bound achieved by jointly measureable assemblages

% Requires: CVX (http://cvxr.com/cvx/), A solver like SDPT3, Sedumi or
% Mosek, etc. Steeringreview-master (https://arxiv.org/abs/1604.00501) (https://github.com/paulskrzypczyk/steeringreview) 
% for the function (SinglePartyArray). 

% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruß
% Code by: Lucas Tendick
% Last Update: 01.09.2022




[d,~,o,m] = size(Max); %Get the dimension d, the number of outcomes o and 
                       %the number of settings m


SingleParty = genSinglePartyArray(o,m);
% get all the deterministic assignements

Ndet = size(SingleParty,3); % get the number of deterministic strategies


cvx_begin % start cvx for the SDP.

cvx_precision best
% note that the precision can sometimes lead to 
% errors depeding on your solver and your local settings
% the standard precision is typically good enough


variable L(d,d) hermitian
% Variable for the classical bound

variable rho_x(d,d,m) hermitian semidefinite
% Variables for the quantum states rho_x

variable Cax(d,d,o,m) hermitian semidefinite
% Variables for the hyperplane coefficients

expression S(1) % place holder for the Objective function

expressions det_strategy(d,d,Ndet) 
% place holder for the constraints for all lambda

% Constraints and setting up the objective

for x = 1:m
    trace(rho_x(:,:,x)) == 1; 
    for a = 1:o
        S = S+p(x)*trace(Max(:,:,a,x)*Cax(:,:,a,x));
        rho_x(:,:,x)-Cax(:,:,a,x) == hermitian_semidefinite(d);
    end
end

% Constraints for all deterministic strategies

for lam = 1:Ndet
    for x = 1:m
        for a = 1:o
            det_strategy(:,:,lam) = det_strategy(:,:,lam)+p(x)*SingleParty(a,x,lam)*Cax(:,:,a,x);
        end
    end
   L-det_strategy(:,:,lam) == hermitian_semidefinite(d);
end

Obj = real(S-trace(L)); % Define the final objective function

maximize real(Obj)

cvx_end

IncompDual = real(Obj);
end

