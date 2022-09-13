function [IncompDual,Cax,rho_x,L,p] = IncompDiamondNormDualOptDistribution(Max)


% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% The function IncompatibilityDiamondNormDualOptDirstribution is
% essentially the same as IncompatibilityDiamondNormDual (see also
% there for more details) but here, we optimize also over the 
% input probability distribution

% Inputs: 
% 1) a 4-D array, containing the POVM elements. The first two
% dimensions contain the POVM elements (the Matrix representation),
% while the remaining two dimensions are (a,x), such that Max(:,:,a,x) = M_a|x.
% The outcome 'a' is the 3rd and the setting 'x' the 4th dimension.



%Outputs:
% 1) IncompDual, the value of the incompatibility of M_a|x
% 2) Cax, 4 dimensional 'coefficients' (matricies in the first two dimensions, 
% the outcomes 'a' in the third dimension, the setting 'x' in the fourth dimension) 
% of the hyperplane separating M_a|x for the set of jointly measureable assemblages
% 3) rho_x, states that constraint the SDP such that Incomp <= 1;
% 4) L, the bound achieved by jointly measureable assemblages
% 5) p, the optimal input distribution

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

variable p(m)
% variable for the input probability distribution


% Constraints

for x = 1:m
    trace(rho_x(:,:,x)) == p(x);
end
sum(p) == 1;

for x = 1:m
    for a = 1:o
        S = S+trace(Max(:,:,a,x)*Cax(:,:,a,x));
        rho_x(:,:,x)-Cax(:,:,a,x) == hermitian_semidefinite(d);
    end
end



for lam = 1:Ndet
    for x = 1:m
        for a = 1:o
            det_strategy(:,:,lam) = det_strategy(:,:,lam)+SingleParty(a,x,lam)*Cax(:,:,a,x);
        end
    end
    L-det_strategy(:,:,lam) == hermitian_semidefinite(d);
end

Obj = real(S-trace(L));

maximize Obj

cvx_end

IncompDual = Obj;
end

