function [InfoDual,Cax,rho_x,L_x] = InformativenessDiamondNormDual(Max,p)


% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)


% The function InformativenessDiamondNormDual quantifies the 
% informativeness according to Eq. (17) respectively SDP Eq. (80) (or Eq.
% (F2).

% Inputs: 
% 1) a 4-D array, containing the POVM elements. The first two
% dimensions contain the POVM elements (the Matrix representation),
% while the remaining two dimensions are (a,x), such that Max(:,:,a,x) = M_a|x.
% The outcome 'a' is the 3rd and the setting 'x' the 4th dimension.

% 2) A probability vector denoting the probabilities for choosing the
% settings x with probabbility p(x). For example p = ones(1,m)/m is
% a uniform distribution for m settings

%Outputs:
% 1) InfoDual, the value of the informativeness of M_a|x
% 2) Cax, hyperplane coefficients (matrcies that certify the
% informativeness in the first two dimension, the outcome 'a' in the
% third dimension, the setting 'x' in the fourth dimension)
% 3) states rho_x that fix the scale of the informativeness inequality
% 4) classical bounds, L_x that together give together the total bound
% \sum_x L_x


% Requires: CVX (http://cvxr.com/cvx/), A solver like SDPT3, Sedumi or
% Mosek, etc.


% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruß
% Code by: Lucas Tendick
% Last Update: 01.09.2022

% Note: The number of outcomes is assumed to be equal for each setting.
% Function can be used within CVX (as a convex constraint)



[d,~,o,m] = size(Max); %Get the dimension d, the number of outcomes o and 
                         %the number of settings m


cvx_begin % start cvx for the SDP.

%cvx_precision best
% note that the precision can sometimes lead to 
% errors depeding on your solver and your local settings
% the standard precision is typically good enough

variable L_x(1,m)
% variables for the local bound

variable rho_x(d,d,m) hermitian semidefinite
% variables for the states that fix the scale

variable Cax(d,d,o,m) hermitian semidefinite
% variables for the hyperplane coefficients 

expression S(1) 
% Place holder for the objective function

expressions class_strategy(o,m)
% Place holder for the classical strategies

% Constraints

for x = 1:m
    trace(rho_x(:,:,x)) == 1;
end

for x = 1:m
    for a = 1:o
        S = S+p(x)*trace(Max(:,:,a,x)*Cax(:,:,a,x));
        rho_x(:,:,x)-Cax(:,:,a,x) == hermitian_semidefinite(d);
    end
    S = S-L_x(1,x);
end

   for x = 1:m
        for a = 1:o
            class_strategy(a,x) = p(x)*trace(Cax(:,:,a,x));
            L_x(1,x)-class_strategy(a,x) >= 0;
        end
  end
    

Obj = real(S);

maximize real(Obj)

cvx_end

InfoDual = real(Obj);
end