function [Info,Fax] = InformativenessDiamondNorm(Max,p)


% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)


% The function InformativenessDiamondNorm quantifies the 
% informativeness according to Eq. (17) respectively SDP Eq. (79) (or Eq.
% (F1).

% Inputs: 
% 1) a 4-D array, containing the POVM elements. The first two
% dimensions contain the POVM elements (the Matrix representation),
% while the remaining two dimensions are (a,x), such that Max(:,:,a,x) = M_a|x.
% The outcome 'a' is the 3rd and the setting 'x' the 4th dimension.

% 2) A probability vector denoting the probabilities for choosing the
% settings x with probabbility p(x). For example p = ones(1,m)/m is
% a uniform distribution for m settings

%Outputs:
% 1) Info, the value of the informativeness of M_a|x
% 2) Fax, the closest uninformative measureable assemblage F_a|x


% Requires: CVX (http://cvxr.com/cvx/), A solver like SDPT3, Sedumi or
% Mosek, etc., QETLAB (http://www.qetlab.com) for the functions (PermuteSystems, DiamondNorm)


% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bru?
% Code by: Lucas Tendick
% Last Update: 01.09.2022

% Note: The number of outcomes is assumed to be equal for each setting.
% Function can be used within CVX (as a convex constraint)


[d,~,o,m] = size(Max); %Get the dimension d, the number of outcomes o and 
                       %the number of settings m

Id = eye(o); % Define an o-dimensional identity to use it for the register states |a><a|

for n = 1:o
    ket_a(:,n) = Id(:,n); % Get kets of the comp. bases
end

C = zeros(o*d,o*d,m); % Choi states for all m settings

for x = 1:m
for a = 1:o
            C(:,:,x) = C(:,:,x) + Tensor(ket_a(:,a)*ket_a(:,a)',transpose(Max(:,:,a,x)));
end
C(:,:,x) = PermuteSystems(C(:,:,x),[2,1],[o,d]);
end

%Qetlab uses Choi states on second subsystem, therefore we need to
%permute the Systems 1 and 2

cvx_begin % start cvx for the SDP.

%cvx_precision best
% note that the precision can sometimes lead to 
% errors depeding on your solver and your local settings
% the standard precision is typically good enough


variable Fax(d,d,o,m) hermitian semidefinite
% Define the variables for the uninformative assemblage

expressions Phi(o*d,o*d,m) hermitian semidefinite
% Place holder for the Choi states of the uninformative measurements

variable pax(o,m) 
% probability distributions

validPOVMs(Fax) == 1;

% constraint Fax to be valid (fixes the normalization of pax here)

Obj = 0; % initial value of the objective function

% Constraints and setting up the Choi states

for x = 1:m
for a = 1:o
        Fax(:,:,a,x) == pax(a,x)*eye(d);
        Phi(:,:,x) = Phi(:,:,x) + Tensor(ket_a(:,a)*ket_a(:,a)',transpose(Fax(:,:,a,x)));
end
Phi(:,:,x) = PermuteSystems(Phi(:,:,x),[2,1],[o,d]);
Obj = Obj + 1/2*p(x)*DiamondNorm(C(:,:,x)-Phi(:,:,x),[d,o]);
end

minimize real(Obj);

cvx_end

Info = real(Obj);

end