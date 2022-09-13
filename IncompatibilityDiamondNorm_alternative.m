function [Incomp, Fax] = IncompatibilityDiamondNorm_alternative(Max,p)

% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% The function IncompatibilityDiamondNorm_alternative quantifies the 
% incompatibility according to Eq. (21) respectively SDP Eq. (69) (or Eq.
% (E1) and is an alternative formulation of the function IncompatibilityDiamondNorm
% which does not require the DiamondNorm and the PermuteSystems functions from Qetlab, as it is
% implemented directly

% Inputs: 
% 1) a 4-D array, containing the POVM elements. The first two
% dimensions contain the POVM elements (the Matrix representation),
% while the remaining two dimensions are (a,x), such that Max(:,:,a,x) = M_a|x.
% The outcome 'a' is the 3rd and the setting 'x' the 4th dimension.

% 2) a probability vector denoting the probabilities for choosing the
% settings x with probabbility p(x). For example p = ones(1,m)/m is
% a uniform distribution for m settings

%Outputs:
% 1) Incomp, the value of the incompatibility of M_a|x
% 2) Fax, the closest jointly measureable assemblage F_a|x


% Requires: CVX (http://cvxr.com/cvx/), A solver like SDPT3, Sedumi or
% Mosek, etc., 
% Steeringreview-master (https://arxiv.org/abs/1604.00501) (https://github.com/paulskrzypczyk/steeringreview) 
% for the function (JMPOVMs). 

% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruß
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

% Set up the Choi states C for all m settings

C = zeros(o*d,o*d,m); 

for x = 1:m
    for a = 1:o
    C(:,:,x) = C(:,:,x) + Tensor(ket_a(:,a)*ket_a(:,a)',transpose(Max(:,:,a,x)));
    end
end


cvx_begin  % start cvx for the SDP.

cvx_precision best
% note that the precision can sometimes lead to 
% errors depeding on your solver and your local settings
% the standard precision is typically good enough



variable Fax(d,d,o,m) hermitian semidefinite
% Define a 4 dimensional array for the jointly measureable
% assemblage, which has the same size as Max.

variable Zx(o*d,o*d,m) hermitian semidefinite

% Additional variable needed to define the diamond norm

expressions Phi(o*d,o*d,m) hermitian semidefinite
% Define the variables for respective Choi states

% Constraints

JMPOVMs(Fax) == 1;
% constraint the assemblage F_a|x to be jointly measureable

% validPOVMs(Fax) == 1; 
% This constraint is redundant if 
% JMPOVMs is already an active constraint, you can still turn it on.


Obj = 0; % set the objective function to an initial value of 0

% Set up the Choi states for F_a|x and the objective function

for x = 1:m
for a = 1:o
    Phi(:,:,x) = Phi(:,:,x) + Tensor(ket_a(:,a)*ket_a(:,a)',transpose(Fax(:,:,a,x)));
end

Zx(:,:,x)-(C(:,:,x)-Phi(:,:,x)) == hermitian_semidefinite(o*d); 
% Z_x >=  C_x-Phi_x
Obj = Obj + p(x)*norm(PartialTrace(Zx(:,:,x),1,[o,d]));
% norm is the spectral norm (largest singular value)
end

minimize real(Obj)

cvx_end

Incomp = real(Obj);


end

