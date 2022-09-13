function [Dist] = AssemblageDiamondDistance(Max,Nax,p)

% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)

% The function AssemblageDiamondDistance computes
% the diamond distance according to Eq. (11) 

% Inputs: 
% 1) a 4-D array, containing the POVM elements of the assemblage M. The first two
% dimensions contain the POVM elements (the Matrix representation),
% while the remaining two dimensions are (a,x), such that Max(:,:,a,x) = M_a|x.
% The outcome 'a' is the 3rd and the setting 'x' the 4th dimension.

% 2) a 4-D array, containing the POVM elements of the assemblage N. The first two
% dimensions contain the POVM elements (the Matrix representation),
% while the remaining two dimensions are (a,x), such that Nax(:,:,a,x) = N_a|x.
% The outcome 'a' is the 3rd and the setting 'x' the 4th dimension.


% 3) a probability vector denoting the probabilities for choosing the
% settings x with probabbility p(x). For example p = ones(1,m)/m is
% a uniform distribution for m settings

%Outputs
% 1) Dist, the distance between the assemblages M and N



% Requires: CVX (http://cvxr.com/cvx/), A solver like SDPT3, Sedumi or
% Mosek, etc.


% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruß
% Code by: Lucas Tendick
% Last Update: 12.09.2022

% Note: The number of outcomes is assumed to be equal for each setting.
% Function can be used within CVX (as a convex constraint)

[d,~,o,m] = size(Max); %Get the dimension d, the number of outcomes o and 
                       %the number of settings m

Id = eye(o); % Define an o-dimensional identity to use it for the register states |a><a|

for n = 1:o
    ket_a(:,n) = Id(:,n); % Get kets of the comp. bases
end

% Set up the Choi states C for all m settings

C = zeros(o*d,o*d,m,2); 

for x = 1:m
    
for a = 1:o
    C(:,:,x,1) = C(:,:,x,1) + Tensor(ket_a(:,a)*ket_a(:,a)',transpose(Max(:,:,a,x)));
    C(:,:,x,2) = C(:,:,x,2) + Tensor(ket_a(:,a)*ket_a(:,a)',transpose(Nax(:,:,a,x)));
end


end

% calculate the distance

Dist = 0;

cvx_begin  % start cvx for the SDP.

cvx_precision best
% note that the precision can sometimes lead to 
% errors depeding on your solver and your local settings
% the standard precision is typically good enough

variable Zx(o*d,o*d,m) hermitian semidefinite

% Additional variable needed to define the diamond norm

% Constraints

% Set up the Choi states for F_a|x and the objective function

for x = 1:m

Zx(:,:,x)-(C(:,:,x,1)-C(:,:,x,2)) == hermitian_semidefinite(o*d); 

Dist = Dist + p(x)*norm(PartialTrace(Zx(:,:,x),1,[o,d]));
% norm is the spectral norm (largest singular value)
end

minimize Dist

cvx_end





end

