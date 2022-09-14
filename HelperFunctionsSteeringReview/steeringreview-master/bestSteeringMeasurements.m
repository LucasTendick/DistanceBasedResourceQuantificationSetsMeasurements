function [beta, Max] = bestSteeringMeasurements(rho,F)
%BESTSTEERINGMEASUREMENTS finds best measurements for steering functional
%and state
% This function has two required arguments:
%  rho: a 2-D array containing a valid bipartite quantum state.
%  F: a 4-D array containing a steering functional. The first two
%  dimensions contain the dB x dB hermitian members of the functional, and
%  the last two dimensions are (a,x).
%
% beta = bestSteeringMeasurements(rho,F) returns the optimal violation of
% the steering functional F, with for the quantum state rho. 
%
% [beta, Max] = bestSteeringMeasurements(rho,F) also returns the optimal
% measurements Max that achieve the optimal violation. Max is a 4-D array,
% containing the oa dA x dA POVM elements of the ma measurements. The first
% two dimensions contain the POVM elements, and the last two dimensions are
% (a,x).
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti
%   last updated: March 17, 2016


[dB,~,oa,ma] = size(F);
% dB = dim. of Bob, oa = # outcomes for Alice, ma = # inputs for Alice
dA = length(rho)/dB;
% dA = dim. of Alice.

sigma = zeros(dA,dA,oa,ma);
% sigma stores the 'assemblage' that Bob would create for Alice by applying
% the operators F_ax on rho. NOTE: This is not an assemblage, since the
% F_ax do not need to be positive operators, and do not need to sum to the
% identity. 

for a = 1:oa
    for x = 1:ma
        sigma(:,:,a,x) = PartialTrace(Tensor(eye(dA),F(:,:,a,x))*rho,2,[dA,dB]);
    end
end

cvx_begin sdp quiet

    variable Max(dA,dA,oa,ma) hermitian semidefinite
    % the POVM elements of Alice
    
    maximise real(sum(reshape(conj(sigma).*Max,1,[])))
    % NOTE: we assume that the steering inequality is of the form
    % F.sigma^LHS <= beta, hence the maximisation.
    
    subject to
    
    validPOVMs(Max) == 1;
    % Max should be a valid set of POVMs
    
 cvx_end
 
 beta = cvx_optval;
 
end