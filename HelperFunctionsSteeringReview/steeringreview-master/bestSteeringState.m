function [beta, rho] = bestSteeringState(Max,F)
%BESTSTEERINGSTATE finds best state for steering functional and set of
%measurements of Alice
% This function has two required arguments:
%  Max: a 4-D array containing a valid set of POVMs. The first two
%  dimensions contain the dA x dA POVM elements, and the last two
%  dimensions are (a,x)
%  F: a 4-D array containing a steering functional. The first two
%  dimensions contain the dB x dB hermitian members of the functional, and
%  the last two dimensions are (a,x).
%
% beta = bestSteeringState(Max,F) returns the optimal violation of
% the steering functional F, with for the set of POVMs Max 
%
% [beta, rho] = bestSteeringState(Max,F) also returns the optimal state rho
% that achieve the optimal violation. rho is a 2-D array, containing the
% bipartite quantum state of dimension dA*dB x dA*dB.
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti
%   last updated: March 17, 2016


[dB,~,~,~] = size(F);
% dB = dim. of Bob
[dA,~,oa,ma] = size(Max);
% dA = dim. of Alice, oa = # outcomes for Alice, ma = # inputs for Alice

cvx_begin sdp quiet

    variable rho(dA*dB,dA*dB) hermitian semidefinite
    % the state
    
    maximise real(sum(reshape(genAssemblage(rho,Max).*conj(F),1,[])))
    % NOTE: we assume that the steering inequality is of the form
    % F.sigma^LHS <= beta, hence the maximisation.
    
    subject to
    
    trace(rho) == 1;
    % rho should be normalised
    
 cvx_end
 
 beta = cvx_optval;
 
end