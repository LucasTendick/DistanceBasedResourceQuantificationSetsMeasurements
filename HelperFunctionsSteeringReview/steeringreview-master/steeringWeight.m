function [SW, F] = steeringWeight(sigma,varargin)
%STEERINGWEIGHT calculates the steering weight of an assemblage
% This function has one required argument:
%  sigma: a 4-D array, containing the members of the assemblage. The first 
%   two dimensions contain the (unnormalised) quantum states, while the
%   remaining two dimensions are (a,x), such that sigma(:,:,a,x) =
%   \sigma_a|x. 
%
% SW = steeringWeight(sigma) returns the steering weight SW of the
% assemblage sigma.
%
% [SW, F] = steeringWeight(sigma) also returns the steering functional F
% that certifies that the steering weight is SW. F is a 4-D array, with the
% first two dimensions containing the members of the steering functional,
% and the last two labelling (a,x).
%
% This function has one optional argument:
%   consistent: (default 0)
%
% [SW, F] = steeringWeight(sigma,consistent) calculates the standard
% steering weight when consistent = 0, and the consistent steering weight
% when consistent = 1.
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti
%   last updated: March 17, 2016

[consistent] = opt_args({0},varargin{:});
% if unspecified, it is assumed the standard steering weight is required.

[dB,~,oa,ma] = size(sigma);
% dB = dim. of Bob, oa = # outcomes for Alice, ma = # inputs for Alice

Ndet = oa^ma;
% number of determinstic probability distributions for Alice

SingleParty = genSinglePartyArray(oa,ma);
% generate array containing the single party distributions

% check that the assemblage is valid
if  NSAssemblage(sigma) == 0
    error('assemblage is not valid')
end

sigR = sum(sigma(:,:,:,1),3);
% sigR = reduced state of Bob, necessary for consistent steering robust.

% NOTE: Here we use the dual formulation of the steering weight.

cvx_begin sdp quiet
    
    variable F(dB,dB,oa,ma) hermitian semidefinite
    % members of the steering functional

    if consistent == 1
        variable X(dB,dB) hermitian
        % if consistent steering weight is required, we need an additional
        % variable.
    else
        X = zeros(dB,dB);
        % if we want the standard steering weight, we can set this variable
        % equal to the zero matrix.
    end    
    
    maximise 1 - real(sum(reshape(F.*conj(sigma),1,[])))
    % max 1 - sum_ax trace(F_ax*sigma_a|x)
    
    subject to
    
    for i = 1:Ndet
        sum(sum(permute(repmat(SingleParty(:,:,i),[1,1,dB,dB]),[3,4,1,2]).*F,3),4) ...
            - X + trace(X*sigR)*eye(dB) - eye(dB) == hermitian_semidefinite(dB);
        % standard SW: sum_ax F_ax D(a|x,lam) - eye(dB) >= 0, forall lam
        % consist. SW: sum_ax F_ax D(a|x,lam) - X + tr(X*sigR)eye(dB) -
        % eye(dB) >= 0, forall lam
    end

cvx_end

SW = cvx_optval;
    
end
