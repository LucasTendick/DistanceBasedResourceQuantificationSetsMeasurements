function [is_fLHS,F,siglam] = fLHSTripartite2Unt(sigma)
%FLHSTRIPARTITE2UNT determines if a triparite assemblage with 2 untrusted
%parties is is fully-LHS (comes from a fully separable state)
% This function has one required argument:
%  sigma: a 6-D array containing the tripartite assemblage produced for
%  Charlie by Alice and Bob. The first two dimensions contain the dC x dC
%  unnormalised quantum states that are prepared. The last four dimensions
%  contain (a,b,x,y), the labels corresponding to the inputs and outcomes
%  and Alice and Bob, in the standard nonlocality ordering.
%
% is_fLHS = fLHSTripartite2Unt(sigma) is the indicator function for whether
% sigma could have arisen from local measurements on a fully separable
% quantum state. It returns 1 if this is the case, and zero otherwise. 
%
% [is_fLHS,F] = fLHSTripartite2Unt(sigma) returns the steering
% functional F that certifies that the assemblage sigma does not have an
% LHS model in this instance. F is a 6-D array. The first two dimensions
% contain the dC x dC elements of the steering functional measured by
% Charlie. The last four dimensions contain (a,b,x,y). A negative value
% certifies that a given assemblage demonstrates multipartite steering.
%
% [is_fLHS,F,siglam] = fLHSTripartite2Unt(sigma) returns the LHS model
% siglam in the case that the assemblage sigma has such a model. siglam is
% a 4-D array. The first two dimensions contain the dC x dC unnormalised
% quantum state that comprise the model. The last two dimensions are
% (mu,nu), the hidden variable respectively for Alice, Bob. If the
% assemblage demonstrates steering, the empty array siglam = [] is
% returned.
%
% requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
% authors: Paul Skrzypczyk, Daniel Cavalcanti
% last updated: March 17, 2016


tol = 1e-8; 
% this is the tolerance used to decide if steering has been demonstrated or
% not. This figure appears safe, but it can be changed with caution.

[dC,~,oa,ob,ma,mb] = size(sigma);
% dC = dim. of Charlie, oa = # outcomes for Alice, ma = # inputs for Alice,
% ob = # outcomes for Alice, mb = # inputs for Alice

NdetA = oa^ma; % number of deterministic strategies for Alice
SPA = genSinglePartyArray(oa,ma); % generate array containing 
                                            %single party strategies 
                                            
NdetB = ob^mb; % number of deterministic strategies for Alice
SPB = genSinglePartyArray(ob,mb); % generate array containing 
                                            %single party strategies 
            
SPAB = permute(repmat(SPA,[1,1,1,dC,dC,NdetB,ob,mb]),[4,5,3,6,1,7,2,8]).*...
    permute(repmat(SPB,[1,1,1,dC,dC,NdetA,oa,ma]),[4,5,6,3,7,1,8,2]);
% SPAB = ones(dC,dC)D(a|x,mu)D(b|y,nu). 
                                            
cvx_begin sdp quiet
    
    variable siglam(dC,dC,NdetA,NdetB) hermitian semidefinite
    % siglam are the members of the LHS model
    variable slack(dC,dC,oa,ob,ma,mb) hermitian semidefinite
    % slack variable in order to include an equality constraint
    dual variable F
    
    minimise real(sum(reshape(repmat(eye(dC),[1,1,oa,ob,ma,mb]).*slack,1,[])))
    % minimise sum_abxy tr slack_ab|xy 
    
    subject to
    
    F : sigma + slack == squeeze(sum(sum(repmat(siglam,[1,1,1,1,oa,ob,ma,mb])...
        .*SPAB,3),4));
    % sig_ab|xy  + slack_ab|xy == \sum_mu,nu D(a|x,mu)D(b|y,nu) sig_mu,nu 
    
cvx_end

% CVX will return a positive value if sigma does not have an fLHS model. It
% will return approx 0 (<= tol) if it is fLHS. This converts {>tol, <= tol}
% to {0,1}
is_fLHS = 1-(cvx_optval > tol);

% if the assemblage is fLHS, there is no inequality to return, so set it
% equal to the empty array.
if is_fLHS == 1
    F = [];
else % otherwise there is no model to return
    siglam = [];
end

end