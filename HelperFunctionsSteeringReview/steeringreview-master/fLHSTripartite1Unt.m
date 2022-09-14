function [is_fLHS,F,siglam] = fLHSTripartite1Unt(sigma,k,varargin)
%FLHSTRIPARTITE1UNT determines if a triparite assemblage with 1 untrusted
%party is is fully-LHS (comes from a fully separable state)
% This function has two required arguments:
%  sigma: a 4-D array, containing the tripartite assemblage prepared for
%  Bob-Charlie by Alice. The first two dimensions contain the dB*dC x dB*dC
%  (unnormalised) quantum states of Bob-Charlie. The last two dimensions
%  contain (a,x), the outcome and input of Alice respectively.
%  k: An integer, which specifies what size symmetric PPT extension each
%  member of the assemblage should have. That is, if k = 2, then each
%  member of the LHS model will be required to have an (unnormalised)
%  2-symmetric PPT extension (as a relaxation of separability).
%  
% This function has one optional argument:
%  dims: a 2-D array containing the dimensions of Bob and Charlie. The
%  default is dB = dC = sqrt(dB*dC) (determined from sigma directly). This
%  argument only needs to be specified if Bob and Charlie do not hold
%  states of the same dimension.
%
% is_fLHS = fLHSTripartite1Unt(sigma,k) is the indicator function for
% whether sigma could have arisen from measurements on a fully-separable
% state or not. This is an outer-approximation, and uses the set of
% k-symmetric PPT extendible states as an outer approximation to the set of
% separable states. 
%
% is_fLHS = fLHSTripartite1Unt(sigma,k,dims) sets the local dimension of
% Bob to dims(1), and the local dimension of Charlie to dims(2). 
%
% [is_fLHS,F] = fLHSTripartite1Unt(sigma,k) returns in F the steering
% functional that certifies that the assemblage sigma could not have arisen
% from a fully separable state. F is a 4-D array. The first two dimensions
% contain the dB*dC x dB*dC operators for Bob-Charlie. The last two
% dimensions are (a,x). A negative value certifies that a given assemblage
% demonstrates multipartite steering. In the case that the assemblage has
% an LHS model, the empty array F = [] is returned.
%
% [is_fLHS,F,siglam] = fLHSTripartite1Unt(sigma,k) returns the LHS model
% siglam in the case that sigma has an fLHS model. siglam is a 3-D array.
% The first two dimensions contain the dB*dC x dB*dC unnormalised quantum
% states which comprise the model. The third dimension gives the value of
% the hidden variable. If the assemblage demonstrates steering, then the
% empty array is returned siglam = [].
%
% requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
% authors: Paul Skrzypczyk, Daniel Cavalcanti
% last updated: March 17, 2016

tol = 1e-8; 
% this is the tolerance used to decide if steering has been demonstrated or
% not. This figure appears safe, but it can be changed with caution.

[dBdC,~,oa,ma] = size(sigma);
% dBdC = dim. of Bob x dim. of Charlie, 
% oa = # outcomes for Alice, ma = # inputs for Alice

[dims] = opt_args({[round(sqrt(dBdC)),round(sqrt(dBdC))]},varargin{:}); 
%fix values of optional inputs: defaults: dB = dC = sqrt(dB*dC), nm = 0
%(corresponding to unnormalised)
dB = dims(1); dC = dims(2);
% dB = dim. of Bob, dC = dim. of Charlie.

Ndet = oa^ma; % number of deterministic strategies for Alice
SingleParty = genSinglePartyArray(oa,ma); % generate array containing 
                                            %single party strategies 

if isa(sigma,'cvx') == 0 % if sigma isn't a CVX variable  
    % check that the assemblage is valid
    if  NSAssemblage(sigma) == 0
        error('assemblage is not valid')
    end 
end

cvx_begin sdp quiet
    
    variable siglam(dB*dC,dB*dC,Ndet) hermitian semidefinite
    % siglam are the members of the LHS model
    variable slack(dB*dC,dB*dC,oa,ma) hermitian semidefinite
    % slack variable in order to include an equality constraint
    dual variable F
    
    minimise real(sum(reshape(repmat(eye(dB*dC),[1,1,oa,ma]).*slack,1,[])))
    % minimise sum_ax tr slack_a|x 
    
    subject to
    
    F : sigma + slack == squeeze(sum(repmat(siglam,[1,1,1,oa,ma])...
        .*permute(repmat(SingleParty,[1,1,1,dB*dC,dB*dC]),[4,5,3,1,2]),3));
    % sig_a|x  + slack_a|x == \sum_lam D(a|x,lam) sig_lam 

    for i = 1:Ndet
        SymmetricExtension(siglam(:,:,i),k,[dB,dC],1,1) == 1;
    end
    % siglam(:,:,i) should have a k-symmetric PPT extension.
    
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
