function [is_bLHS,F] = bLHSTripartite1Unt(sigma,k,varargin)
%BLHSTRIPARTITE1UNT determines if a triparite assemblage with 1 untrusted
%party is biseparable-LHS (comes from a biseparable state)
% This function has two required arguments:
%  sigma: a 4-D array, containing the tripartite assemblage prepared for
%  Bob-Charlie by Alice. The first two dimensions contain the dB*dC x dB*dC
%  (unnormalised) quantum states of Bob-Charlie. The last two dimensions
%  contain (a,x), the outcome and input of Alice respectively. 
%
% k: An integer, which specifies what size symmetric PPT extension to use
% in the corresponding part of the LHS model. That is, if k = 2, then these
% parts will be be required to have an (unnormalised) 2-symmetric PPT
% extension (as a relaxation of separability).
%  
% This function has one optional argument:
%  dims: a 2-D array containing the dimensions of Bob and Charlie. The
%  default is dB = dC = sqrt(dB*dC) (determined from sigma directly). This
%  argument only needs to be specified if Bob and Charlie do not hold
%  states of the same dimension.
%
% is_bLHS = bLHSTripartite1Unt(sigma,k) is the indicator function for
% whether sigma could have arisen from measurements on a bi-separable
% state or not. This is an outer-approximation, and uses the set of
% k-symmetric PPT extendible states as an outer approximation to the set of
% separable states. 
%
% is_bLHS = bLHSTripartite1Unt(sigma,k,dims) sets the local dimension of
% Bob to dims(1), and the local dimension of Charlie to dims(2). 
%
% [is_bLHS,F] = bLHSTripartite1Unt(sigma,k) returns in F the steering
% functional that certifies that the assemblage sigma could not have arisen
% from a bi-separable state. F is a 4-D array. The first two dimensions
% contain the dB*dC x dB*dC operators for Bob-Charlie. The last two
% dimensions are (a,x). A negative value certifies that a given assemblage
% demonstrates genuine multipartite steering. In the case that the
% assemblage has an LHS model, the empty array F = [] is returned.
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

if isa(sigma,'cvx') == 0 % if sigma isn't a CVX variable  
    % check that the assemblage is valid
    if  NSAssemblage(sigma) == 0
        error('assemblage is not valid')
    end 
end

cvx_begin sdp 

    variable sig(dB*dC,dB*dC,oa,ma) hermitian semidefinite
    % the assemblage for the A|BC parition
    variable pii(dB*dC,dB*dC,oa,ma) hermitian semidefinite
    % the assemblage for the B|AC partition
    variable gam(dB*dC,dB*dC,oa,ma) hermitian semidefinite
    % the assemblage for the C|AB partition
    variable slack(dB*dC,dB*dC,oa,ma) hermitian semidefinite
    % a slack variable, used to force an equality constraint below
    
    expression trBpii(dC,dC,oa,ma)
    % an expression holder for tr_B(Pi_a|x)
    expression trCgam(dB,dB,oa,ma)
    % an expression holder for tr_C(gamma_a|x)
    
    dual variable F
    % steering functional, the dual variable
    
    minimise real(sum(reshape(repmat(eye(dB*dC),[1,1,oa,ma]).*slack,1,[])))
    % minimise sum_ax tr slack_a|x 
    
    subject to
    
    F : sigma + slack == sig + pii + gam;
    % sigma_a|x + slack_a|x == sig_a|x + pi_a|x + gam_a|x the bLHS decomp
    
    for x = 1:ma
        for a = 1:oa
            trBpii(:,:,a,x) = PartialTrace(pii(:,:,a,x),1,[dB,dC]);
            % initialise the expression for tr_B(pi_a|x)
            trCgam(:,:,a,x) = PartialTrace(gam(:,:,a,x),2,[dB,dC]);
            % initialise the expression for tr_C(gam_a|x)
            
            SymmetricExtension(pii(:,:,a,x),k,[dB,dC],1,1) == 1;
            % pi_a|x should be separable. Relax to k-sym PPT extension
            SymmetricExtension(gam(:,:,a,x),k,[dB,dC],1,1) == 1;
            % gam_a|x should be separable. Relax to k-symm PPT extension
        end
    end

    LHSAssemblage(sig) == 1;
    % sig_a|x should have an LHS decomposition
    
    LHSAssemblage(trBpii) == 1;
    % tr_B(pi_a|x) should have an LHS decomposition
    
    LHSAssemblage(trCgam) == 1;
    % tr_C(gam_a|x) should have an LHS decomposition
    
cvx_end

% CVX will return a positive value if sigma does not have an fLHS model. It
% will return approx 0 (<= tol) if it is fLHS. This converts {>tol, <= tol}
% to {0,1}
is_bLHS = 1-(cvx_optval > tol);

% if the assemblage is fLHS, there is no inequality to return, so set it
% equal to the empty array.
if is_bLHS == 1
    F = [];
end

end