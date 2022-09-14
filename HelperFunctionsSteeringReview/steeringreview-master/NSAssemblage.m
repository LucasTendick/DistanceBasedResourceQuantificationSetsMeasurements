function is_NS_assemblage = NSAssemblage(sigma,varargin)
%NSASSEMBLAGE Determines whether an assemblage is valid or not
%  This function has one required argument:
%   sigma: a 4-D array, containing the members of the assemblage. The first 
%   two dimensions contain the (unnormalised) quantum states, while the
%   remaining two dimensions are (a,x), such that sigma(:,:,a,x) =
%   \sigma_a|x. 
%
%  is_NS_assemblage = NSAssemblage(sigma) is the indicator function for
%  non-signalling assemblages. It returns 1 if the assemblage is a valid
%  unnormalised assemblage, and 0 otherwise.
%
% This function has one optional argument:
%   nm: (default 0)
%
% is_NS_assemblage = NSAssemblage(sigma,nm) is the indicator function for
% non-signalling assemblages that are additionally required to be
% normalised when nm = 1. In this case, trace(sum(sigma(:,:,:,x),3)) = 1. 
%
% is_NS_assemblage = NSAssemblage(sigma,nm) can also be used inside CVX as 
% a partially specified problem, to enforce the constraint that the CVX
% variable sigma should be an unnormalised assemblage if nm = 0, and a 
% normalised one otherwise. 
%
% EXAMPLE:
%   cvx_begin
%
%       variable sigma(dB,dB,oa,ma) hermitian semidefinite
%
%       subject to
%
%       NSAssemblage(sigma) == 1
%
%   cvx_end
%
%  Inside CVX, sigma is an assemblage with ma inputs and oa outcomes for
%  Alice, preparing quantum states of dimension dB x dB for Bob. CVX
%  enforces that sigma should be a valid non-signalling assemblage, i.e.
%  that it should satisfy all the no-signalling constraints. 
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti
%   last updated: March 17, 2016


tol = 1e-10; 
% numerical tolerance used

[nm] = opt_args({0},varargin{:}); 
%if unspecified, it is assumed that the assemblage is unnormalised.

[dB,~,oa,ma] = size(sigma);
% dB = dim. of Bob, oa = # outcomes for Alice, ma = # inputs for Alice


if isa(sigma,'cvx') == 0 % if sigma isn't a CVX variable
    
    % first check that each member of the assemblage is
    % positive-semidefinite
    for a = 1:oa
        for x = 1:ma
            if ~all(eig(sigma(:,:,a,x))>= -tol)
                is_NS_assemblage = 0;
                display('assemblage not positive semidefinite')
                return
            end
        end
    end
    
    % now check that the assemblage satisfies the no-signalling constraints
    if ~all(reshape(abs(squeeze(sum(sigma(:,:,:,2:ma),3))...
            -repmat(sum(sigma(:,:,:,1),3),[1 1 ma-1])) <= tol,1,[]))
        is_NS_assemblage = 0;
        display('assemblage is signalling')
        return
    end
    
    % if the assemblage should be normalised, check that the trace of the
    % reduced state of Bob is equal to unity
    if nm == 1
        if ~(abs(trace(sum(sigma(:,:,:,1),3)) - 1) <= tol)
            display('assemblage is not normalised')
            is_NS_assemblage = 0;
            return
        end
    end
    
    % if all the above are satisfied, the assemblage is valid
    is_NS_assemblage = 1;
    return

elseif isa(sigma,'cvx') == 1 % if sigma is a CVX variable
    % begin the partial specification
    
    cvx_begin
    
    subject to
        
        % the assemblage should satisfy the no-signalling constraints
        for x = 2:ma
            squeeze(sum(sigma(:,:,:,x),3)) == squeeze(sum(sigma(:,:,:,1),3))
        end
        
        % enforce normalisation if it is required.
        if nm == 1
            trace(sum(sigma(:,:,:,1),3)) == 1
        end
        
    cvx_end
    
    % CVX will return +inf if the problem is infeasible, and 0 if feasible
    % this maps {+inf,0} to {0,1}
    is_NS_assemblage = 1-min(cvx_optval,1);
    
end

end


    