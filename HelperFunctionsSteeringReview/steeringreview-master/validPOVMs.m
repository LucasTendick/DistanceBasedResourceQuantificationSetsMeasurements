function are_valid_POVMs = validPOVMs(Max)
%VALIDPOVMS Determines whether a set of POVMs is valid or not
%  This function has one argument:
%   Max: a 4-D array, containing the POVM elements. The first two
%   dimensions contain the POVM elements, while the remaining two
%   dimensions are (a,x), such that Max(:,:,a,x) = M_a|x.
%
%  are_valid_POVMs = validPOVMs(Max) is the indicator function for valid
%  POVMs. It returns 1 if the measurements are valid , and 0 otherwise.
%
% are_valid_POVMs = validPOVMs(Max) can also be used inside CVX as 
% a partially specified problem, to enforce the constraint that the CVX
% variable Max should be a valid collection of POVMs.
%
% EXAMPLE:
%   cvx_begin
%
%       variable Max(dA,dA,oa,ma)
%
%       subject to
%
%       validPOVMs(Max) == 1
%
%   cvx_end
%
%  Inside CVX, Max is a set of measurements with ma inputs and oa outcomes
%  of dimension dA x dA. CVX enforces that Max should be a valid set of
%  measurements, i.e. that for each x, sum_a M_a|x == identity(dA).
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti
%   last updated: March 17, 2016

tol = 1e-10; % numerical tolerance used

[dA,~,oa,ma] = size(Max);
% dA = dim., oa = # outcomes, ma = # inputs for Alice

if isa(Max,'cvx') == 0 % if Max isn't a CVX variable
    
    % first check that each POVM element is
    % positive-semidefinite
    for a = 1:oa
        for x = 1:ma
            if ~all(eig(Max(:,:,a,x))>= -tol)
                are_valid_POVMs = 0;
                display('POVM elements are not positive semidefinite')
                return
            end
        end
    end
    
    % now check that the POVMs sum to the identity
    if ~all(reshape(abs(squeeze(sum(Max,3))-repmat(eye(dA),[1,1,ma])) ...
            <= tol,1,[]))
        are_valid_POVMs = 0;
        display('POVMs do not sum to identity')
        return
    end
    
    % if all the conditions are satisfied, the POVMs are valid.
    are_valid_POVMs = 1;
    return

elseif isa(Max,'cvx') == 1 % if sigma is a CVX variable
    
    cvx_begin
    
    subject to
        
        squeeze(sum(Max,3)) == repmat(eye(dA),[1,1,ma]);
        % sum_a M_a|x == identity(dA) forall x
        
    cvx_end
    
    % CVX will return +inf if the problem is infeasible, and 0 if feasible
    % this maps {+inf,0} to {0,1}
    are_valid_POVMs = 1-min(cvx_optval,1);
    
end


    