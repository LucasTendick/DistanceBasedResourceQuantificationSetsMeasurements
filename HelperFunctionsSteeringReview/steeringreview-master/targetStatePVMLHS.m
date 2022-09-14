function is_PVM_LHS = targetStatePVMLHS(rho,Max)
%TARGETSTATEPVMLHS determine whether a qubit-qudit state has an LHS model
%for all projective measurements
% This function has two required inputs:
%   rho: a 2-D array containing a 2 x dB bipartite quantum state.
%   Max: a 4-D array containing the POVM elements of a collection of qubit
%   projective measurements. The first two dimensions contain the
%   projectors, while the last two label (a,x).
%
%  is_PVM_LHS = targetStatePVMLHS(rho,Max) returns 1 if the set of
%  measurements Max are able to demonstrate that the qubit-qudit state rho
%  has an LHS model for all qubit projective measurements. It returns 0 if
%  this is not the case. Note that in the latter case one cannot draw any
%  conclusions about whether or not rho has an LHS model.
%
%  requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com),
%    vert2lcon (http://www.mathworks.com/matlabcentral/fileexchange/30892)
%  authors: Paul Skrzypczyk, Daniel Cavalcanti last updated: March 17,
%  2016

[dA,~,oa,ma] = size(Max);
% dA = dim. of Alice, oa = # outcomes for Alice, ma = # inputs for Alice
dB = length(rho)/dA;
% dB = dim. of Bob;

r = findRadiusPolytopeInBlochSphere(Max);
% r = radius of largest ball that fits inside the bloch sphere. 
% it is the shrinking factor that we have to apply to the state.

cvx_begin sdp quiet

    variable O(dA*dB,dA*dB) hermitian
    % O is the quasi-state that we optimise over
    
    subject to
    
    rho == r*O + (1-r)*Tensor(eye(dA)/dA,PartialTrace(O,1,[dA dB]));
    % rho_AB = r*O_AB + (1-r) Id_A/dA otimes O_B

    LHSAssemblage(genAssemblage(O,Max),1) == 1;
    % O_AB should have an LHS model for the measurements Max
        
 cvx_end
 
 % CVX will return +inf if the problem is infeasible, and 0 if feasible
 % this maps {+inf,0} to {0,1}
 is_PVM_LHS = 1-min(cvx_optval,1);
 
end