function rho = findPVMLHSStateGivenWitness(W,Max)
%FINDPVMLHSSTATEGIVENWITNESS searches for qubit-qudit states that have an
%LHS model for all projective measurements and violate a given entanglement
%witness
% This function has two required inputs:
%   W: an entanglement witness for a 2 x dB quantum state, such that
%   tr[W*rho^sep] >= 0, and such that tr[W*rho] < 0 for some entangled 2 x
%   dB quantum state.
%   Max: a 4-D array containing the POVM elements of a collection of qubit
%   projective measurements. The first two dimensions contain the
%   projectors, while the last two label (a,x).
%
%   rho = findPVMLHSStateGivenWitness(W,Max) returns a 2 x dB quantum state
%   that is both entangled and has a LHS model for all projective
%   measurements on Alice. If no such state can be found, the program
%   returns the empty array rho = [];
%
%  requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com),
%    vert2lcon (http://www.mathworks.com/matlabcentral/fileexchange/30892)
%  authors: Paul Skrzypczyk, Daniel Cavalcanti last updated: March 17,
%  2016

[dA,~,oa,ma] = size(Max);
% dA = dim. of Alice, oa = # outcomes for Alice, ma = # inputs for Alice
dB = length(W)/dA;
% dB = dim. of Bob;

r = findRadiusPolytopeInBlochSphere(Max);
% r = radius of largest ball that fits inside the bloch sphere. 
% it is the shrinking factor that we have to apply to the state.

cvx_begin sdp quiet

    variable O(dA*dB,dA*dB) hermitian
    % O is the quasi-state that we optimise over 
    expression Op(dA*dB,dA*dB)
    % Op is an expression which stores O after the application of the noise
    
    Op = r*O + (1-r)*Tensor(eye(dA)/dA,PartialTrace(O,1,[dA dB]));
    % Op = r*O_AB + (1-r) id_A/dA otimes O_B
    
    minimise real(trace(W*Op))
    % minimise trace[W*Op]

    LHSAssemblage(genAssemblage(O,Max),1) == 1;
    % O should have an LHS model for the measurements Max
    
    trace(Op) == 1;
    % Op should be normalised
    
    Op == hermitian_semidefinite(dA*dB);
    % Op should be PSD
    
 cvx_end
 
 % if the state Op is not entangled (or the entanglement is not detected by
 % W, then we have not found an example.
 if cvx_optval >= 0
     rho = [];
 end
 
end
    