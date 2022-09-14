function sigma = genAssemblage(rho,Max)
%GENASSEMBLAGE Create an assemblage from a state and a set of measurements
%  This function has two required inputs:
%   rho: a 2-D square array, containing a valid bipartite quantum state
%   (i.e. a positive semidefinite operator). The state can be unnormalised
%   if desired.
%   Max: a 4-D array, containing the set of POVMs that Alice will measure
%   on the state. The first two dimensions contain the POVM elements, and
%   the last two dimensions label (a,x). That is, Max(:,:,a,x) = M_a|x, the
%   dA x dA POVM element corresponding to input x and outcome a. 
%
%  sigma = NSAssemblage(sigma) generates the assemblage prepared for Bob,
%  when Alice performs measurements Max on the bipartite state rho
%
% sigma = genAssemblage(rho,Max) can also be used inside CVX as 
% a partially specified problem, to generate the assemblage sigma from the
% CVX variable rho and the set of measurements Max
%
% EXAMPLE:
%   cvx_begin
%
%       variable rho(dA*dB,dA*dB) hermitian semidefinite
%       expression sigma(dB,dB,oa,ma) 
%
%       subject to
%
%       sigma = genAssemblage(rho,Max);
%
%   cvx_end
%
%  Inside CVX, sigma is a variable corresponding to a quantum state, and
%  sigma is an expression which contains the variables corresponding to the
%  members of the assemblage.
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti
%   last updated: March 17, 2016

[dA,~,oa,ma] = size(Max);
% dA = dim. of Alice , oa = # outcomes, ma = # inputs for Alice
dB = length(rho)/dA;
% dB = dim. of Bob

if isa(rho,'cvx') == 0 % if rho isn't a CVX variable
    
    sigma = zeros(dB,dB,oa,ma);
    
    for a = 1:oa
        for x = 1:ma
            sigma(:,:,a,x) = ...
                PartialTrace(Tensor(Max(:,:,a,x),eye(dB))*rho,1,[dA,dB]);
        end
    end
    
else
    
    cvx_begin sdp quiet
    
        expression sigma(dB,dB,oa,ma)
        % specificy sigma as an expression (since we don't want to increase
        % the number of variables inside the SDP)
        
        for a = 1:oa
            for x = 1:ma
                sigma(:,:,a,x) ...
                    = PartialTrace(Tensor(Max(:,:,a,x),eye(dB))*rho,1,[dA,dB]);
            end
        end
    
    cvx_end
    
end


end