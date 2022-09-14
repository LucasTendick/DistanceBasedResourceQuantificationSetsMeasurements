function Max = genRandProjMeas(dA,ma)
% GENRANPROJMEAS Generates a set of random rank-1 projective measurements
%  This function has two required inputs:
%   dA: the dimension of the Hilbert space on which the measurements act
%   ma: the number of measurements
%
% Max = genRandProjMeas(dA,ma) generates a 4-D array containing the rank-1
% projective elements of the ma measurements. The first two dimensions
% contain the measurement operators, while the last two label (a,x). That
% is, Max(:,:,a,x) = Pi_a|x, where (Pi_a|x)^2 = Pi_a|x, and rank(Pi_a|x) =
% 1. The measurements are drawn from the Haar measure.
%
%   requires: QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti
%   last updated: March 17, 2016

Max = zeros(dA,dA,dA,ma);
%initalise Max

for x = 1:ma
    temp = RandomUnitary(dA);
    % pick a dA x dA unitary matrix from the Haar measure
    [V,~] = eigs(temp); % find its eigenvectors
    for a = 1:dA-1
        vec = V(:,a); 
        Max(:,:,a,x) = vec*vec'; %projectors onto eigenvectors of U
    end
    Max(:,:,dA,x) = eye(dA) - sum(Max(:,:,1:dA-1,x),3);
    % for stability, we pick the last POVM element to be identity - sum of
    % the others, as sometimes numerical errors cause problems if we assign
    % each projector directly from the random unitary
end

end
        
    

