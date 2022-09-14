function r = findRadiusPolytopeInBlochSphere(Max)
%FINDRADIUSPOLYTOPEINBLOCHSPHERE finds radius of biggest ball which fits inside the polytope
% This function has two required inputs:
%   Max: a 4-D array, containing the set of qubit projective measurements.
%   That is, Max is a 2 x 2 x 2 x ma array, comprising ma measurements,
%   with 2 outcomes each, and each projector being a 2 x 2 matrix. 
%
% r = findRadiusPolytopeInBlochSphere(Max) calculates the radius r of the
% largest ball which fits inside the polytope defined by the bloch vectors
% corresponding to each of the projectors of the set of measurements Max.
%
%   requires: QETLAB (http://www.qetlab.com), 
%     vert2lcon (http://www.mathworks.com/matlabcentral/fileexchange/30892)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti
%   last updated: March 17, 2016

[dA,~,oa,ma] = size(Max);
% dA = dim., oa = # outcomes, ma = # inputs for Alice

if dA ~= 2 || oa ~= 2
    error('This only works for qubits')
end

r = zeros(ma*oa,3); %initialise array of Bloch vectors
cnt = 1;
for x = 1:ma
    for a = 1:oa
        for i = 1:3
            r(cnt,i) = real(trace(Pauli(i)*Max(:,:,a,x)));
            % components of Bloch vector are inner product with Paulis. 
        end
        cnt = cnt+1;
    end
end

[~,b,~,~]=vert2lcon(r); 
% find list of distances from the origin of the facets of the polytope
r = min(b); % the radius of the ball is the distance to the closest facet

end