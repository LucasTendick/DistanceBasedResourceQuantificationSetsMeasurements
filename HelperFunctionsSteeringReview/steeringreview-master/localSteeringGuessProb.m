function [Pg, F] = localSteeringGuessProb(sigma)
%LOCALSTEERINGGUESSPROB one-sided-device-independent local guessing
%probability
% This function has one required argument:
%  sigma: a 4-D array containing a bipartite assemblage. The first two
%  dimensions contain the dB x dB quantum states, while the last two
%  dimensions are (a,x).
%
% Pg = localSteeringGuessProb(sigma) returns the one-sided
% device-independent local guessing probability that Eve will guess the
% outcomes a, when Alice performs measurement x=0.
%
% [Pg, F] = localSteeringGuessProb(sigma) also returns the steering
% function F which certifies that the one-sided device-independent local
% guessing probability is Pg. F is a 4-D array. The first two dimensions
% contain the dB x dB operators which form the steering function. The last
% two dimensions are (a,x).
%
% requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
% authors: Paul Skrzypczyk, Daniel Cavalcanti
% last updated: March 17, 2016

[dB,~,oa,ma] = size(sigma);
% dB = dim. of Bob, oa = # outcomes for Alice, ma = # inputs for Alice

% check that the assemblage is valid
if  NSAssemblage(sigma) == 0
    error('assemblage is not valid')
end

PgObj = permute(repmat(eye(oa),[1,1,dB,dB]),[3,4,1,2])...
    .*repmat(eye(dB),[1,1,oa,oa]);
% PgObj = delta_(a,e)*eye(dB)

cvx_begin sdp quiet
    
    variable sigaex(dB,dB,oa,oa,ma) hermitian semidefinite
    % sigma_a|x^e is the assemblage prepared by Eve.
    
    dual variable F
    % the steering function is the dual variable.
    
    maximise real(sum(reshape(PgObj.*sigaex(:,:,:,:,1),1,[])))
    
    subject to
    
    F : sigma == squeeze(sum(sigaex,4));
    % sigma_a|x = sum_e sigma_a|x^e. F is the dual variable to this
    % constraint.
        
    for e = 1:oa
        NSAssemblage(squeeze(sigaex(:,:,:,e,:))) == 1;
        % sigma_a|x^e is a valid assemblage for all e,e'
    end
    
cvx_end

Pg = cvx_optval;
    
end