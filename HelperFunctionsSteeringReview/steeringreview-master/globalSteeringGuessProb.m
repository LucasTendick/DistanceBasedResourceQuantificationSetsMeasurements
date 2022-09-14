function [Pg, F] = globalSteeringGuessProb(sigma,Mb)
%GLOBALSTEERINGGUESSPROB one-sided-device-independent global guessing
%probability
% This function has two required arguments:
%  sigma: a 4-D array containing a bipartite assemblage. The first two
%  dimensions contain the dB x dB quantum states, while the last two
%  dimensions are (a,x).
%  Mb: a 3-D array containing the POVM elements of the (trusted)
%  measurement of Bob. The first two dimensions contain the dB x dB POVM
%  elements, while the last dimension is b, the outcome label.
%
% Pg = globalSteeringGuessProb(sigma,Mb) returns the one-sided
% device-independent global guessing probability that Eve will guess the
% pair of outcomes (a,b), when Alice performs measurement x=0, and Bob
% performs the measurement Mb.
%
% [Pg, F] = globalSteeringGuessProb(sigma,Mb) also returns the steering
% function F which certifies that the one-sided device-independent global
% guessing probability is Pg. F is a 4-D array. The first two dimensions
% contain the dB x dB operators which form the steering function. The last
% two dimensions are (a,x).
%
% requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
% authors: Paul Skrzypczyk, Daniel Cavalcanti
% last updated: March 17, 2016

[dB,~,oa,ma] = size(sigma);
% dB = dim. of Bob, oa = # outcomes for Alice, ma = # inputs for Alice

[~,~,ob] = size(Mb);
% ob = # outcomes for Bob

% check that the assemblage is valid
if  NSAssemblage(sigma) == 0
    error('assemblage is not valid')
end

% check that the measurement is valid
if validPOVMs(Mb) == 0
    error('Measurement is not valid')
end

PgObj = permute(repmat(eye(oa),[1,1,dB,dB,ob]),[3,4,1,2,5])...
    .*permute(repmat(Mb,[1,1,1,oa,oa]),[1,2,4,5,3]);
% PgObj = delta_(a,e)*M_(e')

cvx_begin sdp quiet
    
    variable sigaeex(dB,dB,oa,oa,ob,ma) hermitian semidefinite
    % sigma_a|x^ee' is the assemblage prepared by Eve.
    
    dual variable F
    % the steering function is the dual variable.
    
    maximise real(sum(reshape(conj(PgObj).*sigaeex(:,:,:,:,:,1),1,[])))
    
    subject to
    
    F : sigma == squeeze(sum(sum(sigaeex,5),4));
    % sigma_a|x = sum_ee' sigma_a|x^ee'. F is the dual variable to this
    % constraint.
        
    for e = 1:oa
        for ep = 1:ob
            NSAssemblage(squeeze(sigaeex(:,:,:,e,ep,:))) == 1;
            % sigma_a|x^ee' is a valid assemblage for all e,e'
        end
    end
    
cvx_end

Pg = cvx_optval;
    
end