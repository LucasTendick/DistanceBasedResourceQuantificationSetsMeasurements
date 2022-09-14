function [Pgrho,optMax] = globalSteeringGuessProbState(rho,Max,Mb,varargin)
%GLOBALSTEERINGGUESSPROBSTATE estimate the one-sided device-independent
%global guessing probability of a quantum state
% This function has three required arguments:
%  rho: a 2-D array containing a dA*dB x dA*dB bipartite quantum state rho.
%  Max: a 4-D array containing an initial set of POVMs to start the
%  algorithm with. The first two dimensions contain the dA x dA POVM
%  elements. The last two dimensions label (a,x). 
%  Mb: a 3-D array containing the POVM elements of the (trusted)
%  measurement of Bob. The first two dimensions contain the dB x dB POVM
%  elements, while the last dimension is b, the outcome label.
%
% This function has one optional argument:
%  improv_tol: (default: 1E-6) a scalar which sets the threshold
%  improvement in one round of the see-saw iteration before the algorithm
%  terminates.
%
% Pgrho = globalSteeringGuessProbState(rho,Max,Mb) returns an upper bound
% on the one-sided device-independent global guessing probability of the
% state rho, when measured with ma measurements of oa outcomes each, taking
% as an initial seed for the measurements Max. It employs a see-saw
% algorithm to try and find the optimal set of measurements for the state
% rho.
%
% [Pgrho,optMax] = globalSteeringGuessProbState(rho,Max,Mb) also returns
% the final set of measurements in optMax, that achieve the global guessing
% prob Pgrho. optMax is a 4-D array. The first two dimenions contain the dA
% x dA POVM elements. The last two dimensions contain (a,x).
%
% [Pgrho,optMax] = globalSteeringGuessProbState(rho,Max,Mb,improv_tol) uses
% the tolerance improv_tol in its stopping criteria for the algorithm: only
% if the improvement in one round of the see-saw algorithm drops below
% improv_tol will the algorithm terminate.
%
% requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
% authors: Paul Skrzypczyk, Daniel Cavalcanti
% last updated: March 17, 2016

[improv_tol] = opt_args({1E-6},varargin{:});
% improv_tol: tolerance in see-saw. default is 1e-6.

[Pgi,F] = globalSteeringGuessProb(genAssemblage(rho,Max),Mb);
% find the initial global guess prob Pgi, and the initial steering
% inequality F which certifies Pgi

improv = 1;
% initalise improv > improv_tol.

while improv >= improv_tol % as long as the improvement is big enough
    [~, Max] = bestSteeringMeasurements(rho,-F);
    % find the best set of POVMs given the state rho, and functional F
    [Pgf,F] = globalSteeringGuessProb(genAssemblage(rho,Max),Mb);
    % find the new global guessing prob and function F given rho and Max
    improv = Pgi - Pgf;
    % calculate the improvement in the Pg
    Pgi = Pgf;
    % set initial = final, and iterate.
end

Pgrho = Pgf;
optMax = Max;

end