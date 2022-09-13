function [CoherenceDual,Cax,rho_x,L_xi] = CoherenceDiamondNormDual(Max,p)


% This function is based on the the work 'Distance-based resource
% quantification for sets of quantum measurements' (see: arXiv:2205.08546)


% The function CoherenceDiamondNormDual quantifies the 
% coherence according to Eq. (19) respectively SDP Eq. (83) (or Eq.
% (F5).

% Inputs: 
% 1) a 4-D array, containing the POVM elements. The first two
% dimensions contain the POVM elements (the Matrix representation),
% while the remaining two dimensions are (a,x), such that Max(:,:,a,x) = M_a|x.
% The outcome 'a' is the 3rd and the setting 'x' the 4th dimension.

% 2) a probability vector denoting the probabilities for choosing the
% settings x with probabbility p(x). For example p = ones(1,m)/m is
% a uniform distribution for m settings

%Outputs:
% 1) CoherenceDual, the value of the coherence of M_a|x
% 2) Cax, hyperplane coefficients (matrcies that certify the
% coherence)
% 3) states rho_x that fix the scale of the coherence inequality
% 4) classical bounds, L_xi that together give together the total bound
% \sum_{x,i} L_xi = \sum_{x,i} L_{x,i}


% Requires: CVX (http://cvxr.com/cvx/), A solver like SDPT3, Sedumi or
% Mosek, etc., QETLAB (http://www.qetlab.com) 


% Authors: Lucas Tendick, Martin Kliesch, Hermann Kampermann, Dagmar Bruß
% Code by: Lucas Tendick
% Last Update: 01.09.2022

% Note: The number of outcomes is assumed to be equal for each setting.
% Function can be used within CVX (as a convex constraint)


[d,~,o,m] = size(Max); %Get the dimension d, the number of outcomes o and 
                       %the number of settings m

% Set uo the incoherent basis
Id = eye(d); 

for n = 1:d
    ket_i(:,n) = Id(:,n); % Get kets of the comp. bases
end

                       
                       
cvx_begin

variable L_xi(d,m)
% Variables for the classical bounds

variable rho_x(d,d,m) hermitian semidefinite
% Variables for the states that fix the scale of the 
% coherence inequality

variable Cax(d,d,o,m) hermitian semidefinite
% Variables for the hyperplane coefficients (matrices) 
% that certify the coherence

expression S(1)
% Place holder for the objective function

expressions class_strategy(o,d,m)
% Place holder for the classical/incoherent strategies

%  Constraints

for x = 1:m
    trace(rho_x(:,:,x)) == 1;
end

for x = 1:m
    for a = 1:o
        S = S+p(x)*trace(Max(:,:,a,x)*Cax(:,:,a,x));
        rho_x(:,:,x)-Cax(:,:,a,x) == hermitian_semidefinite(d);
    end
end

   for x = 1:m
        for a = 1:o
            for i = 1:d
            class_strategy(a,i,x) = p(x)*trace(Cax(:,:,a,x)*ket_i(:,i)*ket_i(:,i)');
            L_xi(i,x)-class_strategy(a,i,x) >= 0;
            end
        end
  end
    

Obj = real(S-sum(sum(L_xi)));

maximize real(Obj)

cvx_end

CoherenceDual = real(Obj);
end