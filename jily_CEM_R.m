function [wcem]=jily_CEM_R(R,d)

%Input: R is an L*L correlation matrix, 
%       d is the spectrum of the target
%Output: wcem: the weight of CEM

AA=inv(R)*d; 
wcem=AA/(d'*AA); 
% CEM=wcem'*Y;

end

