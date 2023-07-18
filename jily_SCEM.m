function SCEM=jily_SCEM(Y,D)
% sum constrained energy minimization (SCEM), for hyperspectral/multispectral muti-target target detection

%Ren, H.; Du, Q.; Chang, C.I.; Jensen, J.O. Comparison between constrained energy minimization based approaches for 
%hyperspectral imagery. In Proceedings of the Advances in Techniques for Analysis of Remotely Sensed Data, 2003 IEEE 
%Workshop on, Oct 2003, pp. 244¨C248. https://doi.org/10.1109/WARSD.2003.1295199.
%
%Input: Y: the data matrix, L*N, L is the number of bands, N is the number
%of pixels
%       D:the target matrix, L*M,Mis the number of targets

%Output: SCEM: SCEM output, 1*L

% Luyan Ji, jily@mail.ustc.edu.cn
%2023.7.18
[L,M]=size(D);
[~,nN]=size(Y);

SCEM=zeros(1,nN);

for i=1:M
    d=D(:,i);
    CEM=jily_CEM(Y,d);
    SCEM=SCEM+CEM;
end

end