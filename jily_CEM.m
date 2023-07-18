function [CEM,wcem,energy]=jily_CEM(Y,d)
%constrained energy minimization(CEM), for hyperspectral/multispectral target detection

%Geng, Xiurui, Luyan Ji, and Kang Sun. 2016. "Clever eye algorithm for target detection of remote sensing imagery." 
%ISPRS Journal of Photogrammetry and Remote Sensing 114:32-9. doi: https://doi.org/10.1016/j.isprsjprs.2015.10.014.

%Input: Y: the data matrix, L*N, L is the number of bands, N is the number
%of pixels
%       d:the target vector, L*1
%Output: CEM: CEM output, 1*L
%        wcem: cem detector
%        energy: average filter output energy of cem
% Luyan Ji, jily@mail.ustc.edu.cn
%2023.7.18


[~,nN]=size(Y);
R=Y*Y'/nN;
%(d'*inv(R)*d)
wcem=inv(R)*d/(d'*inv(R)*d);
CEM=wcem'*Y;
%'CEM erengy';
energy=CEM*CEM'/nN;

end