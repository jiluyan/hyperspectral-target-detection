function [MF,wmf,energy]=jily_MF(Y,d)
% matched filter, for hyperspectral/multispectral target detection
%Geng, Xiurui, Luyan Ji, and Kang Sun. 2016. "Clever eye algorithm for target detection of remote sensing imagery." 
%ISPRS Journal of Photogrammetry and Remote Sensing 114:32-9. doi: https://doi.org/10.1016/j.isprsjprs.2015.10.014.
%
%Input: Y: the data matrix, L*N, L is the number of bands, N is the number
%of pixels
%       d:the target vector, L*1
%Output: MF: MF output, 1*L
%        wmf: mf detector
%        energy: average filter output energy of mf
% Luyan Ji, jily@mail.ustc.edu.cn
%2023.7.18

m=mean(Y,2);

[~,nN]=size(Y);

Y1=Y-repmat(m,[1,nN]);
K=Y1*Y1'/nN;

wmf=inv(K)*(d-m)/((d-m)'*inv(K)*(d-m));
MF=wmf'*Y1;
%'MF erengy'
energy=MF*MF'/nN;


end