function [VMR,Mean,Var] = Vector_TL_DF_square_mark(Points,Scale,N)
% This function computes the ratio between variance and 
% mean of vector data based on the method of random boxes (square box) 
% This function is coded by Mingji Huang, SJTU.
% This function is used for two dimensional data, 
% and it is easy to extend it to higher dimensions.
% This fuction can handle boundary well.
% input:
% Data, rows are samples and column are properties, a matrix
% Scale, the sizes of the boxes, must be a column vector
% N, the number of times of throwing boxes, a number
% output,
% VMR, the ratio between variance and mean, a vector
% Var: Taylor's law, variance
% Mean: Taylor's law, mean

Data = Points(:,1:2) ;
DMax = max(Data) ;
DMin = min(Data) ;
L = DMax - DMin ; %system size
% [Samples,~] = size(Data) ; % sample number
v = zeros(length(Scale),N) ;
% rho0 = Samples/(L(1)*L(2)) ; % point density of this area
for ii = 1:length(Scale)
    for jj = 1:N
        %tp = abs(bsxfun(@minus,Data,Scale(ii)/2+(L-Scale(ii))*rand(1,2))) ;
        tp = abs(bsxfun(@minus,Data,Scale(ii)/2+[(L(1)-Scale(ii))*rand,(L(2)-Scale(ii))*rand])) ; % elegant way of setting closed boundary
%         v(1,ii,jj) = (sum(prod(tp<Scale(ii)/2,2).*Points(:,3))/Scale(ii)^2-rho0)^2 ; % tp<Scale(ii)/2: two direction is shorter than box scale
        v(ii,jj) = sum(prod(tp<Scale(ii)/2,2).*Points(:,3)) ;
    end
end
Square_scales = repmat(Scale,1,N) ;
VMR = var(v./Square_scales.^2,0,2) ;
Mean = mean(v,2) ;
Var = var(v,0,2) ;