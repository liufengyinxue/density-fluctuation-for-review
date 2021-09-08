function [VMR,Var,Mean] = Raster_TL_DF_Convolution(Data,Scale,Pbc,GPU) % ii jj
% This function computes density flucutation and Taylor's law
% based on the method of convolution (square box).
% This function can be executed on GPU.
% Input:
% Data: a raster, a matrix, and the elements of the matrix must be >= 0
% Scale: the sizes of the boxes, a vector
% Pbc: Indicates whether periodic boundary conditions are adopted or not, 
% generally, it is necessary to adopt for simulation data with periodic
% boundary conditions, turn on: 1 (periodic boundary condition), 
% turn off: 0 (open boundary condition).
% GPU: Indicates whether it is executed on the GPU, turn on: 1, turn off: 0
% Output:
% VMR: the density fluctuation, variance
% Var: Taylor's law, variance
% Mean: Taylor's law, mean

if GPU == 1
    DataTemp = gpuArray(Data) ;
elseif GPU == 0 
    DataTemp = Data ;
end

Var = zeros(length(Scale),1); 
Mean = zeros(length(Scale),1);
VMR = zeros(length(Scale),1); 
for i = 1:length(Scale)
    kernel = ones(Scale(i)) ; 
    kernel = single(kernel) ; 
    kersum = sum(kernel(:)) ;
    if GPU == 1
        kerneltemp = gpuArray(kernel) ; 
    elseif GPU == 0
        kerneltemp = kernel ;
    end
    if Pbc == 0
        DesityTemp = conv2(DataTemp,kerneltemp,'valid') ;
        if GPU == 1
            DesityTemp_C = gather(DesityTemp) ;
            Var(i,1) = var(DesityTemp_C(:)) ;
            Mean(i,1) = mean(DesityTemp_C(:)) ;
            VMR(i,1) = var(DesityTemp_C(:)/kersum) ;
        elseif GPU == 0
            Var(i,1) = var(DesityTemp(:)) ;
            Mean(i,1) = mean(DesityTemp(:)) ;
            VMR(i,1) = var(DesityTemp(:)/kersum) ;
        end
    elseif Pbc == 1
        padsize = Scale(i) - 1 ;
        DT = padarray(DataTemp,[padsize padsize],'circular') ;
        DesityTemp = conv2(DT,kerneltemp,'valid') ;
        if GPU == 1
            DesityTemp_C = gather(DesityTemp) ;
            Var(i,1) = var(DesityTemp_C(:)) ;
            Mean(i,1) = mean(DesityTemp_C(:)) ;
            VMR(i,1) = var(DesityTemp_C(:)/kersum) ;
        elseif GPU == 0
            Var(i,1) = var(DesityTemp(:)) ;
            Mean(i,1) = mean(DesityTemp(:)) ;
            VMR(i,1) = var(DesityTemp(:)/kersum) ;
        end
    end
end
% Path = 'E:\ZhenpengGe\New_FFM\' ;
% Filename = ['FFM_TL_DF_' num2str(ii) '_' num2str(jj) '_' num2str(5) '.mat'] ;
% save([Path Filename],'Var','Mean','VMR') ;
end