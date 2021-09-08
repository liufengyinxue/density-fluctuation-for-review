function Ricker_model(L,T,Burntime,lambda,cb)
%%%
%%%

% Parameters
N = L ; % spatial length
LF = 0.1 ;
z = 4 ;
r = 2.3 ;
B = lambda ;
GF = 1 - LF ;

% Data save
Path = '/mnt/9768e637-414c-423b-9fa1-ed0a52e40b62/RickerL512/' ;
Sample_interval = 1000 ;
c0 = 0 ;
c1 = 0 ;

% Initial state
Ricker = zeros(N) ;
MIU = 0.5 ;
SIGMA = 0.1 ;
[Ricker] = Initial_state(Ricker,MIU,SIGMA) ;

% For Neighbours, tricky
North = [1;0;0] ;
South = [0;1] ;
East = [0,1] ; 
West = [1,0,0] ;

% Loop
for k = 1:(Burntime + T + 1)
    % generate noise
    Noise = randn(N) ;
    % extract Neighbours with periodic boundary condition
    Ricker_North = imfilter(Ricker,North,'circular') ;
    Ricker_South = imfilter(Ricker,South,'circular') ;
    Ricker_East = imfilter(Ricker,East,'circular') ;
    Ricker_West = imfilter(Ricker,West,'circular') ;
    Noise_North = imfilter(Noise,North,'circular') ;
    Noise_South = imfilter(Noise,South,'circular') ;
    Noise_East = imfilter(Noise,East,'circular') ;
    Noise_West = imfilter(Noise,West,'circular') ;
    % Global influence
    [Global] = Noisy_Ricker_map(Ricker,Noise,r,B) ;
    % Local influence
    [Local_North] = Noisy_Ricker_map(Ricker_North,Noise_North,r,B) ;
    [Local_South] = Noisy_Ricker_map(Ricker_South,Noise_South,r,B) ;
    [Local_East] = Noisy_Ricker_map(Ricker_East,Noise_East,r,B) ;
    [Local_West] = Noisy_Ricker_map(Ricker_West,Noise_West,r,B) ;
    Ricker = GF*Global + (LF/z)*(Local_North + Local_South + Local_East + Local_West) ;
    if (k > Burntime) && (rem(k - Burntime,Sample_interval) == 0)
        % pre
        c0 = c0 + 1 ;
        Filename0 = ['Ricker_' num2str(cb) '_' num2str(c0) '_0.mat'] ;
        save([Path Filename0],'Ricker') ;
    elseif (k > Burntime + 1) && (rem(k - Burntime,Sample_interval) == 1)
        % post
        c1 = c1 + 1 ;
        Filename1 = ['Ricker_' num2str(cb) '_' num2str(c1) '_1.mat'] ;
        save([Path Filename1],'Ricker') ;
    end
end
end

    function [S_post] = Initial_state(S_pre,miu,sigma)
        [m,n] = size(S_pre) ;
        S_post = zeros(m,n) ;
        for i = 1:m
            for j = 1:n
                R1 = normrnd(miu,sigma) ; 
                while R1 < 0
                    R1 = normrnd(miu,sigma) ; 
                end
                S_post(i,j) = R1 ;
            end
        end
    end

    function [MatC] = Noisy_Ricker_map(MatA,MatB,Par1,Par2)
        MatC = MatA.*exp(Par1*(1 - MatA)).*(1 + Par2*MatB) ;
    end