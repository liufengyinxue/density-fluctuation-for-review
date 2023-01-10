Scale = (1:10)' ;
MC = 1000 ;

load('Data.mat') ;
Data2003(Data2003(:,3) == -1,:) = [] ;
Data2003(Data2003(:,3) == -2,:) = [] ;
Data2004(Data2004(:,3) == -1,:) = [] ;
Data2004(Data2004(:,3) == -2,:) = [] ;
Data2005(Data2005(:,3) == -1,:) = [] ;
Data2005(Data2005(:,3) == -2,:) = [] ;
Data2006(Data2006(:,3) == -1,:) = [] ;
Data2006(Data2006(:,3) == -2,:) = [] ;
Data2007(Data2007(:,3) == -1,:) = [] ;
Data2007(Data2007(:,3) == -2,:) = [] ;

VMR_2003 = zeros(MC,length(Scale)) ;
VMR_2004 = zeros(MC,length(Scale)) ;
VMR_2005 = zeros(MC,length(Scale)) ;
VMR_2006 = zeros(MC,length(Scale)) ;
VMR_2007 = zeros(MC,length(Scale)) ;
Var_2003 = zeros(MC,length(Scale)) ;
Var_2004 = zeros(MC,length(Scale)) ;
Var_2005 = zeros(MC,length(Scale)) ;
Var_2006 = zeros(MC,length(Scale)) ;
Var_2007 = zeros(MC,length(Scale)) ;
Mean_2003 = zeros(MC,length(Scale)) ;
Mean_2004 = zeros(MC,length(Scale)) ;
Mean_2005 = zeros(MC,length(Scale)) ;
Mean_2006 = zeros(MC,length(Scale)) ;
Mean_2007 = zeros(MC,length(Scale)) ;
for i = 1:MC
    [VMR_2003(i,:),Mean_2003(i,:),Var_2003(i,:)] = Vector_DF_square_mark(Data2003,Scale,1000) ;
    [VMR_2004(i,:),Mean_2004(i,:),Var_2004(i,:)] = Vector_DF_square_mark(Data2004,Scale,1000) ;
    [VMR_2005(i,:),Mean_2005(i,:),Var_2005(i,:)] = Vector_DF_square_mark(Data2005,Scale,1000) ;
    [VMR_2006(i,:),Mean_2006(i,:),Var_2006(i,:)] = Vector_DF_square_mark(Data2006,Scale,1000) ;
    [VMR_2007(i,:),Mean_2007(i,:),Var_2007(i,:)] = Vector_DF_square_mark(Data2007,Scale,1000) ;
end
save('Withoutmale_TL_DF.mat') ;
