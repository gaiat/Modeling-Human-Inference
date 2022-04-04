
% This function computes the asymptotic Algorithmic Complexity per estimation for the single Sliding Window.
% NSP=1 --> Neural-Spiking P system implementation of the algorithm.
% NSP=0 --> Each operation counts 1.
% lopt=optimal run length

function[f]=AC_SA_estim(lopt,NSP)

if(NSP)

    Mwabs=2; % writing in memory
    Mrabs=2; % reading from memory
    Mmabs=2; % maintaining in memory
    Aabs=17; % addition
    Sabs=21; % subtraction
    Pabs=62; % product
    Dabs=179; % division
    Eabs=2; % exponential
    Fabs=0; % floor 
    OPPabs=Aabs; % opposite     
    Mw=Mwabs/Aabs;
    Mr=Mrabs/Aabs;
    Mm=Mmabs/Aabs;
    A=Aabs/Aabs;
    S=Sabs/Aabs;
    P=Pabs/Aabs;
    D=Dabs/Aabs;
    E=Eabs/Aabs;
    F=Fabs/Aabs;
    OPP=OPPabs/Aabs;   

    f = AC_SA_pred(lopt,NSP) -3*Mr -2*Mm - (S+A+2*P);

else   
    
    Mw=1; % writing in memory
    Mr=1; % reading from memory
    Mm=1; % maintaining in memory
    S=1; % additions+subtractions
    P=1; % product
    D=1; % division
    E=1; % exponential
    F=1; % floor 
    R=1; % square root
    OPP=1; % opposite 
       
    f = AC_SA_pred(lopt,NSP) -3*Mr -2*Mm - (2*S+2*P);

end


end




