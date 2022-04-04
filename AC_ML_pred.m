
% This function computes the asymptotic Algorithmic Complexity per prediction for the Memoryless model.
% NSP=1 --> Neural-Spiking P system implementation of the algorithm.
% NSP=0 --> Each operation counts 1.


function[f]=AC_ML_pred(NSP)

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

    f = (2*A+2*S+3*P) + 6*Mr + 3*Mm ;

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

    f = (4*S+3*P) + 6*Mr + 3*Mm ;

end


end


