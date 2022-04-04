
% This function computes the asymptotic Algorithmic Complexity per prediction for the Mixture of Sliding Windows with N>1 units.
% NSP=1 --> Neural-Spiking P system implementation of the algorithm.
% NSP=0 --> Each operation counts 1.
% lopt=Nx1 vector of optimal run lengths

function[f]=AC_NSAs_pred(lopt,N,NSP) 

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
    
    f = (2+N)*Mm + Mw + Mm*(max(floor(lopt))+1) + Mm*N + Mm*N + 2*Mm*(N^2) + (6*Mr + 3*S +A + 3*P + Mw + Mm)*N + ...
        (Mm+Mr)*N + (3*Mr+2*S+OPP+5*P+2*D+E+(3*Mr+A+S+P+2*D))*N + (5*N*Mr + (S+A+3*P)*N +A*(N-1) +P)*N +Mw + ...
        (Mr+D+Mw)*N + ((N+1)*Mr+(N-1)*A+D+Mw+Mm)*N + (2*Mr*N+P*N+A*(N-1)) + (3*Mr+S+A+2*P) ;

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
    
    f = (2+N)*Mm + Mw + Mm*(max(floor(lopt))+1) + Mm*N + Mm*N + 2*Mm*(N^2) + (6*Mr + 4*S + 3*P + Mw + Mm)*N + ...
                (3*Mr+S+OPP+5*P+2*D+E+R)*N + (5*N*Mr + (2*S+3*P)*N +S*(N-1) +P)*N +Mw + (Mr+D+Mw)*N + ((N+1)*Mr+(N-1)*S+D+Mw+Mm)*N + (2*Mr*N+P*N+S*(N-1)) + (3*Mr+2*S+2*P) ;
        
end


end


