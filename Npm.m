%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NORMALIZATION FACTOR OF ANGULAR MATHIEU FUNCTIONS
%   y = Npm(KF,mv,nmax)      [p,m = e,o (even,odd)]
%
%   INPUTS:     -mv= matrix of expansion coefficients
%               -nmax= maximum order 
%               -KF= function code:  KF=1 even-even
%                                    KF=2 even-odd
%                                    KF=3 odd-even
%                                    KF=4 odd-odd      
     
%   OUTPUTS:    -y= vector of normalization factors for all 'nmax' orders 
%   'mv' is determined beforehand with function 'eig_Spm'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Npm FUNCTION CALL
function y = Npm(KF,mv,nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-even---KF == 1----------------------------------------

if KF == 1  

    ik=0:24;  vt=2*ik;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);  
    A0=Apm(1);
    
         yc = 2*pi*A0^2;
            for j = 2:ncoeffs
                yc = yc + pi * Apm(j) * Apm(j);
            end
            
    y=[y yc];
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-odd----KF == 2-----------------------------------------
%---------------odd-even ---KF == 3----------------------------------------
%---------------odd-odd ----KF == 4----------------------------------------


elseif KF ~= 1
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);      
    
    yc = pi * sum(Apm.*Apm);

    y=[y yc];
    
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end