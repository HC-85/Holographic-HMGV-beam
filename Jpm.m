%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   RADIAL MATHIEU FUNCTION OF THE FIRST KIND 
%   y = Jpm(KF,u,q,mv,nmax)       [p,m = e,o (even,odd)]    
%   
%   INPUTS:     -u= value of radial coordinate to compute function 
%               -q= elliptical parameter (q > 0)
%               -mv= matrix of expansion coefficients
%               -nmax= maximum order
%               -KF= function code:  KF=1 even-even
%                                    KF=2 even-odd
%                                    KF=3 odd-even
%                                    KF=4 odd-odd                                                        
%   OUTPUTS:    -y= vector of function values for all 'nmax' orders 
%   'mv'is determined beforehand with function 'eig_Spm'
%   This Radial Mathieu Function is approximated by an expansion of
%   products of Bessel functions. It is related to Spm(KF,i*u,mv,nmax) by
%   Spm(KF,i*u,mv,nmax)=sqrt(2*pi)*gpm(KF,q,mv,nmax)*Jpm(KF,u,q,mv,nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Jpm FUNCTION CALL
function y = Jpm(KF,u,q,mv,nmax)

    v1=sqrt(q)*exp(-u);   
    v2=sqrt(q)*exp(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-even---KF == 1----------------------------------------

if KF == 1  

    ik=0:24; vt=2*ik;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);  in=vt(k);
    A0=Apm(1);
    
         yc=A0*besselj(0,v1)*besselj(0,v2);
     
            for j=2:ncoeffs
    jm=fix(j-1);
    yc = yc + ((-1)^jm)*Apm(j)*besselj(jm,v1)*besselj(jm,v2);
           end

    r=fix(in/2);
    coef=((-1)^r)*sqrt(pi/2)/A0;
    yc=yc*coef;
    
    y=[y yc];

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-odd---KF == 2-----------------------------------------

elseif KF == 2
 
    ik=0:24;    vt=2*ik+1;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);  in=vt(k);
    A1=Apm(1);
    
         yc=A1*(besselj(1,v1)*besselj(0,v2)+besselj(1,v2)*besselj(0,v1));
     
            for j=2:ncoeffs
           
    jm=fix(j-1);
    yc = yc + ((-1)^jm)*Apm(j)*(besselj(jm,v1)*besselj(j,v2)+ ...
                      besselj(jm,v2)*besselj(j,v1));
           end
           
    r=fix((in-1)/2);
    coef=((-1)^r)*sqrt(pi/2)/A1;
    yc=yc*coef;
    
    y=[y yc];
    
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------odd-even--KF == 3------------------------------------------

elseif KF == 3
  
    ik=1:25;   vt=2*ik;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);  in=vt(k);
    A2=Apm(1);
    
        yc = A2*(besselj(0,v2)*besselj(2,v1)-besselj(2,v2)*besselj(0,v1));
     
            for j=2:ncoeffs
             
    jm=fix(j-1);
    jp=fix(j+1);

    yc = yc + ((-1)^j)*Apm(j)*(besselj(jp,v2)*besselj(jm,v1)- ...
                      besselj(jm,v2)*besselj(jp,v1));
           end

    r=fix(in/2);
    coef=((-1)^r)*sqrt(pi/2)/A2;
    yc=yc*coef;
    
    y=[y yc];
    
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------odd-odd--KF == 4-------------------------------------------

elseif KF == 4
  
    ik=0:24;   vt=2*ik+1;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);  in=vt(k);
    A1=Apm(1);
    
         yc = A1*(besselj(1,v2)*besselj(0,v1)-besselj(0,v2)*besselj(1,v1));
     
            for j=2:ncoeffs
             
    jm=fix(j-1);
    
    yc = yc + ((-1)^jm)*Apm(j)*(besselj(j,v2)*besselj(jm,v1)- ...
                      besselj(jm,v2)*besselj(j,v1));
           end

    r=fix((in-1)/2);
    coef=((-1)^r)*sqrt(pi/2)/A1;
    yc=yc*coef;
    
    y=[y yc];
    
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
