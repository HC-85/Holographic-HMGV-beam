%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   eig_Spm FUNCTION      [p,m = e,o (even,odd)]
%   [va,mv,vt]=eig_Spm(KF,q)
%
%   INPUTS:      -q= scalar, elliptic parameter value (q > 0)
%               -KF= function code:  KF=1 even-even
%                                    KF=2 even-odd
%                                    KF=3 odd-even
%                                    KF=4 odd-odd
%   OUTPUTS:    -va=  vector of eigenvalues 'a'
%               -mv= matrix of expansion coefficients
%               -vt= values of order n:   
%                          KF=1   vt=[0 2 4 6 ... (2*ncoeffs-2)]
%                          KF=2   vt=[1 3 5 7 ... (2*ncoeffs-1)]
%                          KF=3   vt=[2 4 6 8 ... (2*ncoeffs)]
%                          KF=4   vt=[1 3 5 7 ... (2*ncoeffs-1)]
%   Default value for number of coefficients: ncoeffs=25 
%   'va' and 'mv' are determined for all 25 orders specified in 'vt'
%   (The size of 'va' is [ncoeffs 1], the size of 'vt' is [1 ncoeffs],
%    and the size of 'mv' is [ncoeffs ncoeffs])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   eig_Spm FUNCTION CALL
function [va,mv,vt]=eig_Spm(KF,q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-even---KF == 1----------------------------------------

if KF == 1  
%   COMPUTE THE MATRIX M 
    ik = 0:24;
    off_diag = q*ones(1,length(ik)-1); 
    off_diagd= off_diag; off_diagd(1)=2*off_diag(1);
    M= diag((2*ik).^2) + diag(off_diagd, -1) + diag(off_diag, 1);
%   COMPUTE EIGENVALUES AND EIGENVECTORS
    [u,d]   = eig(M,'nobalance');
    [va,num]= sort(diag(d));   % va is vector of eigenvalues 'a'
    v= u(:,num);               % matrix of eigenvectors
    vt=2*ik;                   % even values of order n 
%   PROCESS COEFFICIENTS    
%   (For a given order n, sum of coefficients = 1 )   
    nc=size(v,2);        % number of columns in v
    sum_v=sum(v);        % sum_vector over columns in v

%   Compute matrix mv of processed eigenvectors 
    mv=[];
    for k=1:nc
      cn=sum_v(k);      % constant of processing on column k
      ncv=v(:,k)/cn;    % processed column vector
      mv=[mv ncv];      % matrix mv of processed eigenvectors   
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-odd---KF == 2 ----------------------------------------

elseif KF == 2  
%   COMPUTE THE MATRIX M
    ik = 0:24;                         
    off_diag = q*ones(1,length(ik)-1);
    M= diag((2*ik+1).^2) + diag(off_diag, -1) + diag(off_diag, 1);
    M(1,1)= M(1,1) + q;
%   COMPUTE EIGENVALUES AND EIGENVECTORS
    [u,d]   = eig(M,'nobalance');
    [va,num]= sort(diag(d));  % va is vector of eigenvalues 'a'
    v= u(:,num);              % matrix of eigenvectors
    vt=2*ik+1;                % odd values of order n
    
%   PROCESS COEFFICIENTS 
%   (For a given order n, sum of coefficients = 1 )   
    nc=size(v,2);    % number of columns in v
    sum_v=sum(v);    % sum_vector over columns in v
    
%   Compute matrix mv of processed eigenvectors      
    mv=[];
    for k=1:nc
    cn=sum_v(k);     % constant of processing on column k 
    ncv=v(:,k)/cn;   % processed column vector
    mv=[mv ncv];     % matrix mv of processed eigenvectors  
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------odd-even--KF == 3 -----------------------------------------
elseif KF == 3  

%   COMPUTE THE MATRIX M    
    ik=1:25;                             
    off_diag = q*ones(1,length(ik)-1);     
    M= diag((2*ik).^2) + diag(off_diag, -1) + diag(off_diag, 1);

%   COMPUTE EIGENVALUES AND EIGENVECTORS
    [u,d]   = eig(M,'nobalance');
    [va,num]= sort(diag(d));   % va is vector of eigenvalues 'a'
    v=u(:,num);                % matrix of eigenvectors
    vt=2*ik;                   % even values of order n

%   PROCESS COEFFICIENTS 
%   [For a given order n, (dSpm/dv)=1 at v=0]   
    nc=size(v,2);              % number of columns in v
    
%   Compute matrix mv of processed eigenvectors  
    mv=[];
    for k=1:nc
    cv=v(:,k).*vt';   
    cn=sum(cv);           % constant of processing on column k
    ncv=v(:,k)/cn;        % processed column vector
    mv=[mv ncv];          % matrix mv of processed eigenvectors
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %---------------odd-odd--KF == 4 -----------------------------------------
 elseif KF == 4  

 %  COMPUTE THE MATRIX M    
    ik=0:24;                            
    off_diag =q*ones(1,length(ik)-1);
    M=diag((2*ik+1).^2) + diag(off_diag, -1) + diag(off_diag, 1);
    M(1,1)=M(1,1)-q;
%   COMPUTE EIGENVALUES AND EIGENVECTORS
    [u,d]   = eig(M,'nobalance');
    [va,num]=sort(diag(d));     % va is vector of eigenvalues 'a'
    v=u(:,num);                 % matrix of eigenvectors
    vt=2*ik+1;                  % odd values of order n
    
%   PROCESS COEFFICIENTS 
%   [For a given order n, (dSpm/dv)=1 at v=0]   
    nc=size(v,2);            % number of columns in v
    
%   Compute matrix mv of processed eigenvectors      
    mv=[];
    for k=1:nc
    cv=v(:,k).*vt';  
    cn=sum(cv);          % constant of processing on column k
    ncv=v(:,k)/cn;       % processed column vector
    mv=[mv ncv];         % matrix mv of processed eigenvectors
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 