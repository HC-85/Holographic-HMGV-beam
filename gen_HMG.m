function [X, Y, hmg] = gen_HMG(n, e, varargin)

%Parse inputs
ip = inputParser;
ip.CaseSensitive = true;
ip.KeepUnmatched = true;
addRequired(ip,'e',@isnumeric)
addRequired(ip,'n',@isnumeric)
addParameter(ip,'N',200,@isnumeric)
addParameter(ip,'L',2,@isnumeric)
addParameter(ip,'a',1,@isnumeric)
addParameter(ip,'w0',1,@isnumeric)
addParameter(ip,'lambda',632.8E-9,@isnumeric)
parse(ip, n, e, varargin{:})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REQUIRED
% Order (n)
% Eccentricity (e)

% OPTIONAL
% Domain resolution (N)
% Cartesian domain side length in beam waists (L)
% Major axis (a)
% Beam Waist (w0)
% Wavelength (lambda)

var = struct();
for ii = 1:length(ip.Parameters)
    par = ip.Parameters{ii};
    var.(par) = ip.Results.(par);
end

kt = 4*var.w0;           %Transverse component of k 
f = var.e*var.a;         %Semi-focal distance
q = (f*kt/2)^2;          %Elliptic parameter (f*kt/2)^2; 
f0 = 2*sqrt(q)/kt;       %f(z=0)    

x = linspace(-var.L*var.w0/2, var.L*var.w0/2, var.N); %Domain size
[X, Y] = meshgrid(x);         
UV = acosh((X+1i*Y)/f0);          %Elliptical transformation
U = real(UV);                     %Radial coordinates
V = imag(UV);                     %Angular coordinates

%Gaussian Beam
k = 2*pi/(var.lambda);   %Wavevector magnitude
z = 0;                   %Transverse plane
zr = k*(var.w0^2)/2;         %Rayleigh range
mu = 1 + 1i*z/zr;        %Parameter
gb = (exp(1i*k*z)/mu)*exp(-(X.^2+Y.^2)/(mu*var.w0^2));

hmg = cell(1,2);
for i = 1:2
    if mod(n, 2)==0
        even = 1;
        odd = 3;
    else
        even = 2;
        odd = 4;
    end

    %Characteristic values and expansion coefficients
    [~, mc_even, vt_even] = eig_Spm(even,q);  
    [~, mc_odd, vt_odd] = eig_Spm(odd,q); 

    %True orders
    t_even = find(vt_even==n);
    t_odd = find(vt_odd==n);

    %Mathieu functions
    R_even = arrayfun(@(u) Jpm_f(even,u,q,mc_even,t_even, t_even), U); %Even radial
    A_even = arrayfun(@(v) Spm_f(even,v,mc_even,t_even, t_even), V);   %Even angular
    
    R_odd = arrayfun(@(u) Jpm_f(odd,u,q,mc_odd,t_odd, t_odd), U);      %Odd radial
    A_odd = arrayfun(@(v) Spm_f(odd,v,mc_odd,t_odd, t_odd), V);        %Odd angular

    %Normalization constants
    N_even = Npm(even,mc_even,t_even);
    N_odd = Npm(odd,mc_odd,t_odd);
    
    % Mathieu-Gauss beam
    mg_even = N_even(t_even)*exp(-1i*(z*kt^2)/2*k*mu).*gb.*R_even.*A_even;     
    mg_even = mg_even/max(abs(mg_even),[],'all'); %Normalization
    mg_odd = N_odd(t_odd)*exp(-1i*(z*kt^2)/2*k*mu).*gb.*R_odd.*A_odd; 
    mg_odd = mg_odd/max(abs(mg_odd),[],'all'); %Normalization

    %Helical Mathieu-Gauss beam
    if i==1
        hmg{i} = mg_even + 1i.*mg_odd;  
    elseif i==2
        hmg{i} = mg_even - 1i.*mg_odd;
    end
end 
