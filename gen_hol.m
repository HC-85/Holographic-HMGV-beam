function gen_hol(profile_name, varargin)
%Parse inputs
ip = inputParser;
ip.CaseSensitive = true;
ip.KeepUnmatched = true;
addRequired(ip,'profile_name',@ischar)
addParameter(ip,'path','', @ischar)
addParameter(ip,'thz',5e-7,@isnumeric)
addParameter(ip,'thxy',5e-7,@isnumeric)
addParameter(ip,'lambda',632.8E-9,@isnumeric)
parse(ip, profile_name, varargin{:})

%Display parameters
disp(['Transverse profile: ', ip.Results.profile_name])
disp(['thz: ', num2str(ip.Results.thz)])
disp(['thxy: ', num2str(ip.Results.thxy)])

%Shortened struct
var = struct();
for ii = 1:length(ip.Parameters)
    par = ip.Parameters{ii};
    var.(par) = ip.Results.(par);
end

fn = fieldnames(ip.Unmatched);
for ii = 1:length(fn)
    par = fn{ii};
    var.(par) = ip.Unmatched.(par);
end

% Helical Mathieu-Gauss
if profile_name == 'hmg' 
    [X, Y, hmg] = gen_HMG(var.n, var.e, varargin{:});
    tr_profile = hmg{1};
end

k = 2*pi/(var.lambda);

for thz = var.thz
    for thxy = var.thxy
        kxy = k*sin(thz);
        kx = kxy*cos(thxy);
        ky = kxy*sin(thxy);

        V = exp(2i*(kx*X + ky*Y)).*tr_profile;

        A = angle(V);
        B = abs(V);
        B = asin(B/max(B(:)));
        T = .5+.5*sign(cos(A)+cos(B));

        filename = join([var.path, 'holo_', profile_name,'_e', num2str(var.e*10),'_n',num2str(var.n),'_thz',num2str(thz),'_thxy',num2str(thxy),'.bmp'], 'bmp');
        imwrite(1-T, filename)
    end
end