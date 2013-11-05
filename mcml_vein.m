function [] = mcml_vein()
%MCML_VEIN Multi-ftn.layered Monte Carlo simulation, line source, vein present
%
%   Illumination width = 50 mm
%   (Illumination length =  mm) (not required information)
%   Detection width = 40 mm
%   Detection length = 40 mm
%   Bin resolution = 0.05 mm
%
%CONSTANTS (contained in params struct)
%   rusnum   number of scatters before Russian roulette engaged
%   rusfrac   fraction of occurances where photon survives roulette
%   kftn   number of groups of thousands of photons
%   g   anisotropy factor
%   g2   g^2
%   nrel   relative index of refraction (FROM tissue TO air)
%   musp   reduced scattering coefficient (mm^-1)
%   mua   absorption coefficient (mm^-1)
%   mutp   reduced total interaction coefficient (mm^-1)
%   mut   total interaction coefficient (mm^-1)
%   albedo   scattering albedo
%   zb   total depth to boundary (mm)
%   dv   depth to top of vein along central axis (mm)
%   rv   vein radius (mm)
%
%VARIABLES (contained in ftn struct)
%   x[3]   Cartesian coordinates
%   ct[3]   direction cosines
%   nrus   total number of scatters
%   wt   photon weight
%   layer   epidermis(1), dermis(2), subcutis(3), vein(4), "muscle"(5)
%
%OTHER
%   s   pathlength (mm)

params = struct; % Defines empty struct

% Seed the random number generator based on the current time
%stream = RandStream('mt19937ar','Seed',sum(100*clock));  %Needed for older matlabs
%RandStream.setDefaultStream(stream); %Needed for older matlabs
rng('shuffle'); % Works for new matlabs

% Read XLSX file containing list of MC sims to run & sim-specific data
[file, nruns] = read_list_sims;
disp(['Processing file ',file])

   
    params = read_param(1, file); % Fill up parameter struct
    
    lbin = zeros(1000,1); % Zero launch bins
    dbin = zeros(800, 800); % Zero detect bins
    
    %     disp(['Currently processing run ', int2str(currun), '/', int2str(nruns), ' ', int2str(1000*params.kftn), ' Photons'])
    tic %Start timer for performance measurements
    

        for ftncount = 1:(1000*params.kftn)
            
            local_lbin = zeros(1000,1);
            local_dbin = zeros(800, 800);
            
            ftn = struct;
            ftn = ftnini();
            
            % Score launch location
            if ftn.x(1) >= -25 && ftn.x(1) < 25 % Photon is within scoring range
                xbin = floor((ftn.x(1) + 25)/0.05) + 1;
                local_lbin(xbin,1) = ftn.wt;
            end
            
            % Main portion of code
            while (ftn.wt > 1e-16)
                [ftn, local_dbin] = tissue(ftn, params, local_dbin);
                
                if ftn.wt <= 1e-6 % Photon weight is too small
                    ftn.wt = 0;
                end
                
                if realsqrt(ftn.x(1)^2+ftn.x(2)^2) >= 50 % Photon is too far away (x & y dir only)
                    ftn.wt = 0;
                end
                
                if ftn.nrus >= params.rusnum % Photon has scattered too many times
                    if rand <= params.rusfrac
                        ftn.nrus = -10*ftn.nrus; % Reset number of scatters
                        ftn.wt = ftn.wt/params.rusfrac; % Correct weight for roulette
                    else
                        ftn.wt = 0;
                    end
                end
                
            end
            
            lbin = lbin + local_lbin;
            dbin = dbin + local_dbin;
        end
    
    toc %Stop timer for performance measurements, output time
    
    str = num2str(1);
    runloc = strcat('Sheet',str);
    xlswrite(file,lbin,runloc,'C20')
    xlswrite(file,dbin,runloc,'E20')
    

end

function [file, nruns] = read_list_sims()
%READ_LIST_SIMS Read in the location of number of simulations and sim-specific data
%   Detailed explanation goes here

[filename, pathname] = uigetfile('*.xlsx', 'Seleftn.ct the XLSX file containing the MC data to run');
file = strcat(pathname,filename);
nruns = xlsread(file,'Sheet0');

end

function [params] = read_param(currun, file)
%READ_PARAM Read in parameters for the current simulation
%   Detailed explanation goes here

str = num2str(currun);
runloc = strcat('Sheet',str);
[num_data, ~] = xlsread(file,runloc);

params.rusnum = num_data(1); % Russian roulette - # scatters between
params.rusfrac = num_data(2); % Russian roulette - survival fraftn.ction
params.kftn = num_data(3); % Thousands of photons in simuation
params.g(1:3) = num_data(4); % Anisotropy parameter
params.g(4) = num_data(15);
params.g2 = 1 - params.g.^2;
params.nrel = 1/num_data(5); % Relative index of refraftn.ction for moving FROM tissue INTO air
params.musp(1:3) = num_data(6);
params.musp(4) = num_data(14);
params.mua = [num_data(7) num_data(9) num_data(11) num_data(13)];
params.mutp = params.mua + params.musp;
params.mut = params.mua + params.musp./(1-params.g);
params.albedo = (params.mut - params.mua)./params.mut;
params.zb = [0 num_data(8) num_data(8)+num_data(10) num_data(8)+num_data(10)+num_data(12)];
params.dv = num_data(16);
params.rv = num_data(17);

end

function ftn = ftnini()
%FTNINI Initialize photon
%   Detailed explanation goes here
ftn.x(1) = (rand - 0.5)*50; % -25 mm <= x(1) <= 25 mm
ftn.x(2) = 0; % Line source along y = 0
ftn.x(3) = 0; % Launched at surface

ftn.ct = [0, 0, 1]; % Launch straight down (+z) into tissue

ftn.nrus = 0;
ftn.wt = 1;
ftn.layer = 1; % ftn.layer = 1 for epidermis

end

function [ftn, local_dbin] = tissue(ftn, params, local_dbin)
%TISSUE Map out next step in tissue
%   Detailed explanation goes here

s = -log(rand)/params.mut(ftn.layer); % Sample for new pathlength
[db, layer2] = distbound(ftn, params); % Determine distance to boundary along current direction vector
[dv, layer3] = distvein(ftn, params); % Determine distance to vein along current direction vector

% Is the step size greater than the distance to the boundary or vein?
if dv >= 0 && dv <= s && dv < db % The photon encounters the vein
    s = dv;
    ftn.layer = layer3;
    %     if ftn.layer == 4
    %         ftn.layer = 3;
    %     else
    %         ftn.layer = 4;
    %     end
    ftn.x = ftn.x + s*ftn.ct; % Move to new position on vein boundary
elseif db > 0 && db <= s && db < dv % The photon encounters the boundary
    s = db;
    if layer2 == 5 % layer = 4 is reserved for inside veing
        ftn.wt = 0;
    elseif layer2 == 0
        ftn.x = ftn.x + s*ftn.ct; % Move to new position at surface
        
        % Determine whether to reflect or refract
        bscwt = fresnel(params.nrel, ftn.ct);
        if rand <= bscwt % Reflects off surface back into tissue
            ftn.ct(3) = -ftn.ct(3);
            ftn.layer = 1;
        else % Successfully transmits and is scored/killed
            
            % Determine bin
            if (ftn.x(1) >= -20 && ftn.x(1) < 20) && (ftn.x(2) >= -20 && ftn.x(2) < 20) % Photon is within scoring range
                xbin = floor((ftn.x(1) + 20)/0.05) + 1;
                ybin = floor((ftn.x(2) + 20)/0.05) + 1;
                local_dbin(xbin,ybin) = ftn.wt; % Add weight to bin
            end
            
            ftn.wt = 0; % Kill photon
        end
    else
        ftn.layer = layer2;
        ftn.x = ftn.x + s*ftn.ct; % Move to new position at boundary of new layer
    end
else % The photon remains in the same layer
    ftn.x = ftn.x + s*ftn.ct; % Move to new position in same layer
    ftn.wt = params.albedo(ftn.layer)*ftn.wt; % Update photon weight
    ftn.nrus = ftn.nrus + 1;
    ftn.ct = scatter(ftn, params); % Sample for new direction vector
    % ftn.layer remains the same
end

end

function [db, layer2] = distbound(ftn, params)
%DISTBOUND Determine distance to boundary along current direction vector
%   Detailed explanation goes here

if ftn.layer == 4
    layer2 = ftn.layer;
    db = 1000;
else
    if ftn.ct(3) > 0
        layer2 = ftn.layer+1;
        z = params.zb(layer2);
        db = (z-ftn.x(3))/ftn.ct(3);
    elseif ftn.ct(3) < 0
        layer2 = ftn.layer-1;
        z = params.zb(ftn.layer);
        db = (z-ftn.x(3))/ftn.ct(3);
    else
        layer2 = ftn.layer;
        db = 1000; % Current photon path runs parallel to boundary. No intersection
    end
    
    if layer2 == 4 % Layer = 4 is reserved for vein
        layer2 = 5;
    end
end

end

function [dv, layer3] = distvein(ftn, params)
%DISTVEIN Determine distance to vein along current direction vector
%   Detailed explanation goes here

if ftn.layer == 3 || ftn.layer == 4 % Special check in layer 3 for interaction with vein
    a = ftn.ct(1)^2 + ftn.ct(3)^2; % Calculate discriminant
    b = 2*ftn.ct(1)*ftn.x(1) + 2*ftn.ct(3)*(ftn.x(3)-(params.dv + params.rv));
    c = ftn.x(1)^2 + (ftn.x(3)-(params.dv+params.rv))^2 - params.rv^2;
    D = b^2 - 4*a*c;
    if D > 0 % 2 possible intersection points (most common occurance)
        dva = (-b + sqrt(D))/(2*a);
        dvs = (-b - sqrt(D))/(2*a);
        
        if dva > 0 && dvs >= 0
            if ftn.layer == 3
                dv = dvs;
                layer3 = 4;
            else % ftn.layer == 4
                dv = dva;
                layer3 = 3;
            end
        elseif dva >= 0 && dvs < 0
            if ftn.layer == 4
                dv = dva;
                layer3 = 3;
            else
                dv = 1000;
                layer3 = 3;
            end
        elseif dva < 0 && dvs < 0
            if ftn.layer == 3
                dv = 1000;
                layer3 = 3;
            else
                dv = 0.0001;
                layer3 = 3;
            end
        else
            dv = 1000;
            layer3 = ftn.layer;
        end
        %         if dva > 0.001 && dvs > 0.001 % Select the smallest, positive distance
        %             dv = dvs;
        %         elseif dva < 0 && dvs < 0 % Both points are "behind" photon
        %             dv = 1000;
        %         else
        %             dv = dva; % Select only positive distance
        %         end
        %
        %         if dv < 0.001 % Check that dv isn't ridiculously small
        %             dv = 1000;
        %         end
        
    elseif D == 0 % 1 possible intersection point
        dv = 1000; % Treat like a glancing blow/no interaction
        layer3 = ftn.layer;
        %dv = -b/(2*a);
    else % 0 possible intersection points
        dv = 1000;
        layer3 = ftn.layer;
    end
    
else
    dv = 1000;
    layer3 = ftn.layer;
end

end

function bscwt = fresnel(nrel, ct)
%BSCWT Calculate the Fresnel refleftn.ction and transmission coefficients
%   Detailed explanation goes here
if nrel == 1
    bscwt = 0;
else
    nrsq = nrel^2;
    snisq = 1 - ct(3)^2;
    snrsq = snisq/nrsq;
    if snrsq >= 1
        bscwt = 1;
    else
        ns = realsqrt(nrsq - snisq);
        np = abs(nrsq*ct(3));
        rs = (abs(ct(3)) - ns)/(abs(ct(3)) + ns);
        rp = (ns - np)/(ns + np);
        bscwt = (rs^2 + rp^2)/2;
    end
end

end

function dir = scatter(ftn, params)
%SCATTER Sample for new scatter angle in tissue
%   Detailed explanation goes here

% Get sin & cos of azimuthal angle
phi = 2*pi*rand;
sinphi = sin(phi);
cosphi = cos(phi);

% Get scatter angle from Henyey-Greenstein
p = 2*params.g(ftn.layer)*rand + 1 - params.g(ftn.layer);
p = params.g2(ftn.layer)/p;
p = p^2;

costheta = (2 - params.g2(ftn.layer) - p)/(2*params.g(ftn.layer));
sintheta = realsqrt(1 - costheta^2);

% Change direction
if abs(ftn.ct(3)) > 0.999 % if theta = 0, then equations simplify (also prevents error if costheta > 1)
    ftn.ct(1) = sintheta*cosphi;
    ftn.ct(2) = sintheta*sinphi;
    ftn.ct(3) = costheta*ftn.ct(3)/abs(ftn.ct(3));
else
    sinnorm = realsqrt(1-ftn.ct(3)^2);
    cttemp1 = ftn.ct(1)*costheta + sintheta*(ftn.ct(1)*ftn.ct(3)*cosphi - ftn.ct(2)*sinphi)/sinnorm;
    cttemp2 = ftn.ct(2)*costheta + sintheta*(ftn.ct(2)*ftn.ct(3)*cosphi - ftn.ct(1)*sinphi)/sinnorm;
    
    ftn.ct(3) = ftn.ct(3)*costheta - sinnorm*sintheta*cosphi;
    ftn.ct(1) = cttemp1;
    ftn.ct(2) = cttemp2;
end

dir = ftn.ct/norm(ftn.ct);

end