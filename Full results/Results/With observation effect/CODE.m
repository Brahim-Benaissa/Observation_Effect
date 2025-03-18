clc;
clear;


%% Main Script

% Read galaxy rotation data
[galaxy_names, observed_radius, observed_speed, e_Vobs] = read_rotation_data('OBSERVED_Rotation.txt');

% Save the rotation data 
save('galaxy_Rotation_data.mat', 'galaxy_names', 'observed_radius', 'observed_speed');

% Read galaxy mass data
[gs_galaxy_names, gs_M2L, e_M2L, gs_M2Lb, gs_L36, e_L36, gs_Lbulge, ...
 gs_Rdisk, gs_Reff, Distance_to_G, e_D] = read_mass_data('OBSERVED_Mass.txt');

% Initialize storage with error containers
total_baryonic_mass = cell(size(galaxy_names));
total_baryonic_mass_err = cell(size(galaxy_names));
enclosed_masses = cell(size(galaxy_names));
enclosed_masses_err = cell(size(galaxy_names));
radius_uncertainties = cell(size(galaxy_names)); % New container for radius errors

% Initialize galaxy distance  
galaxy_distances = NaN(size(galaxy_names));

% Calcualte the masses for each galaxy in order
for g = 1:length(galaxy_names)
    current_galaxy = galaxy_names{g};
    
    % Find index of current galaxy in gs_galaxy_names
    idx = find(strcmp(gs_galaxy_names, current_galaxy), 1);
    
    % Get distance and its uncertainty
    D = Distance_to_G(idx); % Distance in Mpc
    e_D_val = e_D(idx); % Distance uncertainty
    
    % Convert radius uncertainties (δR/R = δD/D)
    R_kpc = observed_radius{g};
    e_R_kpc = R_kpc * (e_D_val/D); % Radius uncertainty in kpc
    radius_uncertainties{g} = e_R_kpc * 3.086e+19; % Convert to meters
    
       % Error propagation for disk luminosity
    L_disk = gs_L36(idx) - gs_Lbulge(idx);
    e_L_disk = sqrt(e_L36(idx)^2 + 0^2); % Assuming no Lbulge error
    
    % Disk mass with error propagation (M2L error only)
    M_disk = gs_M2L(idx) * L_disk * 1e9;
    e_M_disk = M_disk * (e_M2L(idx)/gs_M2L(idx));
    
    % Bulge mass with error propagation (M2Lb error not available)
    M_bulge = gs_M2Lb(idx) * gs_Lbulge(idx);
    e_M_bulge = 0; % No M2Lb uncertainty data
    
    % Total baryonic mass with error
    total_baryonic_mass{g} = M_disk + M_bulge;
    total_baryonic_mass_err{g} = sqrt(e_M_disk^2 + e_M_bulge^2);
    
    
    gs_Rdisk_ordered (g) = gs_Rdisk(idx);
    gs_Reff_ordered (g) = gs_Reff(idx);
    M_disk_ordered (g) = M_disk;
    M_bulge_ordered (g) = M_bulge;
        
        
    % Process enclosed masses with errors
    radii = observed_radius{g};
    enc_mass = zeros(size(radii));
    enc_mass_err = zeros(size(radii));
    
  

    % Loop over each radius to calculate enclosed mass
    for r_idx = 1:length(radii)
        r = radii(r_idx);

        % Enclosed disk mass (assuming exponential profile)
        [enc_disk, e_enc_disk] = compute_enclosed_disk_mass(r, M_disk, e_M_disk, gs_Rdisk(idx));

        % Enclosed bulge mass (assuming Hernquist profile)

        [enc_bulge, e_enc_bulge] = compute_enclosed_bulge_mass(r, M_bulge, e_M_bulge, gs_Reff(idx));

        % Total enclosed mass with error
        enc_mass(r_idx) = enc_disk + enc_bulge;
        enc_mass_err(r_idx) = sqrt(e_enc_disk^2 + e_enc_bulge^2);
    end

    % Store the calculated enclosed masses and galaxy distance
    enclosed_masses{g} = enc_mass;
    enclosed_masses_err{g} = enc_mass_err;
end

% Save results to a .mat file
save('galaxy_mass_properties.mat', 'galaxy_names', 'total_baryonic_mass', 'enclosed_masses', 'galaxy_distances');






for galaxy = 1:1:length(galaxy_names)
        
    plot_option = true;

    % Define algorithm parameters
    Max_Its = 50000;
    Pop = 20;
    EXP = 0.99;
    
    % Define scalar values for lower and upper bounds
    lb_scalar = 0;
    ub_scalar = 1;

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(['%%%% Galaxy: ' num2str(galaxy) ' %%%%%'])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')

    % Load data with uncertainties
    Robs_kpc = observed_radius{galaxy};
    Vobs_km_s = observed_speed{galaxy};
    e_Vobs_km_s = e_Vobs{galaxy};
       
    % Unit conversions
    Robs_m = Robs_kpc * 3.086e+19;
    
    Vobs_m_s = Vobs_km_s * 1e3;
    e_Vobs_m_s = e_Vobs_km_s * 1e3;
    
    % Get Newtonian predictions with errors
    current_enclosed = enclosed_masses{galaxy};
    current_enclosed_err = enclosed_masses_err{galaxy};
    solar_mass_kg = 1.9885e30;
    Enclosed_Mass_kg = current_enclosed * solar_mass_kg;
    e_Enclosed_Mass_kg = current_enclosed_err * solar_mass_kg;
    
    total_baryonic_mass = total_baryonic_mass{galaxy};
    gs_Rdisk = gs_Rdisk_ordered(galaxy);
    gs_Reff = gs_Reff_ordered(galaxy);
    M_bulge = M_bulge_ordered(galaxy);
    M_disk = M_disk_ordered (galaxy);
   
    % Get radius uncertainties
    e_Robs_m = radius_uncertainties{galaxy};
   
    % Replicate scalar bounds into vectors
    lb = repmat(lb_scalar, 1, length(Vobs_km_s)+1);
    ub = repmat(ub_scalar, 1, length(Vobs_km_s)+1);
   
    
    
    lb(1)=0;
    ub(1)=100;
    
    Num_Paramters = length(Vobs_km_s)+1;
    
    % Initialize plotting variables
    bestPositionsHistory = zeros(Max_Its, Num_Paramters);
    bestScoresHistory = zeros(1, Max_Its);
    
    % Optimization loop 
        tic;
        Eva = 0;
    
        % Initialize positions within bounds
        Pos = lb' + rand(Num_Paramters, Pop) .* (ub' - lb');
    
        Fit = zeros(1, Pop);
    
        % Initial population evaluation
        for i = 1:Pop
            [Norm_Error, VNewton, e_VNewton, Vfit, e_Vfit] = analyzeGalaxyOrbitalSpeed(Robs_m, Vobs_m_s, Enclosed_Mass_kg, Pos(:, i), e_Enclosed_Mass_kg, e_Robs_m);
            Fit(i) = Norm_Error;
            Eva = Eva + 1;
        end
    
        % Best solutions initialization
        Best_Poses = Pos;
        Best_Fits = Fit;
        [Optimum_Fit, ind] = min(Fit);
        Center = Pos(:, ind);
    
        % Initialization 
        Best_Poses = Pos;
        Best_Fits = Fit;
        MeanBest = mean(Best_Poses, 2);
        [Optimum_Fit, ind] = min(Fit);
        Center = Pos(:, ind);
        Dist_MeanBest = abs(Center - MeanBest);
    
        % Iterative search
        for It = 1:Max_Its
            % calculations for local search area
            Local_lb = Center - Dist_MeanBest;
            Local_ub = Center + Dist_MeanBest;
    
            % Generate new population within the local search area
            PosLoc = rand(Num_Paramters, Pop) .* (Local_ub - Local_lb) + Local_lb;
    
            % Update positions based on algorithm logic
            for i = 1:Pop
                if rand < EXP
                    Explore = (PosLoc(:, i) - Best_Poses(:, i));
                    Pos(:, i) = PosLoc(:, i) + Explore;
                else
                    Dist_Center = PosLoc(:, i) - Center;
                    Pos(:, i) = PosLoc(:, i) - rand * Dist_Center;
                end
            end
    
            for i = 1:Num_Paramters
                % Put back solutions that went outside the search field for the ith parameter
                Pos(i,:) = (Pos(i,:).*~(Pos(i,:)<lb(i))) + ((rand(1, Pop).*(ub(i)-lb(i)) + lb(i)).*(Pos(i,:)<lb(i))); 
                Pos(i,:) = (Pos(i,:).*~(Pos(i,:)>ub(i))) + ((rand(1, Pop).*(ub(i)-lb(i)) + lb(i)).*(Pos(i,:)>ub(i))); 
            end
    
            % Evaluate and update the best points
            for i = 1:Pop
                [Norm_Error, VNewton, e_VNewton, Vfit, e_Vfit] = analyzeGalaxyOrbitalSpeed( Robs_m, Vobs_m_s, Enclosed_Mass_kg, Pos(:, i), e_Enclosed_Mass_kg, e_Robs_m);
                Fit(i) = Norm_Error;
                Eva = Eva + 1;
                if Fit(i) < Best_Fits(i)
                    Best_Fits(i) = Fit(i);
                    Best_Poses(:, i) = Pos(:, i);
                end
            end
    
            % Update best solutions
            MeanBest = mean(Best_Poses, 2);
            [It_BestFit, ind] = min(Fit);
            It_BestPos = Pos(:, ind);
            Dist_MeanBest = Center - MeanBest;
    
            if It_BestFit < Optimum_Fit
                Optimum_Fit = It_BestFit;
                Center = It_BestPos;
            end
    
            % Display progress
            if mod(It, Pop) == 1
                disp(Optimum_Fit);
            end
    
            % Track best solution and its fitness
            bestPositionsHistory(It, :) = Center;
            bestScoresHistory(It) = Optimum_Fit;
    
        end
  

%   Store data
    All_VNewton{galaxy} = VNewton;
    All_e_VNewton{galaxy} = e_VNewton;
    All_Vobs{galaxy} = Vobs_m_s;
    All_e_Vobs{galaxy} = e_Vobs_m_s;
    All_Vfit{galaxy} = Vfit;
    All_e_Vfit{galaxy} = e_Vfit;
    All_Robs{galaxy} = Robs_m;
    All_e_Robs{galaxy} = e_Robs_m;
    All_Enclosed_Masses{galaxy} =  Enclosed_Mass_kg;
    All_e_Enclosed_Masses{galaxy} =  e_Enclosed_Mass_kg;
    

   %   Store Coefficients  
    All_Coef_value {galaxy} =  Center(:);
    Mass_Coef (galaxy) =  Center(1);
    Space_Scale_Coef {galaxy} =  1./Center(2:end);
    
        
    % Analysis 
    % (Space_Scale/R²)
    Space_Scale_over_R_squared {galaxy} =  Space_Scale_Coef{galaxy}./ All_Robs{galaxy}.^2';

    % Recci spacetime curvature
    [Recci_curve{galaxy},EscapeV_by_radius{galaxy}] = calculate_spacetime_curvature(All_Robs{galaxy}, All_Enclosed_Masses{galaxy});

    % correlation 
    corr_Space_Scale_r_sqrd_Curvature(galaxy) = corr(Space_Scale_over_R_squared{galaxy}, Recci_curve{galaxy}');
    

    
if plot_option  
    disp(['Stopped in: ' num2str(It) ' iterations']);
    fprintf('Fit: %f\n', Optimum_Fit);

  
    % Plotting with error bars
    fig = figure('Position', [100 100 800 600], 'Color', 'w');
    hold on;
    
    % Fit speed with error bars in a clear green
    errorbar(Robs_m, Vfit, e_Vfit, '.', 'MarkerSize', 20, ...
             'CapSize', 8, 'LineWidth', 1.8, 'DisplayName', 'Fit Speed', ...
             'Color', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor', [0.4660, 0.6740, 0.1880]);
   
    % Observed data with error bars in a distinct blue
    errorbar(Robs_m, Vobs_m_s, e_Vobs_m_s, '.', 'MarkerSize', 20, ...
             'CapSize', 8, 'LineWidth', 1.8, 'DisplayName', 'Observed Speed', ...
             'Color', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);

    % Newtonian prediction with error bars in a striking orange
    errorbar(Robs_m, VNewton, e_VNewton, '.', 'MarkerSize', 20, ...
             'CapSize', 8, 'LineWidth', 1.8, 'DisplayName', 'Newtonian(With Obs Efct)', ...
             'Color', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);


    xlabel('Observed Radius (m)', 'FontSize', 14);
    ylabel('Orbital Speed (m/s)', 'FontSize', 14);
    title(sprintf('Galaxy %d : %s', galaxy, galaxy_names{galaxy}), 'FontSize', 16);
    legend('Location', 'best', 'FontSize', 12);
    grid on;
    set(gca, 'GridColor', [0.8, 0.8, 0.8]);  % Set grid color to light gray
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    box on;
    
    % Save high-quality figure
    saveas(fig, sprintf('Galaxy_%d_%s_FullUncertainty.jpg', galaxy, galaxy_names{galaxy}), 'jpg');
 
    
   
     % Plotting Space_Scale vs Observed radius in the same style as Figure A
    fig = figure('Position', [100 100 800 600], 'Color', 'w');
    hold on;

    % Plotting Space_Scale Coefficient data with similar styling to error bars in Figure A
    plot(All_Robs{galaxy}, Space_Scale_Coef{galaxy}, '-o', 'MarkerSize', 5, ...
         'LineWidth', 1.8, 'DisplayName', 'Space_Scale Coefficient', ...
         'Color', [0.4940, 0.1840, 0.5560], 'MarkerFaceColor', [0.4940, 0.1840, 0.5560]);

    xlabel('Observed Radius (m)', 'FontSize', 14);
    ylabel('Space_Scale Scale Factor', 'FontSize', 14);
    title(sprintf('Galaxy %d : %s', galaxy, galaxy_names{galaxy}), 'FontSize', 16);
    grid on;
    set(gca, 'GridColor', [0.8, 0.8, 0.8]);  % Set grid color to light gray
    
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    box on;

    % Save high-quality figure
    saveas(fig, sprintf('Galaxy_%d_%s_Observed_R_vs_Space_Scale.jpg', galaxy, galaxy_names{galaxy}), 'jpg');
 

 
    close all;
end
     

    % Load the data from the MAT file
    load('galaxy_Rotation_data.mat');
    load('galaxy_mass_properties.mat');
end

 

% Save the data to a MAT file
save('Coef.mat', 'All_VNewton', 'All_Vobs', 'All_Vfit', 'All_Robs', 'All_Enclosed_Masses', 'All_Coef_value', 'Mass_Coef', 'Space_Scale_Coef', 'Space_Scale_over_R_squared', 'Recci_curve', 'EscapeV_by_radius', 'corr_Space_Scale_r_sqrd_Curvature');

save('Uncertainty.mat','All_e_VNewton','All_e_Vobs','All_e_Vfit','All_e_Robs','All_e_Enclosed_Masses');





%% objective function 
function [Norm_Error, VNewton, e_VNewton, Vfit, e_Vfit] = analyzeGalaxyOrbitalSpeed(Robs_m, Vobs_m_s, Enclosed_Mass_kg, X, e_Enclosed_Mass_kg, e_Robs_m)
    % Constants
    G = 6.67430e-11; % m³/kg/s²

    Enclosed_Mass_kg = Enclosed_Mass_kg * X(1); 
    
    % Calculate Newtonian velocity using ENCLOSED MASSES
    VNewton = sqrt(G * Enclosed_Mass_kg ./ Robs_m); % m/s
    
    % Error propagation formula: δV/V = 0.5√[(δM/M)² + (δR/R)²]
    relative_mass_error = e_Enclosed_Mass_kg ./ Enclosed_Mass_kg;
    relative_radius_error = e_Robs_m ./ Robs_m;
    combined_relative_error = 0.5 * sqrt(relative_mass_error.^2 + relative_radius_error.^2);
    
    % Absolute velocity uncertainty
    e_VNewton = VNewton .* combined_relative_error;
    
    % Apply scaling factors and calculate error
    Vfit = VNewton .* X(2:end)'; % m/s
    
    % Fit velocity uncertainty
    e_Vfit = Vfit .* combined_relative_error;
    
    
    Norm_Error = norm(Vfit - Vobs_m_s);
end





%% Helper Functions 
function [enc_disk, e_enc_disk] = compute_enclosed_disk_mass(r, M_disk, e_M_disk, R_disk)
    x = r / R_disk;
    enc_disk = M_disk * (1 - (1 + x) * exp(-x));
    e_enc_disk = e_M_disk * (1 - (1 + x) * exp(-x)); % Error scales linearly
end

function [enc_bulge, e_enc_bulge] = compute_enclosed_bulge_mass(r, M_bulge, e_M_bulge, R_eff)
    if R_eff > 0 && M_bulge > 0
        a = R_eff / 1.8153;
        enc_bulge = M_bulge * (r / (r + a))^2;
        e_enc_bulge = e_M_bulge * (r / (r + a))^2; % Error scales with mass
    else
        enc_bulge = 0;
        e_enc_bulge = 0;
    end
end




% Calculate the Ricci scalar curvature and escape velocity
function [R, v_esc] = calculate_spacetime_curvature(r, M)
    % Constants
    G = 6.67430e-11;  % gravitational constant (m^3 kg^-1 s^-2)
    c = 299792458;    % speed of light (m/s)
    
    % Compute derivative of M with respect to r using central differences
    dMdr = gradient(M, r);
    
    % Compute Ricci scalar curvature: R = (2G/c^2) * (dM/dr) / r^2
    R = 2 * G * dMdr ./ (c^2 * r.^2);
    
    % Compute escape velocity: v_esc = sqrt(2GM/r)
    v_esc = sqrt(2 * G * M ./ r);
    
end



% Function to read galaxy rotation data
function [galaxy_names, observed_radius, observed_speed, e_Vobs] = read_rotation_data(filename)
    % Open the file
    fileID = fopen(filename, 'r');
    
    % Initialize output variables
    galaxy_names = {};
    observed_radius = {};
    observed_speed = {};
    e_Vobs = {};
    
    % Read the file line by line
    while ~feof(fileID)
        line = fgetl(fileID);
        parts = strsplit(line);
        
        % Extract data: galaxy name, radius, and speed
        galaxy_name = char(parts(1));
        radius = str2double(parts{2});
        speed = str2double(parts{3});
        err_speed = str2double(parts{4});
        
        % Check if galaxy name exists in the list
        galaxy_index = find(strcmp(galaxy_names, galaxy_name), 1);
        
        % If galaxy name not found, add it to the list
        if isempty(galaxy_index)
            galaxy_names{end+1} = galaxy_name;
            observed_radius{end+1} = [];
            observed_speed{end+1} = [];
            e_Vobs{end+1} = [];
            galaxy_index = length(galaxy_names); % Get index of the new galaxy
        end
        
        % Append the radius and speed data
        observed_radius{galaxy_index} = [observed_radius{galaxy_index}, radius];
        observed_speed{galaxy_index} = [observed_speed{galaxy_index}, speed];
        e_Vobs{galaxy_index}(end+1) = err_speed;
        
    end
    
    % Close the file
    fclose(fileID);
end



% Function to read galaxy mass data
function [gs_galaxy_names, M2L, e_M2L, M2Lb, L36, e_L36, Lbulge,  ...
          Rdisk, Reff, Distance, e_D] = read_mass_data(filename)
    % Open the file
    fileID = fopen(filename, 'r');
    gs_galaxy_names = {};
    M2L = []; e_M2L = []; M2Lb = [];
    L36 = []; e_L36 = []; Lbulge = [];
    Rdisk = []; Reff = []; Distance = []; e_D = [];
    
    while ~feof(fileID)
        line = fgetl(fileID);
        parts = strsplit(strtrim(line));
        
        % Column indices based on SPARC data format
        gs_galaxy_names{end+1} = parts{1};
        M2L(end+1) = str2double(parts{2});
        e_M2L(end+1) = str2double(parts{3});
        M2Lb(end+1) = str2double(parts{4});
        Distance(end+1) = str2double(parts{6});
        e_D(end+1) = str2double(parts{7});
        L36(end+1) = str2double(parts{10});
        e_L36(end+1) = str2double(parts{11});
        Lbulge(end+1) = str2double(parts{12});
        Reff(end+1) = str2double(parts{13});
        Rdisk(end+1) = str2double(parts{17});
    end
    
    % Close the file
    fclose(fileID);
end





