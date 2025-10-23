function sigma_vm = vonMisesStress(tau_exhaust, tau_coolant, sigma_long)
% VONMISESSTRESS Compute von Mises equivalent stress (2D plane-stress or 3D)
%
% Usage:
%   sigma_vm = vonMisesStress(sxx, syy, sxy)                % 2D plane stress
%   sigma_vm = vonMisesStress(sxx, syy, szz, sxy, syz, sxz) % full 3D
%
% Optional: pass a final options struct to compute missing components from
% flow/pressure parameters. Example:
%   opts.v = 10; opts.D = 0.01; opts.rho = 1000; opts.mu = 1e-3; opts.eps_rough = 1e-5;
%   opts.p = 2e5; opts.r = 0.02; opts.t = 0.002; % pressure/radius/thickness
%   sigma_vm = vonMisesStress(NaN, syy, NaN, opts) % will compute sigma_x and tau_xy
%
% Inputs may be scalars, vectors or arrays of identical size (vectorized).
% Units: Pascals (Pa)

% Detect optional options struct at end
%
tau_coolant
n = nargin;
opts = [];
if n >= 1 && isstruct(varargin{end})
    opts = varargin{end};
    effectiveN = n - 1;
else
    effectiveN = n;
end

if effectiveN == 3
    sxx = varargin{1}; syy = varargin{2}; sxy = varargin{3};
    szz = 0; syz = 0; sxz = 0;
elseif effectiveN == 6
    sxx = varargin{1}; syy = varargin{2}; szz = varargin{3};
    sxy = varargin{4}; syz = varargin{5}; sxz = varargin{6};
else
    error('vonMisesStress:InvalidInput','Call with 3 (plane stress) or 6 (3D) components. Optionally pass an options struct as final argument.');
end

% If options provided, compute missing components (sxx or sxy) from flow/pressure
isMissing = @(x) isempty(x) || (isnumeric(x) && all(isnan(x(:))));
if ~isempty(opts)
    % compute shear from flow params if requested/possible
    % support two separate shear sources: regen and chamber (thrust)
    tau_regen = [];
    tau_thrust = [];

    % generic single-source fields (backwards compatible)
    if isfield(opts,'v') && isfield(opts,'D') && isfield(opts,'rho') && isfield(opts,'mu')
        epsr = 0; if isfield(opts,'eps_rough'), epsr = opts.eps_rough; end
        L = []; if isfield(opts,'L'), L = opts.L; end
        tau_generic = darcyWallShear(L, opts.D, opts.v, opts.rho, opts.mu, epsr);
        tau_regen = tau_generic;
    end

    % regen-specific fields
    if isfield(opts,'v_regen') && isfield(opts,'D_regen') && isfield(opts,'rho_regen') && isfield(opts,'mu_regen')
        epsr = 0; if isfield(opts,'eps_rough_regen'), epsr = opts.eps_rough_regen; end
        Lr = []; if isfield(opts,'L_regen'), Lr = opts.L_regen; end
        tau_regen = darcyWallShear(Lr, opts.D_regen, opts.v_regen, opts.rho_regen, opts.mu_regen, epsr);
    end

    % thrust/chamber-specific fields
    if isfield(opts,'v_thrust') && isfield(opts,'D_thrust') && isfield(opts,'rho_thrust') && isfield(opts,'mu_thrust')
        epst = 0; if isfield(opts,'eps_rough_thrust'), epst = opts.eps_rough_thrust; end
        Lt = []; if isfield(opts,'L_thrust'), Lt = opts.L_thrust; end
        tau_thrust = darcyWallShear(Lt, opts.D_thrust, opts.v_thrust, opts.rho_thrust, opts.mu_thrust, epst);
    end

    % compute longitudinal stress from pressure if requested
    if (isMissing(sxx) || (isfield(opts,'overrideSigmaLong') && opts.overrideSigmaLong)) && ...
            isfield(opts,'p') && isfield(opts,'r') && isfield(opts,'t')
        sxx = longitudinalStressFromPressure(opts.p, opts.r, opts.t);
    end

    % combine shear contributions into single shear term used in VM calculation
    % combination rule: root-sum-square of available shear components
    if isMissing(sxy)
        % if both regen and thrust present, combine; if only one present use that
        if ~isempty(tau_regen) && ~isempty(tau_thrust)
            sxy = sqrt(tau_regen.^2 + tau_thrust.^2);
        elseif ~isempty(tau_regen)
            sxy = tau_regen;
        elseif ~isempty(tau_thrust)
            sxy = tau_thrust;
        end
    else
        % if sxy provided and opts requests combination, combine them
        if isfield(opts,'combineShear') && opts.combineShear
            extrasq = 0;
            if ~isempty(tau_regen), extrasq = extrasq + tau_regen.^2; end
            if ~isempty(tau_thrust), extrasq = extrasq + tau_thrust.^2; end
            sxy = sqrt(sxy.^2 + extrasq);
        end
    end
end

% Basic size checks (vectorized)
if ~isequal(size(sxx), size(syy)) || ~isequal(size(sxx), size(szz)) || ...
   ~isequal(size(sxx), size(sxy)) || ~isequal(size(sxx), size(syz)) || ...
   ~isequal(size(sxx), size(sxz))
    error('vonMisesStress:SizeMismatch','All inputs must have the same size.');
end

% --- Computation: von Mises (engineering shear convention) ---
% 3D formula:
%   sigma_vm = sqrt( 0.5*((sxx-syy).^2 + (syy-szz).^2 + (szz-sxx).^2) + 3*(sxy.^2 + syz.^2 + sxz.^2) );
sigma_vm = sqrt( 0.5*((sxx - syy).^2 + (syy - szz).^2 + (szz - sxx).^2) ...
                 + 3*(sxy.^2 + syz.^2 + sxz.^2) );
end

%% Local helper: compute wall shear from Darcy-Weisbach
function tau_exhaust  = exhaustWallShear(chamberDiameterArray, heightStepNumber)

    eps_rough = 8e-7;

    vm_Data = readmatrix("CEAOutFz_10-22-25.xlsx");

    exhaustDensity = vm_Data(:,11);
    exhaustMach = vm_Data(:,4);
    exhaustSoS = vm_Data(:,12);
    exhaustArea = vm_Data(:,2);
    exhaustViscosity = vm_Data(:,7);

    chamberDiameter = chamberDiameterArray(heightStepNumber);

    exhaustVelocity = exhaustMach * exhaustSoS;

    LStar = 73.4776667586; % THIS IS NOT THE REAL NUMBER

    exhaustReynolds = (exhaustDensity *exhaustVelocity * LStar) / exhaustViscosity;

    exhaustFriction = 1.325/(log(((eps_rough)/(3.7*chamberDiameter))+(5.74)/((exhaustReynolds)^0.9)));

    tau_exhaust = (frictionFactorExhaust/4)*(1/2)*(exhaustDensity)*(exhaustVelocity^2);

end


%% Local helper: compute wall shear from Darcy-Weisbach
function tau_coolant = coolantWallShear(density, flowVelocity, flowPressure, frictionFactor)
    % coolant side wall shear stress

    tau_coolant = (frictionFactor/4)*(1/2)*(density)*(flowVelocity^2);

end

function sigma_long = longitudinalStressFromPressure(flowPressure, chamberDiameterArray)

    vm_Data = readmatrix("CEAOutFz_10-22-25.xlsx");

     exhaustPressure = (vm_Data(:,9))/100000; % Pa

     dP = flowPressure - exhaustPressure;

     chamberDiameter = chamberDiameterArray(heightStepNumber); %?? unit?

     wallThick = 0.0015556992 % m

      % sigma_long = longitudinalStressFromPressure(p, r, t)
    % Thin-wall longitudinal stress: sigma_long = p * r / (2*t)
    % Inputs:
    %   p - internal pressure (Pa)
    %   r - inner radius (m)
    %   t - wall thickness (m)

    sigma_long = (-dP * (chamberDiameter/2)) / (2 * wallThick);
end
