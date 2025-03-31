% =========================================================================
% Function: lkf_OD
% Author:   Jack Fox
% Date:     02/04/25
%     General LKF used for OD. Implementing notation from Born et al.
% =========================================================================
% Inputs:
%          t_0   (1x1) Initial time (s)
%       Xhat_0   (nx1) Initial state guess
%          P_0   (nxn) Initial state covariance
%         meas   (kx(2+p)) True measurements with rows [t ID m1 m2 ... mp]
%       params   Struct containing problem-specific parameters
%     sysFuncs   Struct containing specific dynamic and measurement 
%                  functions. Any helper functions must be on path.
%       - sysFuncs.computeAugStateDot
%       - sysFuncs.computeMeasResidual
%    verbosity   (bool) Flag for output text
% Outputs:
%    lkfOut   Struct containing outputs from LKF
%      - Xhat_hist       (kxn) Full state history for each measurement time
%      - P_hist          (kxn^2) State covariance history for each measurement time
%      - xhat_hist       (kxn) State deviation history for each measurement time
%      - Pbar_hist       (kxn^2) State covariance history prior to measurements
%      - STM_im1toi_hist (kxn^2) Incremental state transition matrices (ref)
%      - Q_hist          (kxn^2) Noise compensation history
% =========================================================================
function lkfOut = lkf_OD(t_0,Xhat_0,P_0,meas,params,sysFuncs,verbosity)

    n = size(Xhat_0,1);
    k = size(meas,1);

    % Measurement noise covariance
    R = params.cov.R;

    % Initialize
    t_im1 = t_0;
    xhat_im1 = zeros(n,1);
    P_im1 = P_0;
    
    Xref_im1 = Xhat_0;
    STM_im1 = eye(n);

    % Create containers to save
    Xref_hist = zeros(k,n);
    xhat_hist = zeros(k,n);
    P_hist = zeros(k,n^2);
    Pbar_hist = zeros(k,n^2);
    STM_im1toi_hist = zeros(k,n^2);
    Q_hist = zeros(k,n^2);

    % Loop over measurements
    for i=1:k

        % Read next observation
        t_i = meas(i,1);
        Y_i = meas(i,:);
        R_i = R;

        % If this is a distinct time:
        if t_i > t_im1
            % Integrate reference trajectory and STM from t_im1 to t_i
            Xref_im1_aug = [Xref_im1; reshape(STM_im1,n^2,1)];
            [~,xah] = params.settings.integrator(@(t,x) sysFuncs.computeAugStateDot(t,x,params),...
                                 [t_im1 t_i],Xref_im1_aug,params.settings.OPTIONS);
            Xref_i = xah(end,1:n)';
            STM_im1toi = reshape(xah(end,n+1:end),n,n);
    
            % Time update
            xbar_i = STM_im1toi * xhat_im1;

            if strcmp(params.settings.noiseComp,"SNC") && (t_i-t_im1) < params.SNC.thresh
                comp = params.SNC.Gam(t_i-t_im1) * getQ_SNC(Xref_im1,params) * params.SNC.Gam(t_i-t_im1)';
                
            elseif strcmp(params.settings.noiseComp,"DMC") && (t_i-t_im1) < params.DMC.thresh
                comp = getQ_DMC(t_i-t_im1,Xref_im1,params);
                
            else
                comp = zeros(n);
            end
            
            Pbar_i = STM_im1toi * P_im1 * STM_im1toi' + comp;
        else
            % At same time - do not perform time update
            Xref_i = Xref_im1;
            xbar_i = xhat_im1;
            Pbar_i = P_im1;
            STM_im1toi = eye(n);
            comp = zeros(n);
        end

        % Compute observation deviation, observation-state matrix, gain matrix
        [y_i,Htilde_i] = sysFuncs.computeMeasResidual(t_i,Xref_i',Y_i,params);

        K_i = Pbar_i * Htilde_i' / (Htilde_i * Pbar_i * Htilde_i' + R_i);

        % Measurement update
        xhat_i = xbar_i + K_i * (y_i - Htilde_i * xbar_i);
        P_i = (eye(n) - K_i * Htilde_i) * Pbar_i * (eye(n) - K_i * Htilde_i)' + K_i * R * K_i';

        % Save to history
        xhat_hist(i,:) = xhat_i';
        P_hist(i,:) = reshape(P_i,n^2,1)';
        Xref_hist(i,:) = Xref_i';
        Pbar_hist(i,:) = reshape(Pbar_i,n^2,1)';
        STM_im1toi_hist(i,:) = reshape(STM_im1toi,n^2,1)';
        Q_hist(i,:) = reshape(comp,n^2,1)';

        % Set new im1 vars
        t_im1 = t_i;
        xhat_im1 = xhat_i;
        P_im1 = P_i;
        
        Xref_im1 = Xref_i;

        % Status
        if verbosity && (mod(i,10)==0 || i==k)
            fprintf('LKF: Processed measurement %d of %d\n',i,k);
        end

    end

    Xhat_hist = Xref_hist + xhat_hist;

    % Package into output struct
    lkfOut.Xhat_hist = Xhat_hist;
    lkfOut.P_hist = P_hist;
    lkfOut.xhat_hist = xhat_hist;
    lkfOut.Pbar_hist = Pbar_hist;
    lkfOut.STM_im1toi_hist = STM_im1toi_hist;
    lkfOut.Q_hist = Q_hist;

end