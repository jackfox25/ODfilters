% =========================================================================
% Function: ekf_OD
% Author:   Jack Fox
% Date:     02/05/25
%     General EKF used for OD. Implementing notation from Born et al.
% =========================================================================
% Inputs:
%          t_0   (1x1) Initial time (s)
%       Xhat_0   (nx1) Initial state guess
%          P_0   (nxn) Initial state covariance
%         meas   (kx(2+p)) True measurements with rows [t ID m1 m2 ... mp]
%       params         Struct containing problem-specific parameters
%     sysFuncs         Struct containing specific dynamic and measurement 
%                        functions. Any helper functions must be on path.
%       - sysFuncs.computeAugStateDot
%       - sysFuncs.computeMeasResidual
%    verbosity   (bool) Flag for output text
% Outputs:
%     ekfOut   Struct containing outputs from EKF
%      - Xhat_hist       (kxn) Full state history for each measurement time
%      - P_hist          (kxn^2) State covariance history for each measurement time
% =========================================================================
function ekfOut = ekf_OD(t_0,Xhat_0,P_0,meas,params,sysFuncs,verbosity)

    n = size(Xhat_0,1);

    % === PROCESS SET OF MEASUREMENTS AS CKF/LKF TO INITIALIZE === %
    nLKF = params.settings.EKFinitThr;
    if size(meas,1) > nLKF && nLKF > 0 && exist("lkf_OD.m","file")
        fprintf('EKF: Initializing %d measurements as LKF\n',nLKF);
        lkfMeas = meas(1:nLKF,:);
        lkfOut = lkf_OD(t_0,Xhat_0,P_0,lkfMeas,params,sysFuncs,verbosity);

        Xhat_hist_LKF = lkfOut.Xhat_hist;
        P_hist_LKF = lkfOut.P_hist;

        % Reset starting conditions for EKF
        t_0 = meas(nLKF,1);
        Xhat_0 = Xhat_hist_LKF(end,:)';
        P_0 = reshape(P_hist_LKF(end,:),n,n);
        meas = meas(nLKF+1:end,:);
    end
    % ============================================================ %

    k = size(meas,1);

    % Measurement noise covariance
    R = params.cov.R;

    % Initialize
    t_im1 = t_0;
    Xref_im1 = Xhat_0;
    P_im1 = P_0;

    STM_im1 = eye(n);

    % Create containers to save
    Xhat_hist = zeros(k,n);
    P_hist = zeros(k,n^2);

    % Loop over measurements
    for i=1:k

        % Read next observation
        t_i = meas(i,1);
        Y_i = meas(i,:);
        R_i = R;

        % If long enough has passed, run one measurement as LKF:
        if t_i - t_im1 > params.settings.EKFmissedMeasTime && exist("lkf_OD.m","file")
            fprintf('EKF: Measurement outage: re-initializing as LKF\n');
            lkfOut = lkf_OD(t_im1,Xref_im1,P_im1,Y_i,params,sysFuncs,verbosity);
            Xh = lkfOut.Xhat_hist;
            Ph = lkfOut.P_hist;

            % Save to history
            Xhat_hist(i,:) = Xh(end,:);
            P_hist(i,:) = reshape(Ph(end,:),n^2,1);
    
            % Set new im1 vars
            t_im1 = t_i;
            Xref_im1 = Xhat_hist(i,:)';
            P_im1 = reshape(P_hist(i,:),n,n);
    
            % Status
            if verbosity
                fprintf('LKF: Switching back to EKF\n');
            end

            continue;

        end


        % If this is a distinct time:
        if t_i > t_im1
            % Integrate reference trajectory and STM from t_im1 to t_i
            Xref_im1_aug = [Xref_im1; reshape(STM_im1,n^2,1)];
            [~,xah] = params.settings.integrator(@(t,x) sysFuncs.computeAugStateDot(t,x,params),...
                                 [t_im1 t_i],Xref_im1_aug,params.settings.OPTIONS);
            Xref_i = xah(end,1:n)';
            STM_im1toi = reshape(xah(end,n+1:end),n,n);
    
            if strcmp(params.settings.noiseComp,"SNC") && (t_i-t_im1) < params.SNC.thresh
                comp = params.SNC.Gam(t_i-t_im1) * getQ_SNC(Xref_im1,params) * params.SNC.Gam(t_i-t_im1)';

            elseif strcmp(params.settings.noiseComp,"DMC") && (t_i-t_im1) < params.DMC.thresh
                comp = getQ_DMC(t_i-t_im1,Xref_im1,params);
                
            else
                comp = zeros(n);
            end

            % Time update
            Pbar_i = STM_im1toi * P_im1 * STM_im1toi' + comp;
        else
            % At same time - do not perform time update
            Xref_i = Xref_im1;
            Pbar_i = P_im1;
        end

        % Compute observation deviation, observation-state matrix, gain matrix
        [y_i,Htilde_i] = sysFuncs.computeMeasResidual(t_i,Xref_i',Y_i,params);
        K_i = Pbar_i * Htilde_i' / (Htilde_i * Pbar_i * Htilde_i' + R_i);

        % Measurement and reference orbit update
        xhat_i = K_i * y_i;
        Xref_i = Xref_i + xhat_i;
        P_i = (eye(n) - K_i * Htilde_i) * Pbar_i * (eye(n) - K_i * Htilde_i)' + K_i * R * K_i';

        % Save to history
        Xhat_hist(i,:) = Xref_i';
        P_hist(i,:) = reshape(P_i,n^2,1)';

        % Set new im1 vars
        t_im1 = t_i;
        Xref_im1 = Xref_i;
        P_im1 = P_i;

        % Status
        if verbosity && (mod(i,10)==0 || i==k)
            fprintf('EKF: Processed measurement %d of %d\n',i,k);
        end

    end

    % Append EKF solution to LKF solution
    if size(meas,1) > nLKF && nLKF > 0
        Xhat_hist = [Xhat_hist_LKF; Xhat_hist];
        P_hist = [P_hist_LKF; P_hist];
    end

    ekfOut.Xhat_hist = Xhat_hist;
    ekfOut.P_hist = P_hist;

end