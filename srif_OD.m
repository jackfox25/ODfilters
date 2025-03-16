% =========================================================================
% Function: srif_OD
% Author:   Jack Fox
% Date:     03/10/25
%     General SRIF used for OD. Implementing notation from Born et al.
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
%    srifOut   Struct containing outputs from SRIF
%      - Xhat_hist       (kxn) Full state history for each measurement time
%      - P_hist          (kxn^2) State covariance history for each measurement time
%      - xhat_hist       (kxn) State deviation history for each measurement time
%      - Rbar_hist       (kxn^2) Square root of state covariance history
%      - STM_im1toi_hist (kxn^2) Incremental state transition matrices (ref)
% =========================================================================
function srifOut = srif_OD(t_0,Xhat_0,P_0,meas,params,sysFuncs,verbosity)
    
    n = size(Xhat_0,1);
    k = size(meas,1);

    xhat_0 = zeros(n,1);
    xbar_im1 = xhat_0;

    % Initialize Information Form
    Rbar_im1 = chol(inv(P_0), 'upper');  % Upper-triangular factor of Information matrix
    bbar_im1 = Rbar_im1 * xhat_0;               % Information vector

    if strcmp(params.settings.noiseComp,"SNC")
        U = chol(params.SNC.Q);
        Ru_inv = U';
        Ru = inv(Ru_inv);
        ubar = zeros(3,1);
    end

    % Initialize reference trajectory and state transition matrix
    t_im1 = t_0;
    Xref_im1 = Xhat_0;  % Initial reference state
    STM_im1 = eye(n);   % Initial state transition matrix

    % Create containers for output
    Xhat_hist = zeros(k,n);
    P_hist = zeros(k,n^2);
    xhat_hist = zeros(k,n);
    Rbar_hist = zeros(k,n^2);
    Xref_hist = zeros(k,n);
    STM_im1toi_hist = zeros(k,n^2);

    for i=1:k
        % Read measurement
        t_i = meas(i,1);
        Y_i = meas(i,:);

        % If this is a distinct time step, propagate reference state and STM
        if t_i > t_im1
            Xref_im1_aug = [Xref_im1; reshape(STM_im1,n^2,1)];
            
            % Propagate reference trajectory and STM using dynamics
            [~, xah] = ode45(@(t,x) sysFuncs.computeAugStateDot(t,x,params),...
                             [t_im1 t_i], Xref_im1_aug, params.settings.OPTIONS);
            
            Xref_i = xah(end,1:n)';
            STM_im1toi = reshape(xah(end,n+1:end),n,n);
                
            % Process noise (SNC):
            if strcmp(params.settings.noiseComp,"SNC") && (t_i-t_im1) < params.SNC.thresh
                Rtilde_i = Rbar_im1 / STM_im1toi;
                A = [Ru zeros(3,n) Ru_inv\ubar;...
                     -Rtilde_i*params.SNC.Gam(t_i-t_im1) Rtilde_i bbar_im1];
                [~,AHH] = qr(A);

                Rbar_i = AHH(4:end,4:end-1);
                bbar_i = AHH(4:end,end);
                xbar_i = Rbar_i \ bbar_i;
            
            % No process noise:
            else
                Rbar_i = Rbar_im1 / STM_im1toi;
                 % Force Rbar to be upper triangular?
                if params.settings.RHH
                    [~,Rbar_i] = qr(Rbar_i);
                end
                xbar_im1 = Rbar_im1 \ bbar_im1;
                xbar_i = STM_im1toi * xbar_im1;
                bbar_i = Rbar_i * xbar_i;
                
            end
        
        else
            % Same time - nothing changes
            Xref_i = Xref_im1;
            STM_im1toi = eye(n);

            xbar_i = xbar_im1;
            Rbar_i = Rbar_im1;
            bbar_i = bbar_im1;
        end

        % Compute pre-whitened measurement residual and Jacobian
        [y_i, Htilde_i] = sysFuncs.computeMeasResidual(t_i,Xref_i',Y_i,params);
        V = chol(params.cov.R,'upper'); % Pre-whitening matrix

        % Apply pre-whitening to measurement
        y_i_W = V \ y_i;
        Htilde_i_W = V \ Htilde_i;

        % Construct augmented matrix for Householder transformation
        A = [Rbar_i, bbar_i; Htilde_i_W, y_i_W];

        % Apply Householder QR decomposition (Section 5.5.2)
        [~, R_aug] = qr(A);

        % Extract updated R and b (ensuring upper triangular form)
        Rbar_i = R_aug(1:n,1:n);
        bbar_i = R_aug(1:n,end);

        % Compute state estimate deviation
        xhat_i = Rbar_i \ bbar_i;

        % Store results
        Xhat_hist(i,:) = (Xref_i + xhat_i)';  % Convert back to full state
        P_hist(i,:) = reshape(inv(Rbar_i)*inv(Rbar_i'),n^2,1)';
        xhat_hist(i,:) = xhat_i';
        Rbar_hist(i,:) = reshape(Rbar_i,n^2,1)';
        Xref_hist(i,:) = Xref_i';
        STM_im1toi_hist(i,:) = reshape(STM_im1toi,n^2,1)';

        % Update im1 vars
        t_im1 = t_i;
        Xref_im1 = Xref_i;

        xbar_im1 = xbar_i;
        Rbar_im1 = Rbar_i;
        bbar_im1 = bbar_i;

        % Display status
        if verbosity && mod(i,10)==0
            fprintf('SRIF: Processed measurement %d of %d\n',i,k);
        end
    end

    % Store outputs in the same format as LKF
    srifOut.Xhat_hist = Xhat_hist;
    srifOut.P_hist = P_hist;
    srifOut.xhat_hist = xhat_hist;
    srifOut.Rbar_hist = Rbar_hist;
    srifOut.Xref_hist = Xref_hist;
    srifOut.STM_im1toi_hist = STM_im1toi_hist;

end
