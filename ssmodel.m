classdef ssmodel < handle
    %State space model
    %   an one dimensional state-space model.
    %   x_{t+1} = A x_t + B N(epsilon)
    %   y_t     = C x_t + D N(nu)
    
    properties
        x0 = 0.0;
        epsilon = 1.0;
        nu = 1.0;
        A = 1.0;
        B = 1.0;
        C = 1.0;
        D = 1.0;
    end
    
    properties (SetAccess = private)
        xt = 0.0;
        ep_cov_flag = false;
        nu_cov_flag = false;
    end
    
    methods
        function obj = ssmodel(x0, epsilon, nu, A, B, C, D)
            %SSMODEL Construct an instance of this class
            %   x0: the initial state.
            %   epsilon: covariance matrix for the state.
            %   nu: covariance matrix for the measurement.
            %   trans: four parameters, A, B, C, D
            if size(x0,2) ~= 1
                error('The input `x0` should be a column vector');
            end
            if size(epsilon,2) ~= 1 && ( size(epsilon,1) ~= size(epsilon,2) )
                error(['The input `epsilon` should be a column vector', ...
                    'or a square matrix']);
            end
            if size(nu,2) ~= 1 && ( size(nu,1) ~= size(nu,2) )
                error(['The input `epsilon` should be a column vector', ...
                    'or a square matrix']);
            end
            M_x = size(x0,1);
            M_e = size(epsilon,1);
            M_n = size(nu,1);
            if size(A,1) ~= size(A,2) || size(A,1) ~= M_x
                error(['The transmission matrix `A` should be a', ...
                '%d * %d matrix.'], M_x, M_x);
            end
            if size(B,1) ~= M_x || size(B,2) ~= M_e
                error(['The transmission matrix `B` should be a', ...
                '%d * %d matrix.'], M_x, M_e);
            end
            if size(C,2) ~= M_x
                error(['The transmission matrix `C` should be a', ...
                'N * %d matrix.'], M_x);
            end
            if size(D,1) ~= size(C,1) || size(D,2) ~= M_n
                error(['The transmission matrix `D` should be a', ...
                '%d * %d matrix.'], size(C,1), M_n);
            end
            
            obj.x0 = x0;
            
            obj.ep_cov_flag = size(epsilon,2) == 1;
            obj.epsilon = epsilon;
            obj.nu_cov_flag = size(nu,2) == 1;
            obj.nu = nu;
            
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            obj.xt = x0;
        end
        
        function get_obj = copy(obj)
            %get_xt get the current state
            %   return the current state directly.
            get_obj = ssmodel(obj.x0, obj.epsilon, obj.nu, ...
                obj.A, obj.B, obj.C, obj.D);
        end
        
        function xt = get_xt(obj)
            %get_xt get the current state
            %   return the current state directly.
            xt = obj.xt;
        end
        
        function yt = get_y(obj)
            %get_y get the measurement
            %   return the current measurement.
            if obj.ep_cov_flag
                rnd_nu = mvnrnd(zeros([size(obj.nu,1),1]), obj.nu, ...
                    [size(obj.nu,1),1]);
            else
                rnd_nu = obj.nu.*randn(size(obj.nu));
            end
            yt = obj.C * obj.xt + obj.D * rnd_nu;
        end
        
        function xtp = forward(obj)
            %forward forward the simulation for one step.
            %   return the current state directly.
            if obj.ep_cov_flag
                rnd_epsilon = mvnrnd(zeros([size(obj.epsilon,1),1]), ...
                    obj.epsilon, [size(obj.epsilon,1),1]);
            else
                rnd_epsilon = obj.epsilon.*randn(size(obj.epsilon));
            end
            xtp = obj.A * obj.xt + obj.B * rnd_epsilon;
            obj.xt = xtp;
        end
        
        function [xkpkp, skpkp] = estimate(obj, xkk, skk, ykp)
            %estimate estimate the mean and variance based on Kalman
            %filter.
            %   xkk, skk: the estimation for the current step.
            %   ykp: next measurement.
            %   return the estimation for the next step.
            if obj.ep_cov_flag
                Q = obj.epsilon.^2;
            else
                Q = diag(obj.epsilon.^2);
            end
            if obj.nu_cov_flag
                R = obj.nu.^2;
            else
                R = diag(obj.nu.^2);
            end
            xkpk = obj.A * xkk;
            skpk = obj.A * skk * obj.A' + Q;
            K = (skpk * obj.C') / (obj.C*skpk*obj.C + R);
            xkpkp = xkpk + K * (ykp - obj.C*xkpk);
            skpkp = skpk - K * (obj.C*skpk*obj.C + R) * K';
        end
    end
end

