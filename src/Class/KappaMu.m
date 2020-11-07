classdef KappaMu
    %RICE This class holds all the parameters for the Rician fading model.
    % It calculates theoretical envelope and phase PDFs
    % It does a Monte Carlo simulation using the parameters
    
    properties(Constant, Hidden = true)
        NumSamples = 2E6; % number of samples
        r = 0:0.001:6 % envelope range for PDF ploteter
    end 
    
    properties(Access = public)
        kappa; % Rician K factor
        mu; % number of clusters
        r_hat; % root mean square of the signal
    end
    
    properties(Hidden = true) 
        multipathFading; % Found based on the inputs
        envelopeProbability; % Calculated theoretical envelope probability
        xdataEnv; % Simulated envelope density plot x values 
        ydataEnv; % Simlated envelope density plot y valyes
    end
    
    methods(Access = public)
        function obj = KappaMu(kappa,mu,r_hat)
            %ADDITIVESHADOWRICE Construct an instance of this class
            
            %   Assigning input values
            obj.kappa = input_Check(obj,kappa,'\kappa',0,50);
            obj.mu = input_Check(obj,mu,'\mu',1,10);
            obj.r_hat = input_Check(obj,r_hat,'\hat{r}^2',0.5,2.5);
            
            % other calculated properties
            obj.multipathFading = Multipath_Fading(obj);
            obj.envelopeProbability = envelope_PDF(obj);
            [obj.xdataEnv, obj.ydataEnv] = envelope_Density(obj);
        end
    end
    
    methods(Access = private)
        
        function data = input_Check(obj, data, name, lower, upper) 
            % intput_Check checks the user inputs and throws errors
            
            % checks if input is empty
            if isempty(data)
                error(strcat(name,' must be a numeric input'));
            end
            
            % inputs must be a number
            if ~isnumeric(data)
               error(strcat(name,' must be a number, not a %s.', class(data)));
            end
            
            % input must be within the range
            if data < lower || data > upper
               error(strcat(name,' must be in the range [',num2str(lower),', ',num2str(upper),'].'));
            end
            
            % mu must be integer
            if strcmp(name,'\mu') && mod(data, 1) ~= 0
                error(strcat(name,' must be an integer'));
            end
                
        end
        
        function [p_i, q_i] = means(obj)
            %means Calculates the means of the complex Gaussians 
            %representing the in-phase and quadrature components.

            d2 = (obj.r_hat^(2) * obj.kappa)/(1 + obj.kappa);
     
            p_i = sqrt(d2/(2.*obj.mu));
            q_i = p_i;
            
        end
        
        function [sigma2] = scattered_Component(obj)
            %scattered_Component Calculates the power of the scattered 
            %signal component.    
            
            sigma2 = obj.r_hat.^(2) ./(2 * obj.mu .* (1 + obj.kappa));
            
        end
        
        function [gaussians] = generate_Gaussians(obj, mean, sigma) 
            %generate_Gaussians Generates the Gaussian random variables 
            
            gaussians = normrnd(mean,sigma,[1,obj.NumSamples]);
        end
        
        function [multipathFading] = Multipath_Fading(obj) 
            %complex_MultipathFading Generates the Rician fading model 
            
            [p_i, q_i] = means(obj);
            [sigma2] = scattered_Component(obj);
            
            multipathFading = 0;
            for i = 1 : 1 : obj.mu
                X_i = generate_Gaussians(obj, p_i, sqrt(sigma2));
                Y_i = generate_Gaussians(obj, q_i, sqrt(sigma2));

                multipathFading = multipathFading + X_i.^(2) + Y_i.^(2);
            end 
            
        end    
        
        function [eProbTheor] = envelope_PDF(obj)
            %envelope_PDF Calculates the theoretical envelope PDF
              
            A = (2*obj.mu * ((1 + obj.kappa)^((obj.mu+1)/2)))...
                /((obj.kappa^((obj.mu - 1)/2)) * exp(obj.mu * obj.kappa));
            R = obj.r ./ obj.r_hat;
            B = (R.^obj.mu).*exp(- obj.mu.*(1 + obj.kappa).*(R.^ 2));
            C = besseli(obj.mu - 1, 2.*obj.mu.*R.*sqrt(obj.kappa.*(1 + obj.kappa)));
            eProbTheor = A .* B .* C ./ obj.r_hat;
            
        end
        
        function [xdataEnv, ydataEnv] = envelope_Density(obj)
            %envelope_Density Evaluates the envelope PDF
            R = sqrt(obj.multipathFading);

            [f,x] = ecdf(R);
            [ydataEnv, xdataEnv] = ecdfhist(f,x, 0:0.05:max(obj.r));
            
        end
            
    end
end

