function y = Plot_Path(NbPath, Nb_time_step, b_scheme, t_pillars, Param, Intensity_model)

x0 =  Param(1);

T1 = t_pillars(1);
T2 = t_pillars(2);
T3 = t_pillars(3);
T4 = t_pillars(4);

b1 = b_scheme(1);
b2 = b_scheme(2);
b3 = b_scheme(3);
b4 = b_scheme(4);

bp = 0.0001;

%% Simulation
%Nb_time_step = 1000;
Step = (T4-0)/Nb_time_step;
T = 0+(0:Nb_time_step)'.*Step;
%Integral_Approx = zeros(Nb_time_step+1, 1);

ColorSet = varycolor(NbPath);
set(gca, 'ColorOrder', ColorSet);
hold all;

switch Intensity_model
    case 1 %CIR case
                
        
        a = Param(2);
        sigma = Param(3);
        
        disp('Feller Condition : ');
        sigma_square = sigma^2
        ab_double = 2*a*min(b_scheme)
        Feller_Cond = sigma_square < ab_double
        
        for m = 1:NbPath

            X = zeros(Nb_time_step+1, 1);
            X(1) = x0;
            for i = 1:Nb_time_step  
               if (T(i)<T1) 
                   b = b1;
               elseif (T(i)<T2)
                   b = b2;
               elseif (T(i)<T3)
                   b = b3;
               else
                   b = b4;
               end
               
               X(i+1) = X(i) + a*(b - X(i))*Step + sigma*sqrt(X(i)*Step)*normrnd(0,1);
               %Integral_Approx(i+1) = Integral_Approx(i) + X(i+1)*Step;

            end

            plot(T, X/bp);
            title('Intensity Path AIG')
            xlabel('time-to-maturity (year)')
            ylabel('intensity (bp)')
        end        
        
        
    
    case 2  %Gamma OU

        %Param = [x; a; sigma; c; 1/Taille_Saut;]
        a = Param(2);
        sigma = Param(3);
        c = Param(4);
        lambda = Param(5);
        
        for m = 1:NbPath

            X = zeros(Nb_time_step+1, 1);
            X(1) = x0;
            for i = 1:Nb_time_step  
               if (T(i)<T1) 
                   b = b1;
               elseif (T(i)<T2)
                   b = b2;
               elseif (T(i)<T3)
                   b = b3;
               else
                   b = b4;
               end

               X(i+1) = X(i) + a*(b - X(i))*Step + sigma*gamrnd(c*Step,1/lambda);
               %Integral_Approx(i+1) = Integral_Approx(i) + X(i+1)*Step;

            end

            plot(T, X/bp);
            title('Intensity Path AIG')
            xlabel('time-to-maturity (year)')
            ylabel('intensity (bp)')
        end          
        
        
    case 3  % OU Inverse Gaussian 

        %Param = [x; a; sigma; c; 1/Taille_Saut;]
        a = Param(2);
        sigma = Param(3);
        c = Param(4);
        lambda = Param(5);
        
        for m = 1:NbPath

            X = zeros(Nb_time_step+1, 1);
            X(1) = x0;
            for i = 1:Nb_time_step  
               if (T(i)<T1) 
                   b = b1;
               elseif (T(i)<T2)
                   b = b2;
               elseif (T(i)<T3)
                   b = b3;
               else
                   b = b4;
               end

               X(i+1) = X(i) + a*(b - X(i))*Step + sigma*randraw('ig',[c*Step,lambda]);
               %X(i+1) = X(i) + a*(b - X(i))*Step + sigma*randraw('ig',[c*Step,1/lambda]);
               %Integral_Approx(i+1) = Integral_Approx(i) + X(i+1)*Step;

            end

            plot(T, X/bp);
            title('Intensity Path AIG')
            xlabel('time-to-maturity (year)')
            ylabel('intensity (bp)')
        end          
        
        
    otherwise
        exception = MException('VerifyInput:OutOfBounds', ...
       'Unknown model reference');
    throw(exception);     
end


