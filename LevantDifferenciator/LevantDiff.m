%*************************************************
% Author: L.Arturo Torres-Romero
% based on 
% Levant A. (1998). Robust exact differentiation via sliding mode technique. 
% Automatica, 34(3), 379-384
%*************************************************
function LevantDiff()
    clear all;
    dt = 0.01;
    [t, X] = ode45(@LevantDiffOde,[0:dt:15],[0 0]);
    
    figure;
    plot(t,X(:,2),'--',t,cos(t),'-.',t,sin(t),'-'); title('Differential comparison');
    legend('Levant Diff', 'Diff = cos(t)', 'Original = sin(t)', 'Location', 'best');
end

%****************************************
% Levant differentiator
% ***************************************
function dX = LevantDiffOde(t,X) 
    z_0=X(1);    z_1= X(2); 
    
     %Diff Gains
     du_ref=z_1;   
     lambda_u0=20;   lambda_u1=10;  
     u_ref = sin(t);
     
    %Levant Differentiator
     d_z_0=-lambda_u0*sqrt(abs(z_0-u_ref))*sign(z_0-u_ref)+z_1;
     d_z_1=-lambda_u1*sign(z_1-d_z_0);
          
    dX = [d_z_0;d_z_1];
end


