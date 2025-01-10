%*************************************************
% Author: L.Arturo Torres-Romero
% based on 
% Levant A. (1998). Robust exact differentiation via sliding mode technique. 
% Automatica, 34(3), 379-384
%*************************************************
function LevantDiffComparison()
clear all;
dt = 0.01;
[t, X] = ode45(@LevantDiffOde,[0:dt:50],[0 0]);

figure;
diffFnc = sin((2*pi*t.^2)/10);
diffFnc = diff(diffFnc)/dt;
%plot(t,X(:,2),t,cos(t),t,sin(t)); title('Derivative comparation');
figure
plot(t,X(:,2),t(2:end),diffFnc); title('levant diff VS amtlab diff');
end

%*****************************************************
% Levant differentiator
% ***************************************************
function dX = LevantDiffOde(t,X) 
z_0=X(1);    z_1= X(2); 

 %Diff Gains
 du_ref=z_1;   
 %lambda_u0=20;   lambda_u1=10;  % original parameters of the experiments
  lambda_u0=200;   lambda_u1=2000;   % moving this parameters the precision can be increased

 %u_ref = sin(t);
 u_ref = sin((2*pi*t^2)/10);
%Diferenciadores Levant
 d_z_0=-lambda_u0*sqrt(abs(z_0-u_ref))*sign(z_0-u_ref)+z_1;
 d_z_1=-lambda_u1*sign(z_1-d_z_0);
      
dX = [d_z_0;d_z_1];
end


