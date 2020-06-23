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
% 
% 
function dX = LevantDiffOde(t,X) 
eta_u0=X(1);    eta_u1= X(2); 

 %Diff Gains
 du_ref=eta_u1;   
 %lambda_u0=20;   lambda_u1=10;  % original parameters of the experiments
  lambda_u0=200;   lambda_u1=2000;   % moving this parameters the precision can be increased

 %u_ref = sin(t);
 u_ref = sin((2*pi*t^2)/10);
%Diferenciadores Levant
 d_eta_u0=-lambda_u0*sqrt(abs(eta_u0-u_ref))*sign(eta_u0-u_ref)+eta_u1;
 d_eta_u1=-lambda_u1*sign(eta_u1-d_eta_u0);
      
dX = [d_eta_u0;d_eta_u1];
end


