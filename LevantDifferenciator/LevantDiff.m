function LevantDiff()
clear all;
dt = 0.01;
[t, X] = ode45(@LevantDiffOde,[0:dt:15],[0 0]);

figure;
%diffResult = diff(sin(t))/dt;
plot(t,X(:,2),t,cos(t),t,sin(t)); title('Derivative comparation');
%plot(t,X(:,2)); title('Derivative comparation');

end



%*****************************************************
% Levant differentiator
% 
% 
function dX = LevantDiffOde(t,X) 
eta_u0=X(1);    eta_u1= X(2); 

 %Diff Gains
 du_ref=eta_u1;   
 lambda_u0=20;   lambda_u1=10;  

 u_ref = sin(t);
 
%Diferenciadores Levant
 d_eta_u0=-lambda_u0*sqrt(abs(eta_u0-u_ref))*sign(eta_u0-u_ref)+eta_u1;
 d_eta_u1=-lambda_u1*sign(eta_u1-d_eta_u0);
      
dX = [d_eta_u0;d_eta_u1];
end


