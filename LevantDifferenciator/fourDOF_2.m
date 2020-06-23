%4DOF Jann
%x0 = [0 0  20 0.1 0.1 0.1 0.1];K_=[2 2];fourDOF(tspan, x0, K_);
%%
function fourDOF_2(tspan, x0, K)
global k1 k2
parameters

%Poles (control gains)
k1=K(1);k2=K(2);

%References
  xref = 20;       yref = 20;
 
  %Stiff ODEs solvers
  [t, X] = ode113(@eqMotion_2, tspan, [x0 0 0 0 0]);
  
  %Grafico todos los estados y las referencias
  figure;
  subplot(2,1,1); plot(t,X(:,1),'r--',t,t*0+xref,'LineWidth',2);
  subplot(2,1,2); plot(t,X(:,2),'r--',t,t*0+yref,'LineWidth',2);
  figure;
  subplot(2,1,1); plot(t,X(:,3));
  subplot(2,1,2); plot(t,X(:,4));
  figure;
  subplot(2,1,1); plot(t,X(:,5));
  subplot(2,1,2); plot(t,X(:,6));
end

%% 
function dX = eqMotion_2(t, X)
    global k1 k2
    global g Tphi C31 C32 C33 C21h C11h C23h C13h C24h C14h
    
    %states [x, y u, w, phi, psi]
    x = X(1);   y = X(2);   u = X(4);   w = X(5);   phi = X(6);     psi = X(7);     eta_u0=X(8);    eta_u1= X(9);      eta_w0=X(10);     eta_w1= X(11);
    
    %Diff Gains
    du_ref=eta_u1;      dw_ref=eta_w1;
    lambda_u0=20;   lambda_u1=10;   lambda_w0=20;   lambda_w1=10; 
    
    % du
    Va = sqrt(u^2 + w^2); % Air speed
    u_f = Va*C21h*w - Va.*C11h.*u - w./u.*tan(phi).*(g.*sin(phi)+w.*C33-w.*phi./Tphi);
    u_dR = Va.*( C23h.*w - C13h.*u - (w.^2)./u.*C31./Va.*tan(phi) );
    u_dL = Va.*( C24h.*w - C14h.*u + (w.^2)./u.*C32./Va.*tan(phi) );
    %du = u_f + u_dR*dR + u_dL*dL;

    % dw
    w_f = -Va.*C21h.*u - Va.*C11h.*w + tan(phi).*(g.*sin(phi)+w.*C33-w.*phi./Tphi) + g.*cos(phi);
    w_dR = -Va.*( C23h.*u + C13h.*w - w.*C31./Va.*tan(phi) );
    w_dL = -Va.*( C24h.*u + C14h.*w + w.*C32./Va.*tan(phi) );
    %dw = w_f + w_dR*dR + w_dL*dL;

    % dphi
    phi_f = C33 - phi./Tphi;
    phi_dR =  C31;
    phi_dL = -C32;
    %dphi =  phi_f + phi_dR*dR + phi_dL*dL;

    % dpsi
    psi_f = g./u.*tan(phi) + w./u.*C33./cos(phi) - w./u./cos(phi).*phi./Tphi;
    psi_dR =  w./u.*C31./cos(phi);
    psi_dL = -w./u.*C32./cos(phi);
    %dpsi = psi_f + psi_dR*dR + psi_dL*dL;

    %Referencias
    x_ref = 20;             y_ref = 20;
    dx_ref = 0;             dy_ref=0;
    
    % Block 1
    e1 = [x_ref;y_ref]-[x;y];
    M=[cos(psi),sin(psi)*sin(phi);sin(psi),-cos(psi)*sin(phi)];
    Ref=pinv(M)*(k1*e1+[dx_ref;dy_ref ]);
    u_ref =Ref(1);          w_ref = Ref(2);

    % Block 2
    e2 = [u_ref;w_ref] - [u;w];
    G=[u_dL,u_dR;w_dL,w_dR];
    U=pinv(G)*([du_ref;dw_ref]-[u_f;w_f]+k2*e2); 
    dL=U(1);  dR=U(2);
      
    % Perturbations (lambdas = L)  
    %L=[0;0;sin(t);0;0;0];
    L=[0;0;0;0;0;0;0;0;0;0;0];
      
    % ODEs
      dxLNED = u*cos(psi)+w*sin(psi)*sin(phi);
      dyLNED = u*sin(psi)-w*cos(psi)*sin(phi);
      du = u_f + u_dR*dR + u_dL*dL;
      dw = w_f + w_dR*dR + w_dL*dL;
      dpsi = psi_f + psi_dR*dR + psi_dL*dL;
      dphi = phi_f + phi_dR*dR + phi_dL*dL;
      dzLNED = w*cos(phi);
      %Diferenciadores Levant
      d_eta_u0=-lambda_u0*sqrt(abs(eta_u0-u_ref))*sign(eta_u0-u_ref)+eta_u1;
      d_eta_u1=-lambda_u1*sign(eta_u1-d_eta_u0);
      d_eta_w0=-lambda_w0*sqrt(abs(eta_w0-w_ref))*sign(eta_w0-w_ref)+eta_w1;
      d_eta_w1=-lambda_w1*sign(eta_w1-d_eta_w0);
            
      dX = [dxLNED; dyLNED;dzLNED; du; dw; dphi; dpsi;d_eta_u0;d_eta_u1;d_eta_w0;d_eta_w1] + L;
end