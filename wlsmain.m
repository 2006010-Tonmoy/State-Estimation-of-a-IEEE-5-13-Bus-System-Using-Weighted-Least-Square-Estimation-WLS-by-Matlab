% Power System State Estimation using Weighted Least Square Method
clc;
clear all;
num = input('Enter bus number : '); % IEEE - 5 here
p = input('Enter the confidence(%)for bad data : ');
ybus = ybusform(num); % Get YBus
zdata = zdatas(num); % Get Measurement data
bpq = bbusform(num); % Get B data
nbus = num; % Get number of buses
type = zdata(:,2); % Type of measurement, Vi = 1, Pi = 2, Qi = 3, Pij = 4, Qij = 5, Iij = 6
fbus = zdata(:,4); % From bus
tbus = zdata(:,5); % To bus
%% Meter Error Setting
err= (randn(1,length(fbus)))*1e-03; 
met_err = err - sum(err)/length(err);  %% making mean o
met_err = met_err'*met_err;            %% variance formtion
R = zeros(length(err),length(err));
for i = 1:length(err)
    for j = 1:length(err)
        if i==j
            R(i,j) = met_err(i,j);   
        end
    end
end
W= inv(R);  %% weightage matrix
%% Meter Measurement data
z = zdata(:,3)+(err)'; % Measuement values from meter

%% Calculation
%%... Initialize
V = ones(nbus,1);
del = zeros(nbus,1); 

E = [del(2:end); V];   % State Vector(except del1)
G = real(ybus);
B = imag(ybus);

vi = find(type == 1); % Index of voltage magnitude measurements..
rpi = find(type == 2); % Index of real power injection measurements..
qi = find(type == 3); % Index of reactive power injection measurements..
pf = find(type == 4); % Index of real powerflow measurements..
qf = find(type == 5); % Index of reactive powerflow measurements..

nvi = length(vi); % Number of Voltage measurements..
npi = length(rpi); % Number of Real Power Injection measurements..
nqi = length(qi); % Number of Reactive Power Injection measurements..
npf = length(pf); % Number of Real Power Flow measurements..
nqf = length(qf); % Number of Reactive Power Flow measurements..

iter = 1;
tol = 5;

stat_err = zeros(1,1000);
while(tol > 1e-08)
    
    %% Measurement Function, h
    h1 = zeros(nvi,1);
     for i = 1 : nvi
         k = fbus(i,1);
         h1(i,1) = V(k,1);
     end
    h2 = zeros(npi,1);  %% real power inj
    h3 = zeros(nqi,1);  %% react. power inj
    h4 = zeros(npf,1);  %% real power flow
    h5 = zeros(nqf,1);  %% react. power flow
    %% real power injection measurement
    for i = 1:npi
        m = fbus(rpi(i));  %% rpi =real power measurement er index(1,2,3...)
        for k = 1:nbus
            if k==m
                continue
            end    
            h2(i) = h2(i) + V(m)*V(k)*abs(ybus(m,k))*cos(angle(ybus(m,k))+del(k)-del(m)); %% real power
        end
           h2(i) = h2(i)+ G(m,m)*(V(m))^2;
    end
    
    for i = 1:nqi
        m = fbus(qi(i));
        for k = 1:nbus
            if k==m
                continue
            end
            h3(i) = h3(i) + V(m)*V(k)*abs(ybus(m,k))*sin(angle(ybus(m,k))+del(k)-del(m)); %% reactive power  
        end
         h3(i) = - (h3(i)+ B(m,m)*(V(m))^2);
    end
    
    for i = 1:npf              %% real power flow
        m = fbus(pf(i));
        n = tbus(pf(i));
        h4(i) = -V(m)^2*G(m,n) + V(m)*V(n)*abs(ybus(m,n))*cos(angle(ybus(m,n))+del(n)-del(m));
    end
    
    for i = 1:nqf        %% reactive power flow
        m = fbus(qf(i));
        n = tbus(qf(i));
        h5(i) = -(V(m)^2*(-B(m,n)+bpq(m,n)) + V(m)*V(n)*abs(ybus(m,n))*sin(angle(ybus(m,n))+del(n)-del(m)));
    end
    
    h = [h1; h2; h3; h4; h5]; %% self measurement
    
    %% value of e
    e = z - h;  
  
    
    %% Jacobian..
    % H11 - Derivative of V with respect to angles.. All Zeros
    H11 = zeros(nvi,nbus-1);

    % H12 - Derivative of V with respect to V.. 
    H12 = zeros(nvi,nbus);
    for k = 1:nvi
        for n = 1:nbus
            if n == k
                H12(k,n) = 1;
            end
        end
    end

    % H21 - Derivative of Real Power Injections with Angles..
    H21 = zeros(npi,nbus-1);
    for i = 1:npi         %% row
        m = fbus(rpi(i));

        for k = 1:(nbus-1)  %% column   [ del = k+1]
             if k+1 == m
                for n = 1:nbus   
                    if n==m
                        continue
                    end
                    H21(i,k) = H21(i,k) + V(m)* V(n)*abs(ybus(m,n))*sin(angle(ybus(m,n)-del(n)-del(m)));
                end
            else
                H21(i,k) = -V(m)* V(k+1)*abs(ybus(m,k+1))*sin(angle(ybus(m,k+1))+del(k+1)-del(m));
            end
        end
    end
    
    %H22 - Derivative of Real Power Injections with V..
    H22 = zeros(npi,nbus);
    for i = 1:npi
        m = fbus(rpi(i));
        for k = 1:(nbus)     %% for V
            if k == m
                for n = 1:nbus
                    if n==m
                        continue
                    end
                    H22(i,k) = H22(i,k) + V(n)*abs(ybus(m,n))*cos(angle(ybus(m,n))+del(n)-del(m));
                end
                H22(i,k) = H22(i,k) + 2*V(m)*G(m,m); 
            else
                H22(i,k) = V(m)*abs(ybus(m,k))*cos(angle(ybus(m,k))+del(k)-del(m));
            end
        end
    end
    
    % H31 - Derivative of Reactive Power Injections with Angles..
   H31 = zeros(nqi,nbus-1);
    for i = 1:nqi   %% nqi = no of rows
        m = fbus(qi(i));
        for k = 1:(nbus-1)  %% loop for del (no of columns)
            if k+1 == m
                for n = 1:nbus
                    if n==m
                        continue
                    end
                    H31(i,k) = H31(i,k) + V(m)* V(n)*abs(ybus(m,n))*cos(angle(ybus(m,n))+del(n)-del(m));
                end
            else
                H31(i,k) = -V(m)* V(k+1)*abs(ybus(m,k+1))*cos(angle(ybus(m,k+1))+del(k+1)-del(m));
            end
        end
    end
    
    % H32 - Derivative of Reactive Power Injections with V..
    H32 = zeros(nqi,nbus);
    for i = 1:nqi
        m = fbus(qi(i));
        for k = 1:(nbus)
            if k == m
                for n = 1:nbus
                    if n==m
                        continue
                    end
                    H32(i,k) = H32(i,k) + V(n)*abs(ybus(m,n))*sin(angle(ybus(m,n))+del(n)-del(m));
                end
                H32(i,k) = -H32(i,k) - 2*V(m)*B(m,m);
            else
                H32(i,k) = -V(m)*abs(ybus(m,k))*sin(angle(ybus(m,k))+del(k)-del(m));
            end
        end
    end
    
    % H41 - Derivative of Real Power Flows with Angles..
  H41 = zeros(npf,nbus-1);
    for i = 1:npf   
        m = fbus(pf(i));
        n = tbus(pf(i));
        for k = 1:(nbus-1)   %% for del
            if k+1 == m
                H41(i,k) = V(m)* V(n)*abs(ybus(m,n))*sin(angle(ybus(m,n)+del(n)-del(m)));
            else if k+1 == n
                H41(i,k) = -V(m)* V(n)*abs(ybus(m,n))*sin(angle(ybus(m,n)+del(n)-del(m)));
                else
                    H41(i,k) = 0;
                end
            end
        end
    end
    
    % H42 - Derivative of Real Power Flows with V..
  H42 = zeros(npf,nbus);
    for i = 1:npf
        m = fbus(pf(i));
        n = tbus(pf(i));
        for k = 1:nbus
            if k == m
                H42(i,k) = V(n)*abs(ybus(m,n))*cos(angle(ybus(m,n))+del(n)-del(m))-2*G(m,n)*V(m);
            else if k == n
                H42(i,k) = V(m)*abs(ybus(m,n))*cos(angle(ybus(m,n))+del(n)-del(m));
                else
                    H42(i,k) = 0;
                end
            end
        end
    end
    % H51 - Derivative of Reactive Power Flows with Angles..
       H51 = zeros(nqf,nbus-1);
    for i = 1:nqf
        m = fbus(qf(i));
        n = tbus(qf(i));
        for k = 1:(nbus-1)
            if k+1 == m
                H51(i,k) = V(m)* V(n)*abs(ybus(m,n))*cos(angle(ybus(m,n))+del(n)-del(m));
            else if k+1 == n
                H51(i,k) = -V(m)* V(n)*abs(ybus(m,n))*cos(angle(ybus(m,n))+del(n)-del(m));
                else
                    H51(i,k) = 0;
                end
            end
        end
    end
    
    % H52 - Derivative of Reactive Power Flows with V..
      H52 = zeros(nqf,nbus);
    for i = 1:nqf
        m = fbus(qf(i));
        n = tbus(qf(i));
        for k = 1:nbus
            if k == m
                H52(i,k) = -V(n)*abs(ybus(m,n))*sin(angle(ybus(m,n))+del(n)-del(m))-2*V(m)*(-B(m,n)+ bpq(m,n));
            else if k == n
                H52(i,k) = -V(m)*abs(ybus(m,n))*sin(angle(ybus(m,n))+del(n)-del(m));
                else
                    H52(i,k) = 0;
                end
            end
        end
    end
     % Measurement Jacobian, H..
    H = [H11 H12; H21 H22; H31 H32; H41 H42; H51 H52];
    
    % Gain Matrix, Gm..
    Gm = H'*inv(R)*H; 
    
    % State Vector..
    dE = inv(Gm)*(H'*inv(R)*e);  %% error in state var
    E = E + dE;
    del(2:end) = E(1:nbus-1);  
    V = E(nbus:end); 
     stat_err(iter) = tol;
    iter = iter + 1;
    tol = sum(abs(dE));
   
end

Del = 180/pi.*del;
E2 = [V;Del];
%% Estimated error

    h1 = zeros(nvi,1);
     for i = 1 : nvi
         k = fbus(i,1);
         h1(i,1) = V(k,1);
     end
    h2 = zeros(npi,1);  %% real power inj
    h3 = zeros(nqi,1);  %% react p i
    h4 = zeros(npf,1);  %% real p flow
    h5 = zeros(nqf,1);  %% reac p flow
    
    for i = 1:npi
        m = fbus(rpi(i));  %% rpi =real power measurement er index(1,2,3...)
        for k = 1:nbus
            if k==m
                continue
            end    
            h2(i) = h2(i) + V(m)*V(k)*abs(ybus(m,k))*cos(angle(ybus(m,k))+del(k)-del(m)); 
        end
           h2(i) = h2(i)+ G(m,m)*(V(m))^2;
    end
    
    for i = 1:nqi
        m = fbus(qi(i));
        for k = 1:nbus
            if k==m
                continue
            end
            h3(i) = h3(i) + V(m)*V(k)*abs(ybus(m,k))*sin(angle(ybus(m,k))+del(k)-del(m)); 
        end
         h3(i) = - (h3(i)+ B(m,m)*(V(m))^2);
    end
    
    for i = 1:npf              %% real power flow
        m = fbus(pf(i));
        n = tbus(pf(i));
        h4(i) = -V(m)^2*G(m,n) + V(m)*V(n)*abs(ybus(m,n))*cos(angle(ybus(m,n))+del(n)-del(m));
    end
    
    for i = 1:nqf        %% reactive power flow
        m = fbus(qf(i));
        n = tbus(qf(i));
        h5(i) = -(V(m)^2*(-B(m,n)+bpq(m,n)) + V(m)*V(n)*abs(ybus(m,n))*sin(angle(ybus(m,n))+del(n)-del(m)));
    end
    
    h_est = [h1; h2; h3; h4; h5];
    
    %% value of e
    r = z - h_est;
    f = sum(inv(R)*(r.^2));   %% objective function
    Nm = length(z); 
    Ns = length(E);
    dof = Nm-Ns;
    chi = chi2inv(p/100,dof);

    %% bad data identification   
        mat1 = (H*inv(Gm)*H'*W);
        CvE = zeros(length(mat1));
        for i = 1: length(CvE)
            for j = 1:length(CvE)
                if i == j
                    CvE(i,j) = mat1(i,j);  %% Covariance matrix
                end
            end
        end
       I   = eye(length(CvE));
       R_prime = (I-CvE)*inv(W);   %% R' matrix determination
       s_err = abs(inv(sqrt(R_prime))*r);  %% standardized error matrix
       max_st_err = max(s_err);
       bad_data = find(s_err == max_st_err);
       figure(1)
       plot(s_err);
       title('Plot of Standardized error');
       xlabel('Index of measurement data'),ylabel('Value');
       figure(2)
        plot(stat_err),xlim([1,iter]);
        title('Plot of error in state variable');
       xlabel('No of iteration'),ylabel('Value');
     fprintf('\n');
   
  disp('                           .....State Estimation of a power system Using WLS method....');
  disp('                                      ');
disp('                                     | Bus |  V(pu)  |  Angle(degree)  |     ');
for m = 1:nbus
    fprintf('                                     %4g', m); 
    fprintf('   %8.4f', V(m)); 
    fprintf('    %8.4f', Del(m)); 
    fprintf('\n');
end            
disp('                                  ---------------------------------------------');
fprintf('\n');
fprintf('        Comments about bad data : \n')
 if f > chi
        fprintf('        There is bad data\n');
    else
        fprintf('        There is no significant effect of bad data\n')
    end

     fprintf('         The %d no zdata has the max standardized error so it may be a bad data.\n',bad_data)

