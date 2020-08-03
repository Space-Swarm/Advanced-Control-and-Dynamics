%Francesco Labia
%2 August 2020
%Program to analyse the control system of a second order system using state
%feedback and observer design

%Numerator and denominator of transfer function
num = [1];
den = [1 0.6 4];

%shift value in transfer function
shift = 0;

%Convert to state space representation
[A,B,C,D] = controllable(num,den,shift)
[num1, den2, shift] = controllable_to_tf(A,B,C,D)
%step plot and analysis
figure(1);
sys = ss(A,B,C,D);
step(sys);
S = stepinfo(sys);

%checking if completely controllable and observable
state_control = s_control(A,B);
output_control = out_control(A,B,C,1);
observor_control = obs_control(A,C);

%initialising required symbolic notation
syms s; 
syms f1;
syms f2;
    
damping = 0.5; %damping ratio
frequency = 5; %undamped natural frequency
char_eq = [1 2*damping*frequency frequency^2] %defining the characteristic equation coefficients
delta = char_eq(1)*s^2+char_eq(2)*s+char_eq(3); %defining the characteristic equation using symbolic notation
[F F_bar] = Fs(delta,A,B) %determining matrice F by calling function Fs 

%Calculating H and F_bar
H = zgain(B,C,D,F,F_bar)
F_bar = [0 1; -4 -2]

%Closed loop matrices
Ac_temp = A + B.*F
Ac = zeros(size(Ac_temp))
Ac(1,1) = Ac_temp(1,1), Ac(1,2) = Ac_temp(1,2),Ac(2,1) = Ac_temp(2,1),Ac(2,2) = Ac_temp(2,2);
Ac
Bc_temp = B.*H
Bc = zeros(size(Bc_temp));
Bc(1,1) = Bc_temp(1,1), Bc(2,1) = Bc_temp(2,1)
Cc_temp = C+D*F
Cc = zeros(size(Cc_temp))
Cc(1,1) = Cc_temp(1,1), Cc(1,2) = Cc_temp(1,2)
Dc = D*H
Dc = zeros(1,1)
Dc(1) = Dc(1)
G = ((C+D*F)*inv(s*eye(2,2)-F_bar)*B+D)*H

%plotting unit step response for the closed loop control
figure(2);
step(Ac,Bc,Cc,Dc);
sys = ss(Ac,Bc,Cc,Dc);
%plotting unit step response of open loop control on same plot to compare
%the two
hold on;
step(A,B,C,D);

stepinfo(sys)

[L,M,N] = observer(A,B,C,D,delta)
function [A B C D] = controllable(num,den,shift)
    %removing any leading zeros from the numerator
    index = find(num ~= 0, 1, 'first');
    num = num(index:end);
    
    %finding dimensions of numerator and denominator
    num_dim = size(num,2);
    den_dim = size(den,2)-1;%-1 is used because the first term is ignored
    
    %determing the state space variables in controllable form: A,B,C,D
    A = [zeros(den_dim-1,1) eye(den_dim-1,den_dim-1); -flip(den(2:end))];
   
    B = [zeros(den_dim-1,1); 1];
    
    if(num_dim>1)
        C = [flip(num)];
    else
        C = [num zeros(1,den_dim-1)];
    end
  
    D = shift;
end
%Function to convert from state space controllable representation, to
%transfer function representation
function [den,num,shift] = controllable_to_tf(A,B,C,D)
    num = [1 -flip(A(rank(A),:))];
    dim = size(C);
    if C(2) == 0 && dim(2) == 2
      den = flip(C(1));
    else
      den = flip(C);
    end
    shift = D;
end

%Function to convert a transfer function to the observable form in the
%state space
function [A B C D] = observable(num,den,shift)
    
    %removing any leading zeros from the num
    index = find(num ~= 0, 1, 'first');
    num = num(index:end)
    
    %finding dimensions of numerator and denominator
    num_dim = size(num,2);
    den_dim = size(den,2)-1; %-1 is used because the first term is ignored
    
    %determing the state space variables in observable form: A,B,C,D
    A = [zeros(1,den_dim-1) -den(end); eye(den_dim-1,den_dim-1) -flip(den(2:end-1)')];
    if(num_dim>1)
        B = [flip(num)]';
    else
        B = [num zeros(1,den_dim-1)]';
    end

    C = [zeros(1,den_dim-1) 1];
    D = shift;
end

%Function to check if the state is controllable
function state_control = s_control(A,B)
    P = [B A*B];
    if(rank(P) == rank(A))
        state_control = true;
    else
        state_control = false;
    end
    
end

%Function to check if the output is controllable
function output_control = out_control(A,B,C,l)
    Q = [C*B C*A*B];
    if(rank(Q) == l)
        output_control = true;
    else
        output_control = false;
    end
end

%function to check if the system is observable
function observer_control = obs_control(A,C)
    R = [C ;C*A];
    if(rank(R) == rank(A))
        observer_control = true;
    else
        observer_control = false;
    end
end

%Function to calculate F matrix
%Inputs: delta (symbolic equation), matrices A and B
%Output: marix F
function [F F_bar] = Fs(delta,A,B)
    %initialising required symbolic notation
    syms s;
    syms f1;
    syms f2;
    
    F_temp = [f1 f2]; %initialising F in a symbolic matrix
    Fbar = A + B*F_temp %calculating fbar
    
    delta_F = det(s*eye(size(Fbar))-Fbar) %working out characteric equation with sym f1 and f2
    
    coeff_delta_F = flip(coeffs(delta_F,s)) %finding coefficients of the characteristic equation with sym f1 and f2
    coeff_delta = flip(coeffs(delta,s)) %finding coefficients of the characteristic equation
    
    temp_sym = [symvar(coeff_delta_F(2)) symvar(coeff_delta_F(3))]; %setting temporary matrix to see which sym variable is stored in which coefficient
    
    %based on which coefficient has f1 and f2, equating the two characteristic equations and solving for f1 and f2
    if(f2 == temp_sym(1))
        f2_value = solve(coeff_delta_F(2) == coeff_delta(2));
        f1_value = solve(coeff_delta_F(3) == coeff_delta(3));
        F = [f1_value f2_value]; %returning F
    else
        f1_value = solve(coeff_delta_F(2) == coeff_delta(2));
        f2 = solve(coeff_delta_F(3) == coeff_delta(3));
        F = [f1_value f2_value]; %returning F
    end
    F_bar = subs(Fbar,[f1 f2],[f1_value f2_value]); %subbing in the values for F_bar for use to find H
end

%Function to work out the input feedforward gain matrix, H
%Inputs: B, C, D, F matrices
%Outputs: Input feedforward gain matrix, H
function H = zgain(B,C,D,F,F_bar)
    H = -inv((C+D*F)*inv(F_bar)*B-D);
end

%Function to work out matrices L, M and N for an observer controller design
function [L,M,N] = observer(A,B,C,D,delta)
    %initialising required symbolic notation
    syms s;
    syms n1;
    syms n2;
    
    N_temp = [n1; n2]; %initialising N in a symbolic matrix
    L_temp = A - N_temp*C; %calculating symbolic L with N_tempt
    char_eq = det(s*eye(2,2)-L_temp); %working out characteric equation with sym n1 and n2
    
    coeff_delta_char = flip(coeffs(char_eq,s)); %finding coefficients of the characteristic equation with sym n1 and n2
    coeff_delta = flip(coeffs(delta,s)); %finding coefficients of the characteristic equation
    
    temp_sym = [symvar(coeff_delta_char(2)) symvar(coeff_delta_char(3))]; %setting temporary matrix to see which sym variable is stored in which coefficient
    
    %based on which coefficient has n1 and n2, equating the two characteristic equations and solving for n1 and n2
    if(n2 == temp_sym(1))
        n2_value = solve(coeff_delta_char(2) == coeff_delta(2));
        n1_value = solve(subs(coeff_delta_char(3),n2,n2_value) == coeff_delta(3));
        N = [n1_value; n2_value]; %returning N
    else
        n1_value = solve(coeff_delta_char(2) == coeff_delta(2));
        n2_value = solve(subs(coeff_delta_char(3),n1,n1_value) == coeff_delta(3));
        N = [n1_value; n2_value]; %returning N
    end
    L = subs(L_temp,[n1 n2],[n1_value n2_value]); %subbing in the n values for L
    M = B-N*D; %calculating M, based on matrices B, N and D
end
