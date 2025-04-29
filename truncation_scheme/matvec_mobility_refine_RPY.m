function u = matvec_mobility_refine_RPY(tau,q,rin_base_c,rin_base_f,RPYc,RPYcut,RPYfine)
%matvec using the RPY tensor

P = size(q,1); %number of particles
Nc = size(rin_base_c,1); %points per particle on proxy surface, coarse grid 
Nf = size(rin_base_f,1); %points per particle on proxy surface, fine grid 
 
tau_c = tau(1:3*Nc*P);
tau_cut = tau(3*Nc*P+1:6*Nc*P);
tau_fine = tau(6*Nc*P+1:6*Nc*P+3*Nf*P);
U = tau(6*Nc*P+3*Nf*P+1:end);

u = zeros(6*Nc*P+3*Nf*P+6*P,1); 

u(1:3*Nc*P) = RPYc*tau_c;
u(3*Nc*P+1:6*Nc*P) = RPYcut*tau_cut;
u(6*Nc*P+1:6*Nc*P+3*Nf*P) = RPYfine*tau_fine; %evaluated in different points...

Bc = getKmat(rin_base_c,[0,0,0]);
Bf = getKmat(rin_base_f,[0,0,0]);
for k = 1:P
    u(3*Nc*(k-1)+1:3*Nc*k) = u(3*Nc*(k-1)+1:3*Nc*k)-Bc*U((k-1)*6+1:6*k);
    u(3*Nc*P+1+(k-1)*3*Nc:3*Nc*P+k*3*Nc) = u(3*Nc*P+1+(k-1)*3*Nc:3*Nc*P+k*3*Nc)-Bc*U((k-1)*6+1:6*k);
    u(6*Nc*P+1+(k-1)*3*Nf:6*Nc*P+k*3*Nf) = u(6*Nc*P+1+(k-1)*3*Nf:6*Nc*P+k*3*Nf)-Bf*U((k-1)*6+1:6*k);
    u(6*Nc*P+3*Nf*P+1+(k-1)*6:6*Nc*P+3*Nf*P+k*6) = -Bc'*tau_c((k-1)*3*Nc+1:k*3*Nc) +...
        Bc'*tau_cut((k-1)*3*Nc+1:k*3*Nc)-Bf'*tau_fine((k-1)*3*Nf+1:k*3*Nf);
end

end