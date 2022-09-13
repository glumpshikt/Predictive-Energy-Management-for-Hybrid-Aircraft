%% jacobian sparse for IPOPT (MPC)

function g = jacobiansparseMPC(auxdata)
N = cell2mat(auxdata(1));
K = cell2mat(auxdata(22));
I = speye(N);
Psi = sparse(tril(ones(N), -1));

g = zeros(8*N*K+K-1, 6*N*K);

for k=1:K
    %constraint 1
    g(1+(k-1)*8*N: N+(k-1)*8*N, 1+(k-1)*6*N:2*N+(k-1)*6*N) = [I Psi];
    
    %constraint 2
    g(N+1+(k-1)*8*N: 2*N+(k-1)*8*N, 1+(k-1)*6*N:2*N+(k-1)*6*N) = [I I];
    g(N+1+(k-1)*8*N: 2*N+(k-1)*8*N, 3*N+1+(k-1)*6*N: 6*N+(k-1)*6*N) = [I I I];
    
    %constraint 3
    g(2*N+1+(k-1)*8*N: 3*N+(k-1)*8*N, 2*N+1+(k-1)*6*N: 4*N+(k-1)*6*N) = [I Psi];
    
    %constraint 4
    g(3*N+1+(k-1)*8*N: 4*N+(k-1)*8*N, 2*N+1+(k-1)*6*N: 3*N+(k-1)*6*N) = I;
    g(3*N+1+(k-1)*8*N: 4*N+(k-1)*8*N, 5*N+1+(k-1)*6*N: 6*N+(k-1)*6*N) = I;
    
    %constraint 5
    g(4*N+1+(k-1)*8*N: 5*N+(k-1)*8*N, 3*N+1+(k-1)*6*N: 6*N+(k-1)*6*N) = [I I I];
    
    %constraint 6
    g(5*N+1+(k-1)*8*N: 6*N+(k-1)*8*N, 3*N+1+(k-1)*6*N: 6*N+(k-1)*6*N) = [I I I];
    
    %constraint 7
    g(6*N+1+(k-1)*8*N: 7*N+(k-1)*8*N, 2*N+1+(k-1)*6*N: 3*N+(k-1)*6*N) = I;
    g(6*N+1+(k-1)*8*N: 7*N+(k-1)*8*N, 4*N+1+(k-1)*6*N: 5*N+(k-1)*6*N) = I;
    
    %constraint 8
    g(7*N+1+(k-1)*8*N: 8*N+(k-1)*8*N, 3*N+1+(k-1)*6*N: 6*N+(k-1)*6*N) = [I I I];
end

for k=1:K-1
    g(8*K*N+k, 3*N+1) = 1;
    g(8*K*N+k, 3*N+1+k*6*N) = 1;
end

g = sparse(g);
end