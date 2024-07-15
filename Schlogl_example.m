%%% This is an example of bounding parameters in a Schlogl model (see Fig.2(c) in paper)

load('Example_data.mat') % load example data
Y_l = Y_l_20000;
Y_u = Y_u_20000;

k0 = [2,3,1,4]; % true parameters
kn = [2,4]; % which ones are known
n_ukn = 4-length(kn); % number of unknowns 

d = 5; % chosen order of constrained set, at least 3 here

[g, g0, yT, zs, list, missing] = mmt(k0, kn, d);

% all variables in the set with YALMIP, including parameters and created variable to represent cross-terms
k = sdpvar(length(list)+(4-length(kn)),1);

g = eval(char(str2sym(g))); % make the mmt. eqn. constraints with the YALMIP variables

Fr = [k >= 0]; % trivial constraint
Fr = Fr + [g == zeros(length(g),1)]; % mmt. eqn. constraints

%%%%%% k_j * y_l <= z_j <= k_j * y_u (see paper for details)
Fr = Fr + [Y_l(2:(d+1)) <= flip(k((end-length(list)/(n_ukn+1)+1):end)) <= Y_u(2:(d+1))];
for j = 1:n_ukn
    Fr = Fr + [Y_l(2:(d+1)) .* k(n_ukn-j+1) <= flip(k((end-(j+1)*length(list)/(n_ukn+1)+1):(end-j*length(list)/(n_ukn+1)))) <= Y_u(2:(d+1)) .* k(n_ukn-j+1)];
end

ops = sdpsettings('solver','mosek');
% example: maximizing k1, the first unknown parameter
optimize(Fr,-k(1),ops)
value(k(1))

function [g, g0, yT, zs, list, missing] = mmt(k, kn, d)
n = 1;
db = 3; % highest order of b_j, needs manual change
order = d - db + 1;
n_k = length(k);
n_k_kn = length(kn);
ukn = 1:n_k;
ukn(kn) = [];
n_k_ukn = length(ukn);

% moment equations
g = string(zeros(order+1, 1));
for i = 0:order
    syms x k1 k2 k3 k4 
    % Need to manually change the propensities
    g(i+1) = string(expand(k1 * x * (x - 1) * ((x + 1)^i - x ^ i) + ...
        k2 * x * (x - 1) * (x - 2) * ((x - 1)^i - x ^ i) + ...
        k3 * ((x + 1)^i - x ^ i) + k4 * x * ((x - 1)^i - x ^ i)));
end
g = str2sym(g);
for z = 1:n_k_kn
    g = expand(subs(g, "k"+kn(z), k(kn(z))));
end
g = string(g);

yT = string(zeros(1,n+d));
pos = 0;
for i = 0:(n+d-1) % highest order in y is n+d-1
    yT(pos + 1) = "x^"+(i);
    pos = pos + 1;
end
index = 1:length(yT);

y_sym = yT;
for i = 2:length(yT)
    y_sym(i) = replace(y_sym(i), "x^0*", "");
    y_sym(i) = replace(y_sym(i), "x^1", "x");
end

for l = flip(2:length(y_sym))
    g = replace(g, y_sym(l), "y"+l);
end

missing = index(2:end);
for l = 1:length(missing)
    g = replace(g, "y"+missing(l), "z"+l);
    zs(l) = "z"+l;
end

list = repmat(flip(zs(1:(end))), [1,(n_k_ukn+1)]);
for i = 1:(length(zs))
    for j =1:n_k_ukn
        list((length(zs))*(j-1)+i) = "k("+j+")*" + list((length(zs))*(j-1)+i);
    end
end

g = string(g);
for l = flip(2:length(y_sym))
    g = replace(g, y_sym(l), "y"+l);
end
for l = flip(1:length(missing))
    g = replace(g, "y"+missing(l), "z"+l);
end
for l = flip(1:length(missing))
    g = replace(g, "y"+missing(l), "z"+l);
end
for i = 1:length(ukn)
    g = replace(g, "k"+ukn(i),"k("+i+")");
end
g0 = g;
for i = 1:length(list)
    g = replace(g, list(i), "k("+(i+n_k_ukn)+")");
end

end
