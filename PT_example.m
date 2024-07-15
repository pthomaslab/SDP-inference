%%% This is an example of bounding parameters in a post-transcriptional gene regulation model (see Fig.4 in paper)
%%% We demonstrate how to combine two independent datasets in which different species are observed

k = [6,0.8,5,0.5,1]; % true parameters
kn = [1,3,5]; % which ones are known
n = 2; % number of species
n_ukn = 5 - length(kn); % number of unknowns


% call the function with species 1 is unobserved
[g1, M0s1, M1s1, M2s1, yT1, zs1, list1, list_21, missing1] = mmt(k, kn, d, n, 1);
% call the function with species 2 is unobserved
[g2, M0s2, M1s2, M2s2, yT2, zs2, list2, list_22, missing2] = mmt(k, kn, d, n, 2);


%%%%%%% setting up the constrained sets %%%%%%%


% all variables in the set with YALMIP, including parameters and created variable to represent cross-terms
k = sdpvar((n_ukn+length(list1)+length(list_21)),1); 

g1 = eval(char(str2sym(g1))); % make the mmt. eqn. constraints with the YALMIP variables

Fr = [k >= 0]; % trivial constraint
Fr = Fr + [g1 == zeros(length(g1))]; % mmt. eqn. constraints


%%%%%% k_j * y_l <= z_j <= k_j * y_u (see paper for details)
% using mmt. intervals Y_lower1, Y_upper1 from observed data 1 where
% species 1 considered to be unobserved
Fr = Fr + [Y_lower1(setdiff(1:length(yT1),missing1))' ...
    <= flip(k((end-length(list_21)/(n_ukn+1)+1):end)) ...
    <= Y_upper1(setdiff(1:length(yT1),missing1))'];
for j = 1:n_ukn
    Fr = Fr + [Y_lower1(setdiff(1:length(yT1),missing1))' .* k(n_ukn-j+1) ...
        <= flip(k((end-(j+1)*length(list_21)/(n_ukn+1)+1):(end-j*length(list_21)/(n_ukn+1)))) ...
        <= Y_upper1(setdiff(1:length(yT1),missing1))' .* k(n_ukn-j+1)];
end
% using mmt. intervals Y_lower2, Y_upper2 from observed data 2 where
% species 2 considered to be unobserved
Fr = Fr + [Y_lower2(setdiff(1:length(yT2),missing2))' ... 
    <= flip(k((end-length(list_21)/(n_ukn+1)+1):end)) ...
    <= Y_upper2(setdiff(1:length(yT2),missing2))'];
for j = 1:n_ukn
    Fr = Fr + [Y_lower2(setdiff(1:length(yT2),missing2))' .* k(n_ukn-j+1) ...
        <= flip(k((end-(j+1)*length(list_21)/(n_ukn+1)+1):(end-j*length(list_21)/(n_ukn+1)))) ...
        <= Y_upper2(setdiff(1:length(yT2),missing2))' .* k(n_ukn-j+1)];
end
%%%%%%


%%%%% mmt. matrices constraints
m0 = eval(char(str2sym(M0s1{end})));
m1 = eval(char(str2sym(M1s1{end})));
m2 = eval(char(str2sym(M2s1{end})));
Fr = Fr + [m0 >= 0, m1 >= 0, m2 >= 0];
%%%%%


%%%%% bound k(1) as an example (note k(1) is the second parameter k_2 here)
sol_l = solvesdp(Fr, k(1));
k_2_lower_bound = value(k(1));
sol_u = solvesdp(Fr, -k(1));
k_2_upper_bound = value(k(1));




function [g, M0s, M1s, M2s, yT, zs, list, list_2, missing] = mmt(k, kn, d, n, uknspecies)

%%% setting constants %%%

db = 2; % highest order of b_j, needs manual change
order = d - db + 1; % order of mmt. eqn. up to which we consider
md = nchoosek(n+d,n);
n_k = length(k);
n_k_kn = length(kn);
ukn = 1:n_k;
ukn(kn) = [];
n_k_ukn = length(ukn);

%%% moment equations %%%

syms x1 x2 k1 k2 k3 k4 k5
g = string(zeros(order+1, order+1));
for i = 0:order
    for j = 0:(order-i) 
        % need manual change depending on the model
        g(i+1,j+1) = string(expand(k1 * ((x1 + 1)^i - x1^i) * x2^j + ...
            k2 * x1 * ((x1 - 1)^i - x1^i) * x2^j + ...
            k3 * ((x2 + 2)^j - x2^j) * x1^i + ...
            k4 * x2 * (x2 - 1)* ((x2 - 2)^j - x2^j) * x1^i + ...
            k5 * x1 * x2 * ((x1 - 1)^i * (x2 - 1)^j - x1^i * x2^j)));
    end
end
g = str2sym(g);
for z = 1:n_k_kn
    % substituting known parameters
    g = expand(subs(g, "k"+kn(z), k(kn(z))));
end


yT = string(zeros(1,md)); % vector that contains all involved moments
pos = 0;
for i = 0:(n+d-2) 
    for j = 0:i
        yT(pos+j+1) = "x1^"+(i-j)+"*x2^"+j; 
    end
    pos = pos + i + 1;
end

% this is task specific, it lets you to set which species is unobserved
if uknspecies == 1
    index = find(contains(yT, "x1^0"));
elseif uknspecies == 2
    index = find(contains(yT, "x2^0"));
end

% change notations
y_sym = yT;
for i = 2:length(yT)
    y_sym(i) = replace(replace(y_sym(i), "x1^0*", ""), "*x2^0", "");
    y_sym(i) = replace(replace(y_sym(i), "x1^1", "x1"), "x2^1", "x2");
end

% find which moments involve the unobserved species
missing = setdiff(1:length(y_sym), index);
for l = 1:length(missing)
    zs(l) = "z"+l;
end

%%% moment matrices %%%

M0s = cell(d-db, 1);
M1s = cell(d-db, 1);
M2s = cell(d-db, 1);

for i = db:(d-1)
    % multi-index matrix
    alpha0 = zeros(2,1);
    pos = 1;
    for j = 0:floor(i/2)
        for k = 0:(floor(i/2)-j)
            alpha0(:,pos) = [j k];
            pos = pos + 1;
        end
    end
    n0 = width(alpha0);
    M0 = string(zeros(n0, n0));
    % create mmt. matrix
    for j = 1:n0
        for k = 1:n0
            M0(j, k) = "x1^"+(alpha0(1,j)+alpha0(1,k))+"*x2^"+(alpha0(2,j)+alpha0(2,k));
        end
    end
    % multi-index matrix
    alpha1 = zeros(2,1);
    pos = 1;
    for j = 0:floor((i-1)/2)
        for k = 0:(floor((i-1)/2)-j)
            alpha1(:,pos) = [j k];
            pos = pos + 1;
        end
    end
    n1 = width(alpha1);
    M1 = string(zeros(n1, n1));
    M2 = string(zeros(n1, n1));
    % create mmt. matrices
    for j = 1:n1
        for k = 1:n1
            M1(j, k) = "x1^"+(alpha1(1,j)+alpha1(1,k)+1)+"*x2^"+(alpha1(2,j)+alpha1(2,k));
            M2(j, k) = "x1^"+(alpha1(1,j)+alpha1(1,k))+"*x2^"+(alpha1(2,j)+alpha1(2,k)+1);
        end
    end
    
    % change notations & some substitutions
    for j = flip(1:length(yT))
        M0 = replace(M0, yT(j), "y"+j);
        M1 = replace(M1, yT(j), "y"+j);
        M2 = replace(M2, yT(j), "y"+j);
    end
  
    for l = flip(1:length(missing))
        M0 = replace(M0, "y"+missing(l), "z"+l);
        M1 = replace(M1, "y"+missing(l), "z"+l);
        M2 = replace(M2, "y"+missing(l), "z"+l);
    end
    M0(1,1) = 1;

    for j = flip(1:n_k_ukn)
        M0 = replace(M0, "k"+ukn(j), "k("+j+")");
        M1 = replace(M1, "k"+ukn(j), "k("+j+")");
        M2 = replace(M2, "k"+ukn(j), "k("+j+")");
    end

    M0s{i-db+1} = M0;
    M1s{i-db+1} = M1;
    M2s{i-db+1} = M2;
end

% creating variables to substitue the cross terms between mmt.s and par.s
% to make the mmt. eqn. constraints linear
% (this is a bit complicated in this specific task as we are combining data)
list = repmat(flip(zs(1:(end))), [1,(n_k_ukn+1)]);
for i = 1:(length(zs))
    for j =1:n_k_ukn
        list((length(zs))*(j-1)+i) = "k("+j+")*" + list((length(zs))*(j-1)+i);
    end
end
list_2 = repmat(flip("y"+setdiff(1:length(yT),missing)), [1,(n_k_ukn+1)]);
for i = 1:(length(yT) - length(zs))
    for j =1:n_k_ukn
        list_2((length(yT) - length(zs))*(j-1)+i) = "k("+j+")*" + list_2((length(yT) - length(zs))*(j-1)+i);
    end
end
% substituting
for i = 1:length(list)
    for j = 1:length(M0s)
        M0s{j} = replace(M0s{j}, list(i), "k("+(i+n_k_ukn)+")");
        M1s{j} = replace(M1s{j}, list(i), "k("+(i+n_k_ukn)+")");
        M2s{j} = replace(M2s{j}, list(i), "k("+(i+n_k_ukn)+")");
    end
end
for i = 1:length(list_2)
    for j = 1:length(M0s)
        M0s{j} = replace(M0s{j}, list_2(i), "k("+(i+n_k_ukn+length(list))+")");
        M1s{j} = replace(M1s{j}, list_2(i), "k("+(i+n_k_ukn+length(list))+")");
        M2s{j} = replace(M2s{j}, list_2(i), "k("+(i+n_k_ukn+length(list))+")");
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
for i = 1:length(list_2)
    g = replace(g, list_2(i), "k("+(i+n_k_ukn+length(list))+")");
end
end

