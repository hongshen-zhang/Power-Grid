function [resi_res] = predict_resilience(mpc,varargin)
% PREDICT_RESILIENCE creates the resilience prediction for load bus 
%   
%   Examples:
%      
%      PREDICT_RESILIENCE(mpc);
%
%      PREDICT_RESILIENCE(mpc,'REAL_PMAX',Real_Pmax); 
%           Real_Pmax : numel(g_index) * 2 double
%           Real_Pmax(:,1) = g_index;
%           Real_Pmax(:,2) = P addition bound for Generator;
%
%      PREDICT_RESILIENCE(mpc,'LOAD_DIRECTION',Load_direction);
%           Load_direction : numel(load_index) * 2 double
%           Load_direction(:,1) = load_index;
%           Load_direction(:,2) = load perturbation direction; 
%           (100 + 100j & baseMVA = 100 ==> S = -1-1j)
%
%      PREDICT_RESILIENCE(mpc,'REAL_RESI',Real_Resi);
%           Real_Resi : numel(p_index) * 2 double
%           Real_Resi(:,1) = p_index;
%           Real_Resi(:,2) = Real Resilience for load bus
%
tic;
[Y,g_index,p_index,S] = ini_posi(mpc);
[pf_vol] = get_pf_vol(mpc);
Real_Pmax = [];
Real_Resi = [];
consider_self_power = 0;
okargs = {'real_pmax','load_direction','real_resi'};
for j = 1 : 2 : nargin-1
    pname = varargin{j};
    pval = varargin{j+1};
    k = strmatch(lower(pname), okargs);
    if isempty(k)
        error('Bioinfo:UnknownParameterName',...
            'Unknown parameter name: %s.',pname);
    elseif length(k)>1
        error('Bioinfo:AmbiguousParameterName',...
            'Ambiguous parameter name: %s.',pname);
    else
        switch(k) 
            case 1 %Real Pmax
                if numel(g_index) == numel(pval(:,1))
                    Real_Pmax = pval;
                end
            case 2 % Load and Direction
                p_index = pval(:,1);
                S = pval(:,2);
                consider_self_power = 1;
            case 3 % Real Resilience
                Real_Resi = pval;
        end
    end
end
lambda_set = [];
Y1_set = [];
Y2_set = [];
[weight,Y_eq]= Father_matrix(Y,g_index,p_index);
vol_g = pf_vol(g_index);
Y1_set(:,1) = Y_eq;
Y2_set(:,1) = weight * vol_g;
[lambda_set(:,1)] = predict_ultra(weight,Y_eq,vol_g,S);
% Minus initial load
[weight,Y_eq,vol_g,S] = minus_initial_PQ(mpc,weight,Y_eq,vol_g,S,p_index,consider_self_power);
Y1_set(:,2) = Y_eq;
Y2_set(:,2) = weight * vol_g;
[lambda_set(:,2)] = predict_ultra(weight,Y_eq,vol_g,S);
% With Real Pmax
if numel(Real_Pmax)
    [weight,Y_eq,vol_g,S,c2] = ...
        minus_Pmax_effect(mpc,weight,Y_eq,vol_g,S,p_index,Real_Pmax);
    Y1_set(:,3) = Y_eq;
    Y2_set(:,3) = weight * vol_g;
    resi_res.c2 = c2;
    [lambda_set(:,3)] = predict_ultra(weight,Y_eq,vol_g,S);
end
% Compare real result
if numel(Real_Resi)
    result = [];
    [~,load_posi] = ismember(Real_Resi(:,1),p_index);
    result.list = [Real_Resi(:,1),lambda_set(load_posi,:),...
        Real_Resi(:,2)];
    [result.ave1,result.med1] = average_median_error(Real_Resi(:,2),lambda_set(load_posi,1));
    [result.ave2,result.med2] = average_median_error(Real_Resi(:,2),lambda_set(load_posi,2));
    if numel(lambda_set(1,:)) == 3
        [result.ave3,result.med3] = average_median_error(Real_Resi(:,2),lambda_set(load_posi,3));
    end
    resi_res.result = result;
end
t_end = toc;
resi_res.lambda_set = lambda_set;
resi_res.Y1_set = Y1_set;
resi_res.Y2_set = Y2_set;
resi_res.S = S;
resi_res.runtime_second = t_end;
end

function [lam] = predict_ultra(weight,Y_eq,vol_g,S)
% PREDICT ULTRA : Prediction with voltage and shunt
%   1.weight : numel(p_index) * numel(g_index) complex
%   2.Y_eq : numel(p_index) * 1 complex
%   3.vol_g : numel(g_index) * 1 double
%   3.S : numel(p_index) * 1 complex
lam = abs(weight * vol_g).^2 ./ (2*(abs(Y_eq.*S)-real(Y_eq.*S)));
end

function [average_error,median_error] = average_median_error(real,predict)
% AVERAGE MEDIAN ERROR : Calculate percentage error number for two list
% err list = abs(real - predict ./ real)
%   INPUT:
%       1. real n*1 double
%       2. predict n*1 double
%   OUTPUT:
%       1.average_error 1*1 double
%       2.median_error 1*1 double
%
average_error = mean(100 * abs(real - predict) ./ real);
median_error = median(100 * abs(real - predict) ./ real);
end

function [pf_vol] = get_pf_vol(mpc)
% GET PF VOL do two things :
% 1. Whether the case can runpf.
% 2. Get the voltage of runpf of load.
%
%   INPUT : mpc
%   OUTPUT : numel(mpc.bus(:,8)) * 1
%       RUNPF VOLTAGE of all bus
%
mpopt = mpoption('OUT_ALL',0,'VERBOSE',0);
pf = runpf(mpc,mpopt);
if ~pf.success
    error('This mpc can not runpf directly.');
end
pf_vol = pf.bus(:,8);
end

function [weight,Y_eq]= Father_matrix(Y,g_index,p_index)
% FATHER MATRIX creates equivalent system for mpc
%   INPUT : 
%       1. Y : admittance matrix
%       2. g_index : generator index
%       3. p_index : load index 
%   OUTPUT: 
%       1. weight : numel(p_index) * numel(g_index)
%           Equivalent matrix for load to each gen
%
%       2. Y_eq : numel(p_index) 
%           Equivalent admittance for load the whole gen
%
N = size(Y,1);
load_index = setdiff(1:N,g_index);
[~,p_loc] = ismember(p_index,load_index);
B = (full(Y(load_index,load_index)))^(-1);
C = full(Y(load_index,g_index));
weight = -B(p_loc,:) * C ./ (diag(B(p_loc,p_loc)));
Y_eq = diag(eye(numel(p_loc)) ./ diag(B(p_loc,p_loc)));
end

function [Y,g_index,p_index,S] = ini_posi(mpc)
% INI_POSI creates basic information for mpc
%   INPUT : 
%       1. mpc 
%   OUTPUT: 
%       1. Y : Admittance matrix 
%       2. g_index : index of Generator(PV and slack bus)
%       3. p_index : index of positive Load (P > 0 and Q > 0)
%       4. S : Power of load bus (S = - P - Q j)
%
Y = makeYbus(mpc);
g_bus = find(mpc.bus(:,2) ~=1);
g_index = intersect(unique(mpc.gen(:,1)),mpc.bus(g_bus,1));
Positive_bus = intersect(find(mpc.bus(:,3)>0), find(mpc.bus(:,4)>0));
p_index = intersect(find(mpc.bus(:,2)==1), Positive_bus);
S = (-mpc.bus(p_index,3) - mpc.bus(p_index,4) * j)./ mpc.baseMVA;
end


function [weight,Y_eq,vol_g,S] = minus_initial_PQ(mpc,weight,Y_eq,vol_g,S,p_index,csp)
% MINUS_INITIAL_PQ decrease the impact of initial load 
% INPUT :
%   1.mpc 
%   2.weight : numel(p_index) * numel(g_index) 
%       Equivalent admittance from load to each generator
%   3.Y_eq : numel(p_index) * 1
%       Equivalent admittance from load to the set of generator
%       Note that sum(weight,2) ~= Y_eq is due to the impact of
%       shunt(load shunt, branch shunt, angle and ratio)
%   4.vol_g : numel(g_index) * 1
%       Voltage of each generator
%   5.S : numel(p_index) * 1 
%       Indicate that the direction of perturbation load
%   6.p_index : numel(p_index) * 1
%   7.csp : consider_self_power (boolean)
%
% OUTPUT:
%   1.weight : numel(p_index) * numel(g_index)
%       New Equivalent admittance after minus the impact of load
%   2.Y_eq : numel(p_index) * 1
%       New equivalent admittance
%   3.vol_g : numel(g_index) * 1 (Not update)
%   4. S : numel(p_index) * 1 (Not update)
%
[Y,g_index,~,~] = ini_posi(mpc);
for loop = 1 : numel(mpc.gen(:,1))
    mpc.bus(mpc.gen(loop,1),3) = mpc.bus(mpc.gen(loop,1),3) - mpc.gen(loop,2);
    mpc.bus(mpc.gen(loop,1),4) = mpc.bus(mpc.gen(loop,1),4) - mpc.gen(loop,3);
end
load_index = setdiff(mpc.bus(:,1),g_index);
Nonzero_bus = union(find(mpc.bus(:,3)~=0), find(mpc.bus(:,4)~=0));
nonzero_index = intersect(find(mpc.bus(:,2)==1), Nonzero_bus);
[~,p_posi] = ismember(p_index,load_index);
[~,nonzero_posi] = ismember(nonzero_index,load_index);
Y_load_load = Y(load_index,load_index);
inv_Y = full((Y_load_load)^(-1));
diag_Y = diag(inv_Y);
R_eq_load =  diag_Y * ones(1,numel(load_index)) ...
    +  conj((diag_Y * ones(1,numel(load_index)))') ....
    - inv_Y - conj(inv_Y');
clear inv_Y
Y_2_mat = 1 ./ R_eq_load(p_posi,nonzero_posi);
clear R_eq_load
[weight_load,Y_load] = Father_matrix(Y,g_index,load_index);
Y_1_mat = ones(numel(p_index),1) * conj(Y_load(nonzero_posi)');
Y_3_mat = Y_load(p_posi) * ones(1,numel(nonzero_index));
r_star_1_mat = (1 ./ Y_2_mat + 1 ./ Y_3_mat - 1 ./ Y_1_mat) / 2;
r_star_2_mat = (1 ./ Y_1_mat + 1 ./ Y_3_mat - 1 ./ Y_2_mat) / 2;
r_star_3_mat = (1 ./ Y_1_mat + 1 ./ Y_2_mat - 1 ./ Y_3_mat) / 2;
clear Y_1_mat Y_2_mat Y_3_mat
delta_2_mat = 1 ./ (r_star_2_mat + r_star_3_mat + r_star_2_mat .* r_star_3_mat  ./ r_star_1_mat);
delta_3_mat = 1 ./ (r_star_1_mat + r_star_3_mat + r_star_1_mat .* r_star_3_mat  ./ r_star_2_mat);
clear r_star_1_mat r_star_2_mat r_star_3_mat
S_1_line = -(mpc.bus(nonzero_index,3) +mpc.bus(nonzero_index,4) * j) ...
    / mpc.baseMVA;
S_1_mat = conj((S_1_line * ones(1,numel(p_index)))');
S_2_mat = S * ones(1,numel(nonzero_index));
for loop = 1 : numel(Y_eq)
    ind_false = find(isnan(delta_2_mat(loop,:)));
    if ind_false
        delta_2_mat(loop,ind_false) = Y_eq(loop);
        delta_3_mat(loop,ind_false) = Y_eq(loop) * 9e1;
        if ~csp
            S_1_mat(loop,ind_false) = S_1_mat(loop,ind_false) - S(loop);
        end
    end
end
[resi_1_mat,resi_2_mat] = calculate_add_resi_mat...
    (delta_2_mat,delta_3_mat,S_1_mat,S_2_mat);
raw_pre_wide = predict_without_p(Y_eq,S) * ones(1,numel(nonzero_index));
resi_1_mat(isnan(resi_1_mat)) = 0;
resi_2_mat(isnan(resi_2_mat)) = 0;
delta_resi_1_2 = resi_1_mat - resi_2_mat;
zero_elements = find(abs(delta_2_mat) > 1e4 * abs(delta_3_mat));
not_exact = find(resi_1_mat > 2 * raw_pre_wide);
wrong_index = find(abs(delta_3_mat) > 1e2 * abs(delta_2_mat));
delta_resi_1_2(zero_elements) = 0;
delta_resi_1_2(not_exact) = 0;
delta_resi_1_2(wrong_index) = abs(S_1_mat(wrong_index)) ./  abs(S_2_mat(wrong_index));
clear  delta_2_mat delta_3_mat S_1_mat S_2_mat 
ori_weight_vol_g  = weight * vol_g;
split_ratio = abs(weight_load(nonzero_posi,:) ./ (Y_load(nonzero_posi) * ones(1,numel(g_index))));
delta_weight_resi = delta_resi_1_2 * split_ratio;
S_wide = S * ones(1,numel(g_index));
theo_resi_wide = predict_without_p(weight,S_wide) .* ((vol_g * ones(1,numel(p_index))).^2)';
delta_weight_resi = min(delta_weight_resi,theo_resi_wide);
weight = weight .* (theo_resi_wide - delta_weight_resi) ./ theo_resi_wide;
weight(isnan(weight)) = 0;
new_weight_vol_g = weight * vol_g;
Y_eq = Y_eq .* (new_weight_vol_g ./ori_weight_vol_g);
end


function [weight,Y_eq,vol_g,S,c2] = minus_Pmax_effect(mpc,weight,Y_eq,vol_g,S,p_index,max_P)
% MINUS_INITIAL_PQ decrease the impact of initial load 
% INPUT :
%   1.mpc 
%   2.weight : numel(p_index) * numel(g_index) 
%       Equivalent admittance from load to each generator
%   3.Y_eq : numel(p_index) * 1
%       Equivalent admittance from load to the set of generator
%       Note that sum(weight,2) ~= Y_eq is due to the impact of
%       shunt(load shunt, branch shunt, angle and ratio)
%   4.vol_g : numel(g_index) * 1
%       Voltage of each generator
%   5.S : numel(p_index) * 1 
%       Indicate that the direction of perturbation load
%   6.p_index : numel(p_index) * 1
%   7.max_P : numel(g_index) * 1
%       Indicate the direct 1(Add 100 + 0j) bound of generator
%
% OUTPUT:
%   1.weight : numel(p_index) * numel(g_index)
%       New Equivalent admittance after minus the impact of load
%   2.Y_eq : numel(p_index) * 1
%       New equivalent admittance
%   3.vol_g : numel(g_index) * 1 (Not update)
%   4. S : numel(p_index) * 1 (Not update)
%   5. c2(As required)
%
max_P_posi = max_P(:,2);
[Y,g_index,~,~] = ini_posi(mpc);
S_wide = S * ones(1,numel(g_index));
if max(max_P_posi) == 0 
    max_P_posi(max_P_posi==0) = 1e4;
else
    max_P_posi(max_P_posi==0) = max(max_P_posi) *2;
end
diag_Y = diag(full(Y));
Y(logical(eye(size(Y))))=0;
g_ratio = -sum(Y(g_index,:),2) ./ diag_Y(g_index);
max_P_posi = real(max_P_posi ./ (vol_g.^2) ./ g_ratio);
max_P_posi_wide = ones(numel(p_index),1) * max_P_posi';
weight_per = abs(weight ./ Y_eq);
max_P_posi_wide_per = weight_per .*max_P_posi_wide;
clear max_P_wide weight_per
[c2] = c2_calculation(weight,S_wide,vol_g,max_P_posi_wide_per);
[lambda_wide,lambda_origin] = ...
    predict_p_posi_final(weight,S_wide,max_P_posi_wide_per);
clear max_P_wide_per p_slack_wide S_wide
weight_ratio = lambda_wide ./ lambda_origin;
weight_ratio(isnan(weight_ratio)) = 1;
weight_ori = weight * vol_g;
weight = weight .* weight_ratio;
weight_new = weight * vol_g;
Y_eq = Y_eq .* weight_new ./ weight_ori;
end

function [c2] = c2_calculation(Y_eq,S,vol_g,p_posi)
% c2_calculation calculate the list of c2
% where a and b are any natural number
% INPUT : 
%   1.Y_eq : a * b complex 
%   2.S: a * b complex
%   3.vol_g : a * 1 double
%   4.p_posi : a * b double 
%
% OUTPUT:
%   1.c2: a * 1 double
%
[lambda_without_p] = predict_without_p(Y_eq,S);
lambda_with_p = lambda_without_p;

ind_neg_real_S = intersect(find(real(S) < 0),find(real(Y_eq)>=0));
[slack_posi] = calculate_p_slack(Y_eq(ind_neg_real_S),...
    S(ind_neg_real_S)); 
ind_big_posi = find(p_posi(ind_neg_real_S) < slack_posi);
[lam_big_posi,~] = calculate_p_root(Y_eq(ind_neg_real_S(ind_big_posi )),...
    S(ind_neg_real_S(ind_big_posi )),...
    p_posi(ind_neg_real_S(ind_big_posi )));
lambda_with_p(ind_neg_real_S(ind_big_posi )) = lam_big_posi;

ind_small_posi = find(p_posi(ind_neg_real_S) >= slack_posi);
size_kuoda = p_posi(ind_neg_real_S(ind_small_posi))  ./  ...
    slack_posi(ind_small_posi);
lambda_with_p(ind_neg_real_S(ind_small_posi)) = ...
    lambda_with_p(ind_neg_real_S(ind_small_posi))  .*...
    size_kuoda;

ratio =  lambda_without_p ./ lambda_with_p;
ratio(isnan(ratio)) = 1;
old_weight_vol_g = Y_eq * vol_g;
new_weight_vol_g = (Y_eq .* ratio) * vol_g;
c2 = abs(new_weight_vol_g ./ old_weight_vol_g);
end 


function [resi_1,resi_2] = calculate_add_resi_mat(Y1,Y2,S1,S2)
% CALCULATE ADD RESI MAT estimate the addition resi 
% INPUT :
%   1.Y1 : a * b complex (generator to initial load)
%   2.Y2 : a * b complex (initial load to perturbation load)
%   3.S1 : a * b complex (initial load size)
%   4.S2 : a * b complex (perturbation load size
%
% OUTPUT:
%   1.resi_1 : a * b double (predict resilience without initial load)
%   2.resi_2 : a * b double (esitimation resilience with initial load)
Y_eq_ori = (Y1 .* Y2) ./ (Y1 +Y2);
[resi_1] = predict_without_p(Y_eq_ori,S2);
V_2c = 0.5 - j * imag(resi_1 .*S2 ./ conj(Y_eq_ori));
V_1c = (Y1 + Y2 .* V_2c) ./ (Y1 + Y2);
V_ratio = V_1c ./ V_2c;
K = 1 + Y2 ./ Y1 - Y2  ./(Y1 .* V_ratio);
[v1] = solve_standard_v(Y1 ./ conj(K),S1) ./ K;
[resi_2] = abs(resi_1 .* (v1.^2) ./ ((V_1c).^2));
end

function [v] = solve_standard_v(Y,S)
% SOLVE STANDARD V calculate the voltage with gen = 1, Y and initial load S
% INPUT:
%   1. Y : a * b complex (generator to load admittance)
%   2. S : a * b complex (load power)
v = ( abs(Y).^2 + sqrt((2 * real(S.*Y) +abs(Y).^2).^2 - 4 * abs(S.*Y).^2)...
    - 2 * imag(S.*Y) * j) ./ ( 2 * abs(Y).^2);
end

function [p_slack] = calculate_p_slack(Y_eq,S)
% CALCULATE_P_SLACK calculate the slack of generator
% where a and b are any natural number
% INPUT : 
%   1.Y_eq : a * b complex 
%   2.S : a * b complex
%
% OUTPUT:
%   1.p_slack : a * b double
%
S_eq = abs(Y_eq .* S) - conj(Y_eq .* S);
p_slack = real(S_eq .* Y_eq) ./ (2 * real(S_eq));
end

function [lambda_without_p] = predict_without_p(Y_eq,S)
% Calculate predict result without  pmax
% where a and b are any natural number
% INPUT : 
%   1.Y_eq : a * b complex 
%   2.S : a * b complex
%
% OUTPUT:
%   1.lambda_without_p : a * b double
%
lambda_without_p = abs(Y_eq).*abs(Y_eq)./(2*(abs(Y_eq.*S)-real(Y_eq.*S)));
end

function [lambda_with_p,lambda_without_p] = predict_p_posi_final(Y_eq,S,p_posi)
% PREDICT_P_POSI_FINAL calculate the prediction with Pmax
% where a and b are any natural number
% INPUT : 
%   1.Y_eq : a * b complex 
%   2.S : a * b complex
%   3.p_neg : a * b double (For new energy part,p_neg <0)
%
% OUTPUT:
%   1.lambda_with_p : a * b double
%   2.lambda_without_p : a * b double
%
[lambda_without_p] = predict_without_p(Y_eq,S);
lambda_with_p = lambda_without_p;

ind_neg_real_S = intersect(find(real(S) < 0),find(real(Y_eq)>=0));
[slack_posi] = calculate_p_slack(Y_eq(ind_neg_real_S),...
    S(ind_neg_real_S)); 
ind_big_posi = find(p_posi(ind_neg_real_S) < slack_posi);
[lam_big_posi,~] = calculate_p_root(Y_eq(ind_neg_real_S(ind_big_posi )),...
    S(ind_neg_real_S(ind_big_posi )),...
    p_posi(ind_neg_real_S(ind_big_posi )));
lambda_with_p(ind_neg_real_S(ind_big_posi )) = lam_big_posi;
end

function [lam_big,lam_small] = calculate_p_root(Y_eq,S,P)
% CALCULATE_P_ROOT calculate the P root of Yeq, S and P
% INPUT : 
%   1.Y_eq : a * b complex 
%   2.S : a * b complex
%   3.P : a * b double (For new energy part,p_neg <0)
%
% OUTPUT:
%   1.lam_big : a * b double
%   2.lam_small : a * b double
%
S_p_big = real(S) - ...
    (abs(S).^2 - (imag(S) - 2 .* P .* imag(S./conj(Y_eq))).^2).^(0.5) ...
    + 2 .* P .* imag(S./conj(Y_eq)) .* j;  
lam_big = (2 * P .* (P - real(Y_eq))) ./ real(S_p_big .* Y_eq);
S_p_small = real(S) + ...
    (abs(S).^2 - (imag(S) - 2 .* P .* imag(S./conj(Y_eq))).^2).^(0.5) ...
    + 2 .* P .* imag(S./conj(Y_eq)) .* j;  
lam_small = (2 * P .* (P - real(Y_eq))) ./ real(S_p_small .* Y_eq);
end
