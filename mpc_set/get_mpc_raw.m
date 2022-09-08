%% For convenience of add line
%% Only contain the mpc part of 
%% version, baseMVA, bus, gen, branch
function [mpc_raw] = get_mpc_raw(mpc)
%% First ext2int
mpc_temp = ext2int(mpc);
%% Only contain 5 part
mpc_raw.version = mpc_temp.version;
mpc_raw.baseMVA = mpc_temp.baseMVA;
mpc_raw.bus = mpc_temp.bus;
mpc_raw.gen = mpc_temp.gen;
mpc_raw.branch = mpc_temp.branch;
end