%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: innovation_output.m
% Author: Nicholas von Turkovich
% Date: 
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = innovation_output(econparams_bgp, state_space, event_space, dV_leader, dV_follower)
    
    % create value matrices for leader and follower given state space
    val_follower = econparams_bgp.v_vec(econparams_bgp.s_max + 1 - state_space);
    val_leader = econparams_bgp.v_vec(econparams_bgp.s_max + 1 + state_space);
    
    % leader
    leader_innovations = (event_space == 1);
    newval_leader = val_leader.*leader_innovations;
    oldval_leader = val_leader.*([leader_innovations(:,(2:end)), zeros(size(event_space, 1), 1)]);
    dV_leader = [zeros(size(event_space, 1),1), newval_leader(:,2:end) - oldval_leader(:,1:(end-1))];
    dV_ann_leader = reshape(dV_leader, size(dV_leader,1), 12, size(dV_leader,2)/12);
    dV_ann_leader = reshape(sum(dV_ann_leader,2), size(dV_ann_leader, 1), size(dV_ann_leader,3));
    val_ann_leader = reshape(val_leader, size(val_leader, 1), 12, size(val_leader,2)/12);
    val_ann_leader = reshape(val_ann_leader(:,end,:), size(val_ann_leader,1), size(val_ann_leader,3));
    pct_leader = dV_ann_leader./val_ann_leader;
    
    % find leader innovation events
    follower_innovations = (event_space == 2);
    newval_follower = val_follower.*follower_innovations;
    oldval_follower = val_follower.*([follower_innovations(:,(2:end)), zeros(size(event_space, 1), 1)]);
    dV_follower = [zeros(size(event_space, 1),1), newval_follower(:,2:end) - oldval_follower(:,1:(end-1))];
    
    

    

end

