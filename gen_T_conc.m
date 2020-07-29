function [T_conc_KC] = gen_T_conc( Connectivity_matrix , odour_ID, Trial_number )
% This function generates a matrix that represents the amount of
% neurotransmitters experienced by the population of KCs at their dendritic
% arbours. This conc. of neurotransmitters is summed over all PNs that 
% synapse onto a given KC. 
% The output matrix will be the transmitter conc. summed over all input for
% every KC, over the entire duration of the simulation.

% Reading in appropriate PN response files
%This first 'folder_name' is to be used when you are running normal main
%file

folder_name = sprintf('/home/adithya/Documents/MATLAB/Odour_PNresponse_params_300PNs/Odour%d',odour_ID);

% This folder_name should be used when running the V_thresh_oprimisation
% code
%folder_name = '/home/adithya/Documents/MATLAB/mult_trial';

file_name = sprintf('trial%d.csv',Trial_number);
cd(folder_name);
PN_response = csvread(file_name);
NID = length(PN_response(:,1));
TPT = length(PN_response(1,:));
% Calculating the T_conc matrix
% Parameters of the heaviside funciton used to convert spike time to T
A = 0.5 ;
t_max = 0.3 .* ones(300,1);

% Calculate the matrix of T conc of each PN
PN_T_conc = zeros(NID,TPT);

t0 = zeros(300,1);

for t = 1:TPT
    for k = 1:NID
        if PN_response(k,t) == 1
            t0(k,1) = t;
        end
    end    
    PN_T_conc(:,t) = A .* heaviside(t0 + t_max - (t .* ones(300,1))) .* heaviside((t.* ones(300,1))-t0);
end

% Calculate the matrix represnting the T conc viewed by every indivdual KC
KID = length(Connectivity_matrix(1,:));
PNID = length(Connectivity_matrix(:,1));
T_conc_KC = zeros(KID,TPT);

parfor KC_ID = 1:KID
    for PN_ID = 1:PNID
        T_conc_KC(KC_ID,:) = T_conc_KC(KC_ID,:) + PN_T_conc(Connectivity_matrix(PN_ID,KC_ID),:);
    end
end

end

