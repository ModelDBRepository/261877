% This function generates a 900 X 1 vector of length PNnumber. This vector 
% consists of 0s and 1s.This serves as the input for a time-independent 
% version of the KC model. Basically this vector is a representation of PN 
% activity at a single time point /window.

function [] = gen_singleTP_PNact(number_of_reps, PNnumber,p_of_spiking)
% p_of_spiking tells you the probability that a neuron will spike in this
% time window.
% number_of_reps tells you the number of such PN representations to create.

for j = 1:number_of_reps

    PN_activity = zeros(PNnumber,1);
    
    for i = 1:PNnumber
        % Generate the activity by drawing values from a uniform random
        % distribtuion and using p_of_spiking as a threshold to decide if
        % PN i spikes or not.
        if rand(1) < p_of_spiking
            PN_activity(i,1) = 1;
        end
    end   
    
    % Save PN Activity as .csv file
    filename = sprintf('/home/adithya/Documents/MATLAB/single_time_window_PNrep/single_cycle/trial%d',j);
    xlswrite(filename,PN_activity);
end
