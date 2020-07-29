function [] = odour_specific_response_params(varargin)
% This function generates the parameters that are specific to the response
% of the entire PN population for a given odour. These parameters can then
% be used to generate multiple trials of responses to a given odour

    cd('/home/adithya/Documents/MATLAB/Odour_PNresponse_params_300PNs/');
    home = cd('/home/adithya/Documents/MATLAB/Odour_PNresponse_params_300PNs/');
    odour_folders = dir(home);
    number_already_existing_odours = length(odour_folders) - 2;


    % Default parameters. These can be changed upon calling the function via
    % varargin
    odour_ID = number_already_existing_odours + 1;
    number_PNs = 300; % number of Projection Neurons (PNs) in the system 
    trial_duration = 3000; % duration of the trial (ms)
    odour_onset_time = 500; %(ms)
    odour_duration = 1000; %(ms)
    fr_b_mean = 3.87 * trial_duration/1000; % fr_b_mean - mean of the normal distribution of baseline firing rates of PNs (spikes/second)
    fr_b_SD = 2.23; % fr_b_SD - SD of the normal distirbution of baseline firing rates of PNs (spikes/second)
    mean_num_PN_respond = round(0.2*number_PNs); % mean number of PNs that respond to an odour
    SD_num_PN_respond = round(0.05*number_PNs); % SD in number of PNs that respond to an odour


    if nargin > 1
        for jj = 1:2:length(varargin)-1
            if isa(varargin{jj+1},'double');
                eval([varargin{jj}, '=', num2str(varargin{jj+1})]);
            else
                error([varargin{jj} 'must be of type double']);
            end
        end
    end


    % Generating the parameters that define the odour's population PN response
    field1 = 'responding_PNs';
    field2 = 'n_r';
    field3 = 'n_e';
    field4 = 'd_f';
    field5 = 's_c';

    num_PN_respond = round(SD_num_PN_respond * randn(1) + mean_num_PN_respond);
    while num_PN_respond < 0 
        num_PN_respond = round(SD_num_PN_respond * randn(1) + mean_num_PN_respond);
    end    
    
    responding_PNs = zeros(length(num_PN_respond),1);
    
    for i = 1:num_PN_respond
       PN = randi(number_PNs,1);
        
       while ismember(PN,responding_PNs)==1
           PN = randi(number_PNs,1);
       end
       
       responding_PNs(i) = PN;
    end  
    
    n_r = [];
    n_e = [];
    d_f = [];
    s_c = [];

    for k = 1:num_PN_respond
        odour_spike_locations = [];
        n_r(1,k) = 10.67*randn(1) + 19.53; % n_r - number of spikes produced upon odour presentation
        while n_r(1,k) < 0 
            n_r(1,k) = 10.67*randn(1) + 19.53;
        end    
        n_e(1,k) = round(4*randn(1) + 8); % n_e - number of cycles of the LFP in which PN responds
        while n_e(1,k) <= 0 
            n_e(1,k) = round(4*randn(1) + 6);
        end
        d_f(1,k) = randi(odour_duration/50, 1); % d_f number of cycles from initiation of odour  when PN starts responding
        s_c(1,k) = round(n_r(1,k)/n_e(1,k)); % spikes per cycle in which PN responds
    end    

    parameters = struct(field1,responding_PNs,field2,n_r,field3,n_e,field4,d_f,field5,s_c);


    % Saving parameters to file
    filename = sprintf('Odour%d_params', odour_ID);
    foldername = sprintf('Odour%d',odour_ID);
    cd('/home/adithya/Documents/MATLAB/Odour_PNresponse_params_300PNs/');
    mkdir(foldername);
    new_foldername =  strcat('/home/adithya/Documents/MATLAB/Odour_PNresponse_params_300PNs/',foldername);
    cd(new_foldername)
    save(filename, 'parameters')

end