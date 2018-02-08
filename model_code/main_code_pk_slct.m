% main_code.m script called from "main.m" by user.
% 
% *** Do not modify this file. Normally the user should not need to adjust *** 
% *** anything in this script.                                             ***
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
% Switzerland, Environmental Chemistry Modeling Laboratory, 2014. 
% Authors : Jonas Gros and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

addpath ../model_code/

% A few error checks:
if ~exist(['../',input_path,chromatogram_file],'file')
    error(['Chromatogram file (',['../users/input/',chromatogram_file],' ) does not exist.'])
end
if ~exist(['../',output_path],'dir')
    warning(['Indicated output directory, ''',...
        output_path,''' does not exist... Creating it.'])
    mkdir(['../',output_path])
end
% Check if allowed to write to the output_path:
[stat,CcC] = fileattrib(['../',output_path]); 
if ~CcC.UserWrite
    error(['Matlab is not allowed to write to the output_path folder...',...
        ' Try displacing the folder! (e.g. to ''Desktop'' on a Windows computer)'])
end

% Import chromatogram:
[Chromato,file_type] = importChromato([strrep(pwd,'model_code',''),input_path,chromatogram_file],...
    'SR',acquisition_rate,'MP',modulation_period);

if strcmpi(prompt_output,'verbose')
    disp(['Chromatogram file assumed of the type:   ',file_type])
    disp(['( ',chromatogram_file,' )'])
    disp(' ')
end

if ~skip_alignment_point_selection
    global Peak_Pos_c % Might be useful in case of erros. (Peak_Pos_c
%     contains the positions of the peaks already selected, in case an
%     error occurs during the peak selection.

    % Ask the user to select alignment points:
    if strcmpi(prompt_output,'minimal')
        Peak_Pos = give_peak_pos(Chromato,'EXPL',0);
    else
        Peak_Pos = give_peak_pos(Chromato,'EXPL',1);
    end

    

    dot_indices = strfind(chromatogram_file,'.');
    Output_file_name = [strrep(pwd,'model_code',''),output_path,...
        chromatogram_file(1:dot_indices(end)-1),'_alignment_points.csv'];
    % Save to file:
    if or(strcmpi(prompt_output,'verbose'),strcmpi(prompt_output,'normal'))
        disp('Saving selected peaks to file...')
        if strcmpi(prompt_output,'verbose')
            disp(Output_file_name)
            disp(' ')
        end
    end
    csvwrite(Output_file_name,Peak_Pos)
    
elseif skip_alignment_point_selection
%     Import the positions of alignment points (already known, no manual
%     selection needed):
    Peak_Pos = importdata([strrep(pwd,'model_code',''),input_path,...
    alignment_points_to_guess_file]);
    if ~(strcmpi(prompt_output,'minimal'))
        disp('Trying to match the alignment points in file:')
        disp(alignment_points_to_guess_file)
        disp('....')
    end
    if plot_flag
        figure; plotChromato(Chromato,'shading','flat');
        hold on; plot(Peak_Pos(:,1),Peak_Pos(:,2),'ok','linewidth',2)
        title('\bfPositions of the peaks in the original chromatogram','fontsize',18)
    end
end
    
if Look_other_chromatogram
%     If the user changed the format of the variable by mistake:
    if ~iscell(other_chromatogram_file)
        other_chromatogram_file = {other_chromatogram_file};
    end
%     Iterate over each file:
    for k = 1:length(other_chromatogram_file)
        
        % Import chromatogram:
        [Other,file_type] = importChromato([strrep(pwd,'model_code',''),input_path,other_chromatogram_file{k}],...
            'SR',acquisition_rate,'MP',modulation_period);

        if strcmpi(prompt_output,'verbose')
            disp(['Chromatogram file assumed of the type:   ',file_type])
            disp(['( ',other_chromatogram_file{k},' )'])
            disp(' ')
        end
        
%         Guess peak positions:
        Peaks_Other_est = Find_Peaks_Other(Other,Peak_Pos,plot_flag,...
            window_size(1),window_size(2),1,min(5,window_size(2)));
        
        disp('Proceed to saving? (1 = yes, 0 = no)')
        saving_flag = 'nothing yet';
        while ~or(strcmpi(saving_flag,'0'),strcmpi(saving_flag,'1'))
            saving_flag = input('','s');
        end
        if strcmpi(saving_flag,'1')
            current_file = other_chromatogram_file{k};
            dot_indices = strfind(current_file,'.');
            Output_file_name = [strrep(pwd,'model_code',''),output_path,...
                current_file(1:dot_indices(end)-1),'_alignment_points.csv'];
            csvwrite(Output_file_name,Peaks_Other_est)
            
            if strcmpi(prompt_output,'verbose')
                disp('Saving guessed peaks to file...')
                disp(['( ',Output_file_name,' )'])
                disp(' ')
            end
        end
        clear Peaks_Other_est
    end
    
end

