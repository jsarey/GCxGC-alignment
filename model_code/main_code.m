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

% addpath model_code

% Import chromatograms:

disp(' ')
disp(' ')
disp('---------------------------------------------------------')
disp(' ')
disp('Please cite the following article when publishing any results')
disp('obtained by use of this software:')
disp(' ')
disp('Gros, J., Nabi, D., Dimitriou-Christidis, P., Rutler, R., Arey, J. S.,')
disp('“Robust algorithm for aligning two-dimensional chromatograms”,')
disp('Anal. Chem. 2012, vol 84, p 9033-9040.')
disp('<a href = "http://pubs.acs.org/doi/abs/10.1021/ac301367s">(hyperlink)</a>');
disp(' ')
disp(' ')
disp('---------------------------------------------------------')
disp(' ')
disp(' ')

% A few error checks:
if ~exist(['../',input_path,Target_chromatogram_file],'file')
    error(['Target chromatogram file (',['../users/input/',Weathered_chromatogram_file],' ) does not exist.'])
end
if ~exist(['../',input_path,Reference_chromatogram_file],'file')
    error(['Reference chromatogram file (',['../users/input/',Reference_chromatogram_file],' ) does not exist.'])
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

[Ref,Ref_file_type] = importChromato([strrep(pwd,'model_code',''),input_path,...
    Reference_chromatogram_file],...
    'SR',acquisition_rate,'MP',modulation_period,'struct',1);

if strcmpi(prompt_output,'verbose')
    disp(['Reference chromatogram file assumed of the type:   ',Ref_file_type])
    disp(['( ',Reference_chromatogram_file,' )'])
    disp(' ')
end

[Target,Target_file_type] = importChromato([strrep(pwd,'model_code',''),input_path,...
    Target_chromatogram_file],...
    'SR',acquisition_rate,'MP',modulation_period,'struct',1);

if strcmpi(prompt_output,'verbose')
    disp(['Target chromatogram file assumed of the type:   ',Target_file_type])
    disp(['( ',Target_chromatogram_file,' )'])
    disp(' ')
end

% Import the positions of the alignment points:
% (Including a check whether the file do actually exist)
if ~exist([strrep(pwd,'model_code',''),input_path,...
    Reference_alignment_pts_file],'file')
    if exist([strrep(pwd,'model_code',''),output_path,...
    Reference_alignment_pts_file],'file')
        warning(['The file containing reference alignment points',...
            ' was found in the output_path and not input_path.',...
            ' Loading the file in output_path...'])
        Reference_peaks = importdata([strrep(pwd,'model_code',''),output_path,...
        Reference_alignment_pts_file]);
    else
        error(['The file containing reference alignment points (',...
            [strrep(pwd,'model_code',''),input_path,...
        Reference_alignment_pts_file],' ) does not exist.'])
    end
else
    Reference_peaks = importdata([strrep(pwd,'model_code',''),input_path,...
        Reference_alignment_pts_file]);
end

if ~exist([strrep(pwd,'model_code',''),input_path,...
    Target_alignment_pts_file],'file')
    if exist([strrep(pwd,'model_code',''),output_path,...
    Target_alignment_pts_file],'file')
        warning(['The file containing target alignment points',...
            ' was found in the output_path and not input_path.',...
            ' Loading the file in output_path...'])
        Target_peaks = importdata([strrep(pwd,'model_code',''),output_path,...
        Target_alignment_pts_file]);
    else
        error(['The file containing target alignment points (',...
            [strrep(pwd,'model_code',''),input_path,...
        Target_alignment_pts_file],' ) does not exist.'])
    end
else
    Target_peaks = importdata([strrep(pwd,'model_code',''),input_path,...
        Target_alignment_pts_file]);
end

Aligned = alignChromato(Ref.data,Target.data,Reference_peaks,Target_peaks,...
    'Peak_widths',typical_peak_width);

Aligned_temp = Aligned;
clear Aligned
Aligned.data = Aligned_temp;
clear Aligned_temp
Aligned.SR = Ref.SR;
Aligned.MP = Ref.MP;

if plot_flag 
    set(figure,'name','GCxGC chromatograms')
    subplot(3,1,1); plotChromato(Ref,'acquisition_delay',acquisition_delay); c_a_xis = caxis; title('\bfReference'); ylabel(''); xlabel('')
    hold on; ALgnt_pks_Time = Pix2Time([Reference_peaks(:,1),Reference_peaks(:,2)],modulation_period,acquisition_rate);
    plot(ALgnt_pks_Time(:,1),ALgnt_pks_Time(:,2),'ok','linewidth',2)
    
    subplot(3,1,2); plotChromato(Target,'acquisition_delay',acquisition_delay); caxis(c_a_xis); title('\bfTarget'); xlabel('')
    hold on; ALgnt_pks_Time2 = Pix2Time([Target_peaks(:,1),Target_peaks(:,2)],modulation_period,acquisition_rate);
    plot(ALgnt_pks_Time2(:,1),ALgnt_pks_Time2(:,2),'ok','linewidth',2)
    
    subplot(3,1,3); plotChromato(Aligned,'acquisition_delay',acquisition_delay); caxis(c_a_xis); title('\bfAligned'); ylabel('')
    hold on; 
    plot(ALgnt_pks_Time(:,1),ALgnt_pks_Time(:,2),'ok','linewidth',2)
    
    set(figure,'name','Difference chromatograms')
    subplot(3,1,2); diffChromato(Ref,Target,'acquisition_delay',acquisition_delay); c_a_xis = caxis; title('\bfReference - Target'); xlabel('');
    
    subplot(3,1,3); diffChromato(Ref,Aligned,'acquisition_delay',acquisition_delay); caxis(c_a_xis); title('\bfReference - Aligned'); ylabel('')

end

if strcmpi(prompt_output,'minimal')
    save_flag = 1;
else
    disp('Save aligned chromatogram? (1 = yes, 0 = no)')
    save_flag = input('');
end

if save_flag
    dot_indices = strfind(Target_chromatogram_file,'.');
    Output_file_name = [strrep(pwd,'model_code',''),output_path,...
        Target_chromatogram_file(1:dot_indices(end)-1),'_ALIGNED.csv'];
    To_Save = Aligned.data(:);
    
    if length(To_Save)>50000
        tic
        dlmwrite(Output_file_name,To_Save(1:10000),'precision','%20.20f')
        zu = toc; % seconds (to write 1e4 pixels to the file)
        Speed = zu/10000; % seconds per pixel.
        disp(['(estimated time to completion: ',...
        num2str(length(To_Save)*Speed),' seconds)'])
        dlmwrite(Output_file_name,To_Save(10001:end),'-append','precision','%20.20f')
    else
        disp('(estimated time to completion: fast.)')
        dlmwrite(Output_file_name,To_Save,'precision','%20.20f')
    end

    zu = toc;
    disp(' ')
    disp(['Done in ',num2str(zu),' seconds.'])
    disp(' ')
end


