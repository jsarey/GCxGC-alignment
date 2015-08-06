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
Reference_peaks = importdata([strrep(pwd,'model_code',''),input_path,...
    Reference_alignment_pts_file]);

Target_peaks = importdata([strrep(pwd,'model_code',''),input_path,...
    Target_alignment_pts_file]);

Aligned = alignChromato(Ref.data,Target.data,Reference_peaks,Target_peaks,...
    'Peak_widths',typical_peak_width);

warning('off')
Aligned.data = Aligned;
warning('on')
Aligned.SR = Ref.SR;
Aligned.MP = Ref.MP;

if plot_flag 
    set(figure,'name','GCxGC chromatograms')
    subplot(3,1,1); plotChromato(Ref); c_a_xis = caxis; title('\bfReference'); ylabel(''); xlabel('')
    hold on; ALgnt_pks_Time = Pix2Time([Reference_peaks(:,1),Reference_peaks(:,2)],modulation_period,acquisition_rate);
    plot(ALgnt_pks_Time(:,1),ALgnt_pks_Time(:,2),'ok','linewidth',2)
    
    subplot(3,1,2); plotChromato(Target); caxis(c_a_xis); title('\bfTarget'); xlabel('')
    hold on; ALgnt_pks_Time2 = Pix2Time([Target_peaks(:,1),Target_peaks(:,2)],modulation_period,acquisition_rate);
    plot(ALgnt_pks_Time2(:,1),ALgnt_pks_Time2(:,2),'ok','linewidth',2)
    
    subplot(3,1,3); plotChromato(Aligned); caxis(c_a_xis); title('\bfAligned'); ylabel('')
    hold on; 
    plot(ALgnt_pks_Time(:,1),ALgnt_pks_Time(:,2),'ok','linewidth',2)
    
    set(figure,'name','Difference chromatograms')
    subplot(3,1,2); diffChromato(Ref,Target); c_a_xis = caxis; title('\bfReference - Target'); xlabel('');
    
    subplot(3,1,3); diffChromato(Ref,Aligned); caxis(c_a_xis); title('\bfReference - Aligned'); ylabel('')

end

if strcmpi(prompt_output,'minimal')
    save_flag = 1;
else
    disp('Save aligned chromatogram? (1 = yes, 0 = no)')
    save_flag = input('');
end

if save_flag
    csvwrite(strrep([strrep(pwd,'model_code',''),output_path,...
    Target_chromatogram_file],'.csv','_ALIGNED.csv'),Aligned.data(:))
    
end


