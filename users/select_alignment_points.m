% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
% Switzerland, Environmental Chemistry Modeling Laboratory, 2017. 
% Authors : Jonas Gros and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% This file is designed for the selection of alignment points on a GCxGC
% chromatogram. The positions of the alignment points will be saved to a
% file.
% The output file will be saved to a file called like the name of the input
% chromatogram file, with "_alignment_points" added to it.

% -----------------------------------------------------------------------------
% ALL PROGRAM PARAMETERS THAT SHOULD BE SET BY THE USER ARE SHOWN HERE.

% INSTRUMENT PARAMETERS

    % Set the modulation period (units of seconds)
    modulation_period = 12.5;

    % Set the detector acquisition rate (units of Hertz)
    acquisition_rate = 200; 
    
    % Set the acquisition delay (units of seconds)
    acquisition_delay = 0;
    
% INPUT/OUTPUT PARAMETERS
    
    % Set plot_flag to a value of 0 to suppress plots.
    % Set to a value of 1 to see "normal" level of plotting (DEFAULT).
    % Set to a value of 2 to see a lot of plots (for debugging/trouble-shooting).
    plot_flag = 1;

    % Set the output file path
    output_path = 'users/output/';

    % Set the input file path
    input_path = 'users/input/';

    % The file of the chromatogram on which you want to select alignment
    % points:
     chromatogram_file = 'Neat_oil_1000ppm+1ppmDDTP_1.peg_img01.csv';


    % Set Matlab console output level. Choose: 'minimal', 'normal', or 'verbose'.
    prompt_output = 'verbose';

% OPTIONAL:
    % If you only want to look for a set of alignment points (that you already
    % have) and try to guess their position in another chromatogram, set the
    % variable below to 1, and it will skip the selection of alignment points.
    % If doing so, you also need to provide the file with your known alignment
    % point positions ('alignment_points_to_guess_file').
    skip_alignment_point_selection = 0;
    if skip_alignment_point_selection
    %     Give the name of the file containing the position of the alignment
    %     points which positions you want to guess in other chromatograms (please,
    %     do also edit the parameters below).
        alignment_points_to_guess_file = 'Alignment_pts_Reference_OS0.csv';    
    end
    
    % If you want the alignment code to automatically try to identify the
    % corresponding peaks in another chromatogram, by taking the highest peak
    % in a window around the position of the peak (as in the first 
    % chromatogram) in another chromatogram. This may work only for moderate
    % shifts and will work better for high peaks with no neighboring higher
    % peaks. A check is strongly advised.
    % If you want to do this, fill the corresponding fields below.
    % If you want to use this, set the parameter below to 1, 0 otherwise:
    Look_other_chromatogram = 0;

    if Look_other_chromatogram
        %     Name of the additional chromatogram file(s):
        % (of type {'Chromatogram1.csv','Chromatogram2.csv','etc.csv'})
        other_chromatogram_file = {'J8_2X_split_1.peg_img011_norm.csv','J6_2x_1.peg_img01_norm.csv'};

        % Window size (first dimension and second dimension, in units of pixels):
        % (usually the default values ([2,50]) should be sufficient.)
        window_size = [2,50];

    end

% -----------------------------------------------------------------------------
% Do not modify the lines below.

cd('..');

addpath model_code

cd('model_code')

run main_code_pk_slct;

cd('../users');





