
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
% Switzerland, Environmental Chemistry Modeling Laboratory, 2014. 
% Authors : Jonas Gros and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% -----------------------------------------------------------------------------
% ALL PROGRAM PARAMETERS THAT SHOULD BE SET BY THE USER ARE SHOWN HERE.

% INSTRUMENT PARAMETERS

    % Set the modulation period (units of seconds)
    modulation_period = 12.5;

    % Set the detector acquisition rate (units of Hertz)
    acquisition_rate = 200; 

% MODEL CHOICE PARAMETERS

    % What is the typical width of a peak in first and second dimension
    % (in units of pixels) 
    typical_peak_width = [3 120];

% INPUT/OUTPUT PARAMETERS

    % Set plot_flag to a value of 0 to suppress plots.
    % Set to a value of 1 to see "normal" level of plotting (DEFAULT).
    plot_flag = 1;

    % Set the output file path
    output_path = 'users/output/';

    % Set the input file path
    input_path = 'users/input/';

    % Name of input-output files:
        % Reference chromatogram:
        Reference_chromatogram_file = 'Neat_oil_1000ppm+1ppmDDTP_1.peg_img01.csv';

        % Target chromatogram:
        Target_chromatogram_file = 'J8_2X_split_1.peg_img011_norm.csv';

        % Positions of alignment points in the Reference chromatogram:
        Reference_alignment_pts_file = 'Alignment_pts_Reference_OS0.csv';

        % Positions of alignment points in the Target chromatogram:
        Target_alignment_pts_file = 'Alignment_pts_Target_OS2.csv';

    % Set Matlab console output level. Choose: 'minimal', 'normal', or 'verbose'.
    prompt_output = 'normal';

% -----------------------------------------------------------------------------
% Do not modify the lines below.

cd('..');

addpath model_code

cd('model_code')

run main_code;

cd('../users');



