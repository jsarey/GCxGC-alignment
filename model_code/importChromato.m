
function [Chromato,file_type] = importChromato(filename,varargin) % modulation_period,sampling_rate)

% Function to import a GCxGC chromatogram file acquired with univariate
% detector (ECD, FID, etc.).
% 
% Inputs are:
% 
% REQUIRED:
% - filename: the name of the file to import (can include the whole path).
% (         [Chromato,file_type] = importChromato(filename)             )
% 
%   OPTIONAL (depending on file type):
% - modulation_period: the modulation period in seconds.
% - sampling_rate: the sampling rate in Herz (=1/seconds)
% (    [Chromato,file_type] = importChromato(filename,modulation_period,sampling_rate)      )
% 
% ALTERNATIVELY, the sampling rate and modulation period can be inputed
% using the options 'SR' and 'MP' ([Hz] and [s]).
% An additional option is 'struct'. If set to true (1), then the output
% variable 'Chromato' (the chromatogram) will be a structure with fields:
% Chromato.data = the 2-D chromatogram matrix (equivalent to the output of
% the function when this option is not used). In addition, two fields:
% Chromato.MP and Chromato.SR will contain modulation period and sampling
% rate ([s] and [Hz]). May be used to avoid needing to re-specify the MP
% and SR later.
% 
% Outputs are:
% Chromato: a chromatogram, as a 2-D matrix (of size 2nd dimension size
% times 1st dimension size).
% file_type: Optional output, that tells the type of file assumed (GC Image
% csv file, ChromaTOF csv file, Bob Nelson-way txt file, mat file).
% 
% In any case, if the file does not contain a number of elements that is a
% multiple of modulation_period*sampling_rate, then the last elements are
% discharded.
% 
% .csv files are assumed as either coming from ChromaTOF or GC Image.
% ChromaTOF files are identified based on the first line that contain more
% than one heading, separated with "," or ";". The column labeled "S1" or
% "TIC" is imported.
% .csv file not from ChromaTOF are assumed from GC Image (single column
% column vector or 2-D matrix).
% .txt files are assumed of the "Bob Nelson type".
% .mat files are assumed 2-D matrices.
% 
% File written by J.G., last modified on 1st July 2015.
% 
% (Corrected for Bob's txt file to remove the row and the column with the
% pixel numbers. (MArch 2015))
% 1st July 2015: corrected to allow import of ChromaTOF data when the
% column with data is labeled TIC in addition to data when this column is
% labeled S1.
% 
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE,
% Switzerland, Environmental Chemistry Modeling Laboratory, 2014.
% Author : Jonas Gros.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if length(varargin)==2
    modulation_period = varargin{1};
    sampling_rate = varargin{2};
elseif length(varargin)>2

    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'MP') % modulation period [s]
            modulation_period=varargin{i+1};
        elseif strcmpi(varargin{i},'SR') % sampling rate, frequency [Hz]
            sampling_rate=varargin{i+1};
        elseif strcmpi(varargin{i},'struct') % structure output, 1 = true, 0 = false
            struct_output=varargin{i+1};
        else
            error(['The option ''' varargin{i} ''' is unknown'])
        end
    end
end

if ~exist('struct_output','var')
    struct_output = 0;
end

if strcmpi(filename(end-3:end),'.csv')

    fileID = fopen(filename); % open the file
    
    A = fgetl(fileID); % A is a string (first line of the file)

    if or(sum(findstr(A,';')),sum(findstr(A,','))) % ChromaTOF File or GC Image 2-D file
        if findstr(A,';')
            Separator = ';'; % Separator can change depending on computer 
    %         language (French versus English)
        else
            Separator = ','; %
        end
        
        if findstr(A,'S1') % THe data in ChromaTOF is either in the column
            column_heading = 'S1'; % labeled "S1+ or "TIC", apparently.
        elseif findstr(A,'TIC')
            column_heading = 'TIC';
        else
            column_heading = 'S1'; % WE DON'T FIND IT ANYWAY.
        end
        
        if findstr(A,column_heading) % ChromaTOF file
            Column = sum(findstr(A,Separator)<findstr(A,column_heading));  % Usually 3 as 
        %     I experienced
            if or(~exist('modulation_period','var'),~exist('sampling_rate','var'))
                    error('PLEASE, do specify modulation period (MP) and sampling rate (SR).')
            end
            Chromato_nothing = dlmread(filename,Separator,1,Column); % the 1 skips the first line, which 
            % contains headings in ChromaTOF files.
            length_data = size(Chromato_nothing,1);
            clear Chromato_nothing
%             Line below modified to avoir some error:
            Chromato = dlmread(filename,Separator,[1 Column length_data Column]); % the 1 skips the first line, which 
            % contains headings in ChromaTOF files.
%             Chromato = dlmread(filename,Separator,1,Column); % the 1 skips the first line, which 
            % contains headings in ChromaTOF files.
            Chromato = reshape(Chromato(1:end-mod(length(Chromato),...
                modulation_period*sampling_rate)),...
                modulation_period*sampling_rate,[]); % Reshape the chromatogram.
    %         Uncomplete line at the end is neglected.
            file_type = 'ChromaTOF'; % give the file type
        else % GC Image 2-D file
            Chromato = flipud(importdata(filename));
            Chromato = Chromato(1:end-1,:); % There was a last line of NaNs.
            file_type = 'GC Image 2-D'; % give the file type
        end
    else % assume a GC image 1-D file, one big vector:
        if or(~exist('modulation_period','var'),~exist('sampling_rate','var'))
                error('PLEASE, do specify modulation period (MP) and sampling rate (SR).')
        end
        Chromato = importdata(filename);
        Chromato = reshape(Chromato(1:end-mod(length(Chromato),...
            modulation_period*sampling_rate)),...
            modulation_period*sampling_rate,[]); % Reshape the chromatogram.
%         Uncomplete line at the end is neglected.
        file_type = 'GC Image 1-D'; % give the file type
    end

    fclose(fileID); % close the file
    
elseif strcmpi(filename(end-3:end),'.mat') % mat file, assumed 2-D matrix
    Chromato = importdata(filename);
    file_type = '.mat (assumed 2-D)'; % give the file type
    
elseif strcmpi(filename(end-3:end),'.txt')
%     Let's assume a file format like the files sent to Jonas by Bob
%     Nelson.
    Chromato = flipud(importdata(filename));
    Chromato = Chromato(1:end-1,2:end);
    file_type = 'Bob Nelson txt file'; % give the file type
    
end

if struct_output
    warning('off')
    Chromato.data = Chromato;
    warning('on')
    if exist('modulation_period','var')
        Chromato.MP = modulation_period;
    end
    if exist('sampling_rate','var')
        Chromato.SR = sampling_rate;
    end
end

