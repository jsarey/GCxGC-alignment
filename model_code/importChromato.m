
function [Chromato,file_type] = importChromato(filename,varargin) % varargin might be equal to: modulation_period,sampling_rate)

% Function to import a GCxGC chromatogram file acquired with a univariate
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
% (    [Chromato,file_type] = importChromato(filename,modulation_period,sampling_rate)      
% - file_type: 'ChromaTOF','GC Image 1-D','GC Image 2-D','mat',
% 'Bob Nelson txt', 'Karin Lemkau txt' ---> If the automatic detection of
% file type fails, to force the choice. (Bob Nelson txt = including a line
% and a column with retention times, Karin Lemkau txt = without that).
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
% discarded.
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
% File written by J.G., last modified on 20st June 2017.
% 
% (Added the ability to import a second, slightly different Bob-Nelson
% format. (June 2017))
% (Corrected for Bob's txt file to remove the row and the column with the
% pixel numbers. (March 2015))
% 1st July 2015: corrected to allow import of ChromaTOF data when the
% column with data is labeled TIC in addition to data when this column is
% labeled S1.
% 
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE,
% Switzerland, Environmental Chemistry Modeling Laboratory, 2016.
% Author : Jonas Gros.
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
        elseif strcmpi(varargin{i},'file_type') % file type, string
            file_type=varargin{i+1};
        else
            error(['The option ''' varargin{i} ''' is unknown'])
        end
    end
end

if ~exist('struct_output','var')
    struct_output = 0;
end

if ~exist('file_type','var')
    file_type = 'auto'; % automatically guess the file type.
end

if or(and(strcmpi(file_type,'auto'),strcmpi(filename(end-3:end),'.csv')),...
        or(strcmpi(file_type,'GC Image 1-D'),or(strcmpi(file_type,'GC Image 2-D'),...
        strcmpi(file_type,'ChromaTOF'))))

    fileID = fopen(filename); % open the file
    
    A = fgetl(fileID); % A is a string (first line of the file)

    if or(and(strcmpi(file_type,'auto'),or(sum(findstr(A,';')),sum(findstr(A,',')))),...
          or(strcmpi(file_type,'GC Image 2-D'),...
        strcmpi(file_type,'ChromaTOF')))  % ChromaTOF File or GC Image 2-D file
        if findstr(A,';')
            Separator = ';'; % Separator can change depending on computer 
    %         language (French versus English)
        else
            Separator = ','; %
        end
        
        if findstr(A,'S1') % The data in ChromaTOF is either in the column
            column_heading = 'S1'; % labeled "S1+ or "TIC", apparently.
        elseif findstr(A,'TIC')
            column_heading = 'TIC';
        else
            column_heading = 'S1'; % WE DON'T FIND IT ANYWAY.
        end
        
        if or(and(strcmpi(file_type,'auto'),findstr(A,column_heading)),...
               strcmpi(file_type,'ChromaTOF')) % ChromaTOF file
            Column = sum(findstr(A,Separator)<findstr(A,column_heading));  % Usually 3 as 
        %     I experienced
            if or(~exist('modulation_period','var'),~exist('sampling_rate','var'))
                    error('Please, do specify modulation period (MP) and sampling rate (SR).')
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
        elseif or(strcmpi(file_type,'auto'),strcmpi(file_type,'GC Image 2-D')) % GC Image 2-D file
            Chromato = flipud(importdata(filename));
            Chromato = Chromato(1:end-1,:); % There was a last line of NaNs.
            file_type = 'GC Image 2-D'; % give the file type
        end
    elseif or(strcmpi(file_type,'auto'),strcmpi(file_type,'GC Image 1-D')) % assume a GC image 1-D file, one big vector:
        if or(~exist('modulation_period','var'),~exist('sampling_rate','var'))
                error('Please, do specify modulation period (MP) and sampling rate (SR).')
        end
        Chromato = importdata(filename);
        Chromato = reshape(Chromato(1:end-mod(length(Chromato),...
            modulation_period*sampling_rate)),...
            modulation_period*sampling_rate,[]); % Reshape the chromatogram.
%         Uncomplete line at the end is neglected.
        file_type = 'GC Image 1-D'; % give the file type
    end

    fclose(fileID); % close the file
    
elseif or(and(strcmpi(file_type,'auto'),strcmpi(filename(end-3:end),'.mat')),...
        strcmpi(file_type,'mat')) % mat file, assumed 2-D matrix
    Chromato = importdata(filename);
    file_type = '.mat (assumed 2-D)'; % give the file type
    
elseif or(and(strcmpi(file_type,'auto'),strcmpi(filename(end-3:end),'.txt')),...
        strcmpi(file_type,'Bob Nelson txt'))
%     Let's assume a file format like the files sent to Jonas by Bob
%     Nelson.
    Chromato = flipud(importdata(filename));
    if ~isstruct(Chromato) % Apparently Bob Nelson has two slightly different file formats
        Chromato = Chromato(1:end-1,2:end);
%         Column_headers = Chromato(:,1);
%         Row_headers = Chromato(1,:);
        file_type = 'Bob Nelson txt file'; % give the file type
    else % different file format received by Charlie Sharpless from Bob Nelson ('Sharpless_012115_0h.txt')
        Column_headers = zeros(length(Chromato.rowheaders)-1,1);
        for tutut = 1:length(Column_headers)
            % (We discard the first value, that I don't understand)
            Column_headers(tutut) = str2double(Chromato.rowheaders{tutut + 1});    
        end
        Row_headers = Chromato.data(1,:);
        Chromato = Chromato.data(2:end,1:end);
        if Column_headers(1)>Column_headers(2)
            Chromato = flipud(Chromato);
        end
        file_type = 'Bob Nelson type II txt file'; % give the file type
    end
elseif strcmpi(file_type,'Karin Lemkau txt')
%     Used only if chosen by the user for now. Could make the choice
%     between Bob Nelson type and Lemkau type automatic by considering the
%     regularity of values on the line and columns. Not implemented yet.
    Chromato = importdata(filename);
    Chromato = Chromato(1:end,1:end);
    file_type = 'Karin Lemkau txt file'; % give the file type
    
end

if struct_output
    Chromato_temp = Chromato;
    clear Chromato
    Chromato.data = Chromato_temp;
    clear Chromato_temp
    if exist('modulation_period','var')
        Chromato.MP = modulation_period;
    end
    if exist('sampling_rate','var')
        Chromato.SR = sampling_rate;
    end
end

