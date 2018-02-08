function [] = diffChromato(Chromato1,Chromato2,varargin) % alternative use: diffChromato(X,Y,Chromato1,Chromato2,varargin)

% Generates a pcolor plot of a difference (GCxGC) chromatogram. (It plots
% the difference chromatogram Chromato1 - Chromato2)
% 
% Inputs:
% - Chromato1 and Chromato2: only required input, two GCxGC chromatograms
% of size n x m (n for second dimension, m for first dimension).  
% Both chromatograms can either be a 2-D matrix, or a structure with
% compulsory field
% .data (the 2-D matrix), and optional fields .SR and .MP (sampling rate
% and modulation period, [Hz] and [s])
% Alternative use: diffChromato(X,Y,Chromato1,Chromato2,varargin), with X and Y
% matrices of the same size as Chromato, containing the positions on x-axis
% and y-axis. (similar to pcolor(X,Y,Z) versus pcolor(Z))
% - MP: modulation period, seconds. (alternative name is 'Mod')
% - SR: sampling rate, Hertz (= 1/seconds). (alternative name is 'Freq')
% - acquisition_delay: this is the shift for the first dimension first pixel
% (i.e. time is not starting at 0, but at acquisition_delay). Give it as seconds
% if MP and SR are provided, or in pixels otherwise.
% (alternative name of the option: 'FirstDshift')
% - RedBlue: if set to 1, uses a red-blue-white coloraxis. (default)
% - AutoAxis: if 1 (default value), will set automatically the color axis
% limits to minus (minimum) and plus (maximum) the values below which there
% are 97.5% of the
% data in the difference chromatogram. Color axis can always be adapted
% manually using matlab's caxis function.
% - title: a title for the figure. By default, this function display a
% title (name of the variable Chromato1 minus the name of the variable
% Chromato2). To remove this default title, set 'title' to 0.
% - ShadingMode (or simply 'shading'): default is 'interp', can alternatively be set to 'flat'
% (see pcolor help for explanation). (Note that all pixels are ploted if
% 'flat' is selected. (Default matlab [pcolor(Chromato); shading flat]
% would not plot the last pixels on the right & top; this is not the case
% here). Also, the x- and y-axis labels will fall in the middle of the 
% rectangle (versus bottom left with default matlab [pcolor(Chromato);
% shading flat].)
% - XAxisLocation: top or bottom (bottom is default)
% - AlkNBlabel: a column vector with 1st column carbon numbers, and second
% column retention times (unit is pixels/time) depending on chromatogram
% units. This replaces the normal x-axis in pixels or time units with an
% axis in carbon number units.
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
% Switzerland, Environmental Chemistry Modeling Laboratory, 2014. 
% Authors : Jonas Gros and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 
% File written by J.G., last modified on 10th July 2014.
% 

cnt = 1;
name1 = inputname(1); % name of variable 1
name2 = inputname(2); % name of variable 2
if ~isempty(varargin)
    if ~ischar(varargin{1})
        X = Chromato1;
        Y = Chromato2;
        Chromato1 = varargin{1};
        Chromato2 = varargin{2};
        name1 = inputname(3); % name of variable 1
        name2 = inputname(4); % name of variable 2
        cnt = 3;
    end
end

% First, if ever the user used some options, give their values to the
% corresponding variables:
for i=cnt:2:length(varargin)
    if or(strcmpi(varargin{i},'MP'),strcmpi(varargin{i},'Mod')) % modulation period [s]
        MP=varargin{i+1};
    elseif or(strcmpi(varargin{i},'SR'),strcmpi(varargin{i},'Freq')) % sampling rate, frequency [Hz]
        SR=varargin{i+1};
    elseif strcmpi(varargin{i},'RedBlue')
        RedBlue=varargin{i+1}; % 1 = yes, 0 = no
    elseif strcmpi(varargin{i},'AutoAxis')
        AutoAxis=varargin{i+1}; % 1 = yes, 0 = no
    elseif strcmpi(varargin{i},'title')
        titleFig=varargin{i+1};
    elseif or(strcmpi(varargin{i},'ShadingMode'),strcmpi(varargin{i},'shading'))
        ShadingMode=varargin{i+1};
    elseif strcmpi(varargin{i},'fontsize')
        FtSz=varargin{i+1};
    elseif strcmpi(varargin{i},'XAxisLocation')
        XAxisLocation=varargin{i+1};
    elseif strcmpi(varargin{i},'AlkNBlabel')
        AlkNBlabel=varargin{i+1};
    elseif or(strcmpi(varargin{i},'FirstDshift'),strcmpi(varargin{i},'acquisition_delay'))
        FirstDshift = varargin{i+1};
    else
        error(['The option ''' varargin{i} ''' is unknown'])
    end
end

if ~exist('RedBlue','var')
    RedBlue = 1;
end

if ~exist('MP','var')
    mp_does_not_exist = 1;
else
    mp_does_not_exist = 0;
end

% [nb_measY,nb_measX] = size(Chromato); % size for 1st and 2nd D columns.

if isstruct(Chromato1)
    if isfield(Chromato1,'MP')
        MP = Chromato1.MP;
        mp_does_not_exist_but_yes = 1;
    end
    if isfield(Chromato1,'SR')
        SR = Chromato1.SR;
    end
    Chromato1 = Chromato1.data;
end

if isstruct(Chromato2)
    if isfield(Chromato2,'MP')
        MP = Chromato2.MP;
        mp_does_not_exist_but_yes = 1;
    end
    if isfield(Chromato2,'SR')
        SR = Chromato2.SR;
    end
    Chromato2 = Chromato2.data;
end

% The code issues an error if Chromato1 and Chromato2 do not have the same size.
% In case they do not have the same size, let's here complete the smallest
% one with zeros to get the same size:
[Chromato1,Chromato2] = equalize_size(Chromato1,Chromato2);

if and(mp_does_not_exist==1,exist('mp_does_not_exist_but_yes','var'))
    if ~isempty(varargin)
        plotChromato(Chromato1-Chromato2,varargin{:},'MP',MP,'SR',SR)
    else
        plotChromato(Chromato1-Chromato2,'MP',MP,'SR',SR)
    end
else
    if ~isempty(varargin)
        plotChromato(Chromato1-Chromato2,varargin{:})
    else
        plotChromato(Chromato1-Chromato2)
    end
end

if RedBlue
    Coolormap=ones(64,3);
    Coolormap(1:31,2)=((1:31)')/31;
    Coolormap(1:31,3)=((1:31)')/31;
    Coolormap(34:64,1)=(flipud((1:31)'))/31;
    Coolormap(34:64,2)=(flipud((1:31)'))/31;
    colormap(Coolormap)
    
end
c_c = caxis;
caxis([-c_c(2),c_c(2)])

if ~exist('titleFig','var')
    if ~exist('FtSz','var')
        FtSz = 12; %
    end
    title(['\bf',name1,' - ',name2],'fontsize',FtSz)
elseif titleFig==0
    title('')
end
    
    
    
    
    
    

