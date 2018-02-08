function [] = plotChromato(Chromato,varargin) % alternative use: plotChromato(X,Y,Chromato,varargin)

% Generates a pcolor plot of a GCxGC chromatogram.
% 
% Inputs:
% - Chromato: only required input, a GCxGC chromatogram of n x m (n for
% second dimension, m for first dimension).
% Chromato can either be a 2-D matrix, or a structure with compulsory field
% .data (the 2-D matrix), and optional fields .SR and .MP (sampling rate
% and modulation period, [Hz] and [s])
% Alternative use: plotChromato(X,Y,Chromato,varargin), with X and Y
% matrices of the same size as Chromato, containing the positions on x-axis
% and y-axis. (similar to pcolor(X,Y,Z) versus pcolor(Z))
% - MP: modulation period, seconds. (alternative name is 'Mod')
% - SR: sampling rate, Hertz (= 1/seconds). (alternative name is 'Freq')
% - acquisition_delay: this is the shift for the first dimension first pixel
% (i.e. time is not starting at 0, but at acquisition_delay). Give it as seconds
% if MP and SR are provided, or in pixels otherwise.
% - RedBlue: if set to 1, uses a red-blue coloraxis as used for the Gros et
% al. 2014 biodegradation paper.
% - colormap: can be either of matlab's defaults, or the 'RedBlue' defined
% above. Setting colormap to 'RedBlue' has the same effect than setting
% 'RedBlue' to 1. 'RedBlue2' has a lighter blue tone.
% - AutoAxis: if 1 (default value), will set automatically the color axis
% limits to 0 (minimum) and the values below which there are 97.5% of the
% data (maximum). Color axis can always be adapted manually using matlab's
% caxis function.
% - title: a title for the figure
% - ShadingMode (or simply 'shading'): default is 'interp', can alternatively be set to 'flat'
% (see pcolor help for explanation). (Note that all pixels are ploted if
% 'flat' is selected. (Default matlab [pcolor(Chromato); shading flat]
% would not plot the last pixels on the right & top; this is not the case
% here). Also, the x- and y-axis labels will fall in the middle of the 
% rectangle (versus bottom left with default matlab [pcolor(Chromato);
% shading flat].)
% - XAxisLocation: top or bottom (bottom is default)
% - AlkNBlabel: a column vector with 1st column carbon numbers, and second
% column retention times (unit is pixels/time depending on chromatogram
% units). This replaces the normal x-axis in pixels or time units with an
% axis in carbon number units.
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Authors : Jonas Gros and J. Samuel Arey.
% See license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 
% File written by J.G., last modified on 12th June 2017.
% 


cnt = 1;
if ~isempty(varargin)
    if ~ischar(varargin{1})
        X = Chromato;
        Y = varargin{1};
        Chromato = varargin{2};
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
    elseif strcmpi(varargin{i},'colormap')
        if strcmpi(varargin{i+1},'RedBlue')
            RedBlue=1; % 1 = yes, 0 = no
        elseif strcmpi(varargin{i+1},'RedBlue2')
            RedBlue=19; % 1 = yes, 0 = no
        else
            color_map = varargin{i+1};
        end
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
    elseif strcmpi(varargin{i},'acquisition_delay')
        acquisition_delay = varargin{i+1}/60; % converts to units of minutes
    else
        error(['The option ''' varargin{i} ''' is unknown'])
    end
end

if ~exist('RedBlue','var')
    RedBlue = 0;
    if ~exist('color_map','var')
        color_map = 'jet'; % inforce a default when not provided.
    end
end

if ~exist('ShadingMode','var')
    ShadingMode = 'interp'; % default
end

if ~exist('AutoAxis','var')
    AutoAxis = 1; % by default, enable auto-axis
end

if ~exist('FtSz','var')
    FtSz = 12; %
end

if ~exist('XAxisLocation','var')
    XAxisLocation = 'bottom'; % by default, x axis at bottom.
end

if ~exist('acquisition_delay','var')
    acquisition_delay = 0; % by default, no shift.
end

if isstruct(Chromato)
    if isfield(Chromato,'MP')
        MP = Chromato.MP;
    end
    if isfield(Chromato,'SR')
        SR = Chromato.SR;
    end
    Chromato = Chromato.data;
end

[nb_measY,nb_measX] = size(Chromato); % size for 1st and 2nd D columns.

% figure

if and(exist('MP','var'),exist('SR','var'))
    if strcmpi(ShadingMode,'interp')
%         nb_measX=size(Chromato,2); % [pixel]
        if ~exist('X','var')
            Y = ((1:(MP*SR))'/SR)*ones(1,nb_measX);
            X = acquisition_delay+ones(MP*SR,1)*(1:nb_measX)*MP/60; % 60 = 60 [s]
        end
        pcolor(X,Y,Chromato); shading interp
    elseif strcmpi(ShadingMode,'flat')
        nb_measX=size(Chromato,2)+1; % [pixel]
        if ~exist('X','var')
            Y = ((1:(MP*SR+1))'/SR)*ones(1,nb_measX);
            X = acquisition_delay+ones(MP*SR+1,1)*(1:nb_measX)*MP/60; % 60 = 60 [s]
        else
            XBis = ones(size(X)+1);
            XBis(1:end-1,1:end-1) = X;
            XBis(:,end) = XBis(:,end-1) + (XBis(:,end-1)-XBis(:,end-2));
            XBis(end,:) = XBis(end-1,:);
            X = XBis; clear XBis
            YBis = ones(size(Y)+1);
            YBis(1:end-1,1:end-1) = Y;
            YBis(end,:) = YBis(end-1,:) + (YBis(end-1,:)-YBis(end-2,:));
            YBis(:,end) = YBis(:,end-1);
            Y = YBis; clear YBis 
        end
        Y = Y - ((Y(2,1)-Y(1,1))/2);
        X = X - ((X(1,2)-X(1,1))/2); % To have labels appear in middle 
%             of a given flat value box, not at the left of it.
        ChromatoBis = ones(size(Chromato)+1);
        ChromatoBis(1:end-1,1:end-1) = Chromato;
        Chromato = ChromatoBis; clear ChromatoBis
        pcolor(X,Y,Chromato); shading flat
    end
    xlabel('\bfFirst dimension retention time [min]','fontsize',FtSz)
    ylabel('\bfSecond dimension retention time [s]','fontsize',FtSz)
else
    if strcmpi(ShadingMode,'interp')
        if ~exist('X','var')
%             pcolor(Chromato); shading interp
            Y = ((1:nb_measY)')*ones(1,nb_measX); % pixels
            X = acquisition_delay + ones(nb_measY,1)*(1:nb_measX); % pixels
            pcolor(X,Y,Chromato); shading interp
        else
            pcolor(X,Y,Chromato); shading interp
        end
    elseif strcmpi(ShadingMode,'flat')
        ChromatoBis = ones(size(Chromato)+1);
        ChromatoBis(1:end-1,1:end-1) = Chromato;
        Chromato = ChromatoBis; clear ChromatoBis
        if ~exist('X','var')
            nb_measX = size(Chromato,2); % [pixel]
            Y = ((1:size(Chromato,1))')*ones(1,nb_measX);
            X = acquisition_delay+ones(size(Chromato,1),1)*(1:nb_measX);
        else
            XBis = ones(size(X)+1);
            XBis(1:end-1,1:end-1) = X;
            XBis(:,end) = XBis(:,end-1) + (XBis(:,end-1)-XBis(:,end-2));
            XBis(end,:) = XBis(end-1,:);
            X = XBis; clear XBis
            YBis = ones(size(Y)+1);
            YBis(1:end-1,1:end-1) = Y;
            YBis(end,:) = YBis(end-1,:) + (YBis(end-1,:)-YBis(end-2,:));
            YBis(:,end) = YBis(:,end-1);
            Y = YBis; clear YBis 
        end
        Y = Y - ((Y(2,1)-Y(1,1))/2);
        X = X - ((X(1,2)-X(1,1))/2); % To have labels appear in middle 
    %             of a given flat value box, not at the left of it.
        pcolor(X,Y,Chromato); shading flat
        
    end
    xlabel('\bfFirst dimension retention time [pixel]','fontsize',FtSz)
    ylabel('\bfSecond dimension retention time [pixel]','fontsize',FtSz)
end

if AutoAxis
    SortedData = sort(Chromato(:),'descend');
	try % In case there are e.g. many "Inf" data in the Chromatogram
	    % we avoid here to get an error, we would just go back to 
		% matlab default caxis limits.
		MaxAxis = SortedData(round(0.025*length(SortedData)));
		caxis([0,MaxAxis])
	end
end
    
if RedBlue 
    Coolormap=ones(64,3);
    Coolormap(1:41,2)=((1:41)')/41;
    Coolormap(1:41,3)=((1:41)')/41;
    Coolormap(44:64,1)=(flipud((1:21)'))/21;
    Coolormap(44:64,2)=(flipud((1:21)'))/21;
    if RedBlue==19 % little trick.
        Coolormap(44:64,1)=(flipud((11:0.5:21)'))/21;
        Coolormap(44:64,2)=(flipud((11:0.5:21)'))/21;
    end
    % colormap(flipud(Coolormap))
    Co_lor_map=flipud(Coolormap);
    colormap(Co_lor_map)
end 
    
if exist('titleFig','var')
    title(['\bf',titleFig],'fontsize',FtSz)
end

set(gca,'fontsize',FtSz,'XAxisLocation',XAxisLocation)

if exist('AlkNBlabel','var')
%     TickS = AlkNBlabel(:,1);
    TickSstr={};
    for zut=1:size(AlkNBlabel,1)% zut=1:length(TickS)
        TickSstr(zut)={num2str(AlkNBlabel(zut,1))};
    end
%     TickLocation = zeros(size(TickS));
%     for ty = 1:length(TickS)
%         TickLocation(ty) = AlkNBlabel(AlkNBlabel(:,1)==TickS(ty),2);
%     end
    set(gca,'Xtick',AlkNBlabel(:,2),'XTickLabel',TickSstr)
    xlabel('{\bf{\itn}-alkane number}','fontsize',FtSz)
end

if exist('color_map','var')
    colormap(color_map)
end


    
