function [Peak_Pos] = give_peak_pos(Chromato,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a function which allow to interactively select alignment points:
% the user draw a rectangle around it, and the local maximum within this
% rectangle is taken as the position of the alignment point.
% The input is:
% - "Chromato", the chromatogram in which you want to indicate alignment
% points positions.
% There is one option (use matlab usual syntax (,'option name', 'option
% value'))
% - the option 'EXPL', if it is 0 will avoid the display
% in matlab command window of the explanations about how to use the code
% (usefull when you have understood how it works...)
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
% Switzerland, Environmental Chemistry Modeling Laboratory, 2014. 
% Authors : Jonas Gros and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for i=1:2:length(varargin)
    if strcmpi(varargin{i},'EXPL')
        EXPL_opt=varargin{i+1};
    else
        error(['The option ''' varargin{i} ''' is unknown'])
    end
end

if ~exist('EXPL_opt','var')
    EXPL_opt = 0;
end

global Peak_Pos_c

figure; plotChromato(Chromato,'shading','flat'); colormap('jet')

if EXPL_opt==1
disp('    EXPLANATIONS:     ')
disp('On the displayed picture, zoom as much you want,')
disp('then press ''ENTER'' (or press it anyway if you don''t want')
disp('to zoom), then click twice to describe a rectangle')
disp('around your peak. This rectangle should not')
disp('contain any bigger pixel than the peak apex.')
disp('Then, if you want to add other points,')
disp('press ''w''')
disp('(if you want to zoom for the next point you will add,')
disp('if you don''t, just press ''ENTER'').')
disp('If you have introduced all the points you wanted and want to escape,')
disp('then just click ''q''.')
disp('If you are not satisfied with the last point entered,')
disp('press ''c'' to delete it and enter a new point.')
disp('Press ''e'' to delete the last point entered and escape.')
disp('When selecting a point, if the coloraxis doesn''t make you happy,')
disp('You can click a letter + ''enter'' to change it (and then select')
disp('your point normally. The options available are:')
disp('- ''d'' -> doubles the color axis (peaks get more visible)')
disp('- ''h'' -> halves the color axis (peaks get less visible)')
disp('- ''f'' -> increases by a half the color axis')
disp('- ''t'' -> triples the color axis')
disp('- ''z'' -> keeps 2/3 of the color axis')
disp('- ''o'' -> keeps one third of the color axis')
%pause(2)
end

but=1;
while but==1

zoom on; % use mouse button to zoom in or out
% Press Enter to get out of the zoom mode.

% CurrentCharacter contains the most recent key which was pressed after opening
% the figure, wait for the most recent key to become the return/enter key
try % Added because the new matlab version (2014) requires a character
    waitfor(gcf,'CurrentCharacter',char(13)) % as input, the oldest
catch % version required a numeric value (corresponding to the character...)
    waitfor(gcf,'CurrentCharacter',13)
end

zoom reset
zoom off

% The really important line is "[x,y,but] = ginput(2);". But the "while"
% etc. was added so that if the user by mistake presses 'enter' either
% before pressing with his mouse or after pressing it only once instead of
% twice, the code will not stop, diplaying a mistake (and losing all the
% points positions found so far). Instead, it just begins again with the
% introduction of the point.
x=[];
y=[];
while length(x)<2
[x,y,but] = ginput(2);
if size(but)==1
   if but==100 % 'd' -> double
      caxis(caxis/2)
   elseif but==104 % 'h' -> half
       caxis(caxis*2)
   elseif but==102 % 'f' -> fifty percent
       caxis(caxis/1.5)  % = caxis(caxis*2/3)
   elseif but==116 % 't' -> triple
       caxis(caxis/3)
   elseif but==122 % 'z' -> two-third
       caxis(caxis*3/2)
   elseif but==111 % 'o' -> one third
       caxis(caxis*3)
   end
end
end


% Plot the rectangle indicated by the user:
hold on
plot([x(1),x(1)],[y(1),y(2)],'k')
plot([x(2),x(2)],[y(1),y(2)],'k')
plot([x(1),x(2)],[y(1),y(1)],'k')
plot([x(1),x(2)],[y(2),y(2)],'k')

% Verify that the rectangle contains some pixels in each dimension, and if
% yes compute the maximal value of the peak:
% alignment point):
if and(ceil(min(y(1),y(2)))<=floor(max(y(1),y(2))),(ceil(min(x(1),x(2)))<=floor(max(x(1),x(2)))))
Peak_apex_value=max(max(Chromato(ceil(min(y(1),y(2))):floor(max(y(1),y(2))),ceil(min(x(1),x(2))):floor(max(x(1),x(2))))));
else
drawnow; % <- flush the graphic events
commandwindow; % <- Bring to focus
error('Please, your rectangle should contain at least one chromatogram pixel in each dimension!')

end
% Find the position of the maximal pixel of the peak (= the position of the
% alignment point):
for k=ceil(min(y(1),y(2))):floor(max(y(1),y(2)))
for l=ceil(min(x(1),x(2))):floor(max(x(1),x(2)))
if Chromato(k,l)==Peak_apex_value
Peak_apex=[l,k]; % l=1st D, k= 2nd D
end
end
end
% Plot the alignment point position just determined:
plot(Peak_apex(1),Peak_apex(2),'ok')
% Store the result:
if exist('Peak_Pos','var')
Peak_Pos=[Peak_Pos;Peak_apex];
else
Peak_Pos=Peak_apex;
end

% Now the user can chose what he wants to do:
hj = waitforbuttonpress;
% If the user wants to stop now:
key = get(gcf,'CurrentCharacter') ;
if abs(key)==113 % -> click 'q' to quit
but=19;
drawnow; % <- flush the graphic events
commandwindow; % <- Bring to focus
end
% If the user is not happy with the last alignment point position
% determined:
if or(abs(key)==67,abs(key)==99) % key = 'c' or 'C'
Peak_Pos=[Peak_Pos(1:end-1,:)];
plot([x(1),x(1)],[y(1),y(2)],'r')
plot([x(2),x(2)],[y(1),y(2)],'r')
plot([x(1),x(2)],[y(1),y(1)],'r')
plot([x(1),x(2)],[y(2),y(2)],'r')
plot(Peak_apex(1),Peak_apex(2),'or')
end
% If the user is not happy with the last alignment point position and want
% to escape:
if abs(key)==101 % key = 'e'
Peak_Pos=[Peak_Pos(1:end-1,:)];
plot([x(1),x(1)],[y(1),y(2)],'r')
plot([x(2),x(2)],[y(1),y(2)],'r')
plot([x(1),x(2)],[y(1),y(1)],'r')
plot([x(1),x(2)],[y(2),y(2)],'r')
plot(Peak_apex(1),Peak_apex(2),'or')
but=19;
drawnow; % <- flush the graphic events
commandwindow; % <- Bring to focus
end
% Unless the user wants to keep this zooming without doing more additional
% zooming for next step (in which case he presses 'enter'), zoom out:
if key~=13 % 13 = ENTER, 97 = a
axis([0,size(Chromato,2),0,size(Chromato,1)])
end
Peak_Pos_c=Peak_Pos;
end


end