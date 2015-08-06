
% Takes as input two GCxGC chromatograms (Chromato1 and Chromato2),
% and adds zeros where needed to make them of equal size.
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
% Switzerland, Environmental Chemistry Modeling Laboratory, 2015. 
% Authors : Jonas Gros.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [Chromato1_eq,Chromato2_eq] = equalize_size(Chromato1,Chromato2);

% Sizes of input chromatograms:
size1 = size(Chromato1);
size2 = size(Chromato2);

if or(size1(1)~=size2(1),size1(2)~=size2(2))
	% Sizes of output chromatograms:
	max_size = [max([size1(1),size2(1)]),max([size1(2),size2(2)])];

	% Initialise (adds the zeros):
	Chromato1_eq = zeros(max_size);
	Chromato2_eq = zeros(max_size);

	% Add the chromatographic data into it (the zeros will be on the outside):
	Chromato1_eq(1:size1(1),1:size1(2)) = Chromato1;
	Chromato2_eq(1:size2(1),1:size2(2)) = Chromato2;
else
	Chromato1_eq = Chromato1;
	Chromato2_eq = Chromato2;
end
