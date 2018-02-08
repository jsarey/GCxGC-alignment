
% function [x,y,z] = localmax(Z,x_guess,y_guess,x_span,y_span)
%
% where x_guess and y_guess are eyeball estimates of a peak location
% in the NxM array Z, and x_span and y_span designate the size of the
% subarray in the vicinity of x_guess and y_guess which will be searched
% for the true maximum value (peak). x and y are the location of the
% local peak, and z is its value.
% Version 2, 7 May 2013 (improved, search is now limited to possible space,
% avoids errors).
% (Before was searching a rectangle, leading to errors when the rectangle
% extended outside the searching space. Now the rectangle stops at the
% edges of the chromatogram).
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
% Switzerland, Environmental Chemistry Modeling Laboratory, 2014. 
% Authors : Jonas Gros and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [x,y,z] = localmax(Z,xg,yg,m,n)

W = Z(max((yg-n),1):min((yg+n),size(Z,1)),max((xg-m),1):min((xg+m),size(Z,2)));

z = max(max((W)));

[I,J]=ind2sub(size(W),find(W==z,1));
x=max((xg-m),1)+J-1;
y=max((yg-n),1)+I-1;



