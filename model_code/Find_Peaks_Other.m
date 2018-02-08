function [Peaks_Other_est] = Find_Peaks_Other(Other,Peaks_Ref,PlotVar,lm1,lm2,lm1_,lm2_)

% This function tries to determine the position of the list of alignment
% points "Peaks_Ref" (positions in a reference chromatogram), in another 
% chromatogram (called "Other"). It works as follow:
% First it looks for the local maximum value in a certain window around
% the position corresponding to a peak in the reference chromatogram in
% Other. Then, it verifies that this is really a local maximum, and if not
% it looks for it.
% BEWARE that this function should work pretty well if the
% displacements of the peaks are small and if the peaks chosen are far away
% from bigger peaks... In any case, it is highly recommended to verify the
% result obtained!
% 
% This function require the following intputs:
% - "Other" contains the chromatogram in which you want to find the 
% positions of the alignment points which positions are known in a
% reference chromatogram. This should be a two-dimensional matrix
% (for example with 2nd D column as rows and 1st D column as columns (the
% different parameters are optimised for this type of matrix))
% - "Peaks_Ref" contains the position of the alignment points in the
% reference  chromatogram. It should be a vector with 2 columns and as
% many lines as you want, each line giving the position of a point (please
% note that the position of these points should be maximum pixel of chosen
% peaks). If Other is of size [m,n], then the first column of Peaks_Ref
% refers to dimension n.
% - "PlotVar" should be 1 if you want to see a plot of Other with the
% estimated alignment points positions on it (for verification), if you
% don't want, give it a value of 0. Default value is 1.
% - "lm1", "lm2", "lm1_" and "lm2_" defines the window in which the code
% looks for local maximum (1 beeing for dimension corresponding to the
% first column of Peaks_Ref, 2 for the other one). The two first values
% are for the first search, and then the code verifies that the value
% obtained is really a local maximum in the window defined by the two last
% values, and if not continues searching until it finds one. Those values
% are optional, default values are given to it in the code if no values are
% specified. As it works iteratively until finding a local maximum, default 
% values should work in most of the cases. You can possibly try to give
% less big values if you fear that the code will look on a too large area
% and mistake one peak for another. These values depend on the instrument
% and especially sampling rate, and some adaptation of parameters values 
% may be required.
% 
% The output of the function is:
% - "Peaks_Other_est", the estimated positions of the alignment points in
% Other (might not be correct in some cases, a verification is
% recommended). The alignment points positions are in the same order as in
% Peaks_Ref.
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
% Switzerland, Environmental Chemistry Modeling Laboratory, 2014. 
% Authors : Jonas Gros and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values (if not defined by the user):
if nargin<7
    lm2_=10;
    if nargin<6
        lm1_=2;
        if nargin<5
            lm2=100;
            if nargin<4
                lm1=5;
                if nargin<3
                    PlotVar=1;
                end
            end
        end
    end
end

% Initialise "Peaks_Other_est":
Peaks_Other_est=zeros(size(Peaks_Ref,1),size(Peaks_Ref,2));

% Look for the local maximums (in Other) in the vicinity of the
% positions of the alignment points (in the reference chromatogram):
for k=1:size(Peaks_Ref,1)
% See localmax.m for more details about how the function works:
    [x y z]=localmax(Other,Peaks_Ref(k,1),Peaks_Ref(k,2),lm1,lm2);
    Peaks_Other_est(k,:)=[x y ];
end

% Verify that the value found above is really a local maximum, and if not
% look for the local maximum value:
for m=1:size(Peaks_Ref,1)
    pp=1;
    while pp
        [x,y,z]=localmax(Other,Peaks_Other_est(m,1),Peaks_Other_est(m,2),lm1_,lm2_);
        Peaks_Other_est2=[x,y];
        if Peaks_Other_est(m,:)==Peaks_Other_est2
            pp=0;
        else
            Peaks_Other_est(m,:)=Peaks_Other_est2;
        end
    end
end

% If asked, plot:
if PlotVar
    figure;
    plotChromato(Other,'shading','interp')
    colormap('jet')
    title('\bfEstimated position of the peaks...','fontsize',18)
    hold on
    plot(Peaks_Other_est(:,1),Peaks_Other_est(:,2),'om','linewidth',2)
    hold off
%     colorbar
    disp('Please verify the position of the points on the displayed graph!')
end

    