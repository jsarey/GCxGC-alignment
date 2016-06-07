function [Aligned,Displacement,Deform_output] = alignChromato(Ref,Target,Peaks_Ref,Peaks_Target,varargin)

% This file is a two-dimensional chromatogram alignment function. 
% The user has to identify points (called "alignment points")
% both in a reference chromatogram (called "Ref") and in a target
% chromatogram (called "Target" in this function) which will be aligned to 
% the reference chromatogram (the resulting chromatogram will be called 
% "Aligned" in this function).
% 
% Refer to [Gros, J.; Nabi, D.; Dimitriou-Christidis, P.; Rutler,
% R.; Arey, J. S., "Robust algorithm for aligning two-dimensional 
% chromatograms", Analytical Chemistry, 9033-9040, 2012.] for more 
% explanations, examples of results and discussion of the interest 
% of the method.
% 
% This works as follows:
% 
% The alignment points (defined here as maximum pixels of selected peaks)
% are perfectly aligned.
% 
% DISPLACEMENT INTERPOLATION:
% 
% For each other point, the displacement is interpolated.
% Interpolation in the first dimension is done using a simple linear
% interpolation (according only to first dimension postition) between the
% alignment points, and a linear extrapolation based on the information
% given by the first and last alignment points. (This avoids spurious
% displacement estimates on the edges, especially if the two last (or
% first) alignment points in the first dimension are very close, which can
% give exagerated displacement estimates on the edges (especially when a
% one pixel displacement difference results from a sub-pixel variation,
% rounded (the peak position may vary by one pixel for a very small change
% of elution time...)))).
% In second dimension:
% Default interpolation is Sibson natural-neighbor.
% Sibson natural-neighbor interpolation uses Voronoï diagrams (Thiessen
% polygons in other words): a Voronoï diagram is drawn around the alignment
% points, then for each pixel, if a diagram would also be drawn around it,
% it would take some area to neighboring polygons. The weights for each
% alignment point is proportional to the proportion of the area of the new
% polygon coming from the previous area of the polygon corresponding to
% this alignment point. 
% Please note that this interpolation technique works only within the
% convex hull of the alignment points. Thus to extrapolate to the rest of
% the chromatogram, 4 peaks are added (one at each corner (in fact on the
% corners of a rectangle 10% bigger than the chromatogram)), which
% displacement is extrapolated using a weighted mean to the square of the
% inverse of the distance.
% 
% INTENSITY VALUE INTERPOLATION:
% 
% Then for each point in the Aligned chromatogram,
% the code looks where it should come from in the unaligned target
% chromatogram ("Target"), and interpolate the intensity value (if it
% falls between 4 (or 2) different pixels) with a bicubic interpolation
% (by default).
% 
% Two-dimensional chromatograms properties are taken into
% account:
% - The fact that the two-dimensional shape of the chromatogram is only a
% reconstruction a posteriori while the measure is realised continuously at
% the exit of the second column, thus a deformation of a chromatogram can
% lead some compounds from the above part of the chromatogram to the bottom
% of the next column or vice versa. (Beware: the actual version only allows
% a transfer of a pixel no more than a half column away (to limit memory 
% use). However, a greater
% extent of allowed wrap-around is easy to implement (some lines of code to
% modify))
% - The volume below a peak is linked to the quantity of compound present
% in the sample analysed, thus this volume should ideally not be
% modified through alignment. Thus, a correction of the intensity of the
% signal at each pixel in function of deformation (stretching or
% compression) was applied afterwards to the whole chromatogram, as the
% alignment procedure described above would just enlarge or compress a
% peak, without any correction for the height (intensity) of the peak (a
% thinner peak is higher while a broader peak is smaller if the quantity
% of compound remains unchanged).
% It was decided
% to estimate deformation at each pixel of the chromatogram (except on the
% edges) as the difference between the positions (on the unaligned
% chromatogram) from where should have come the pixel before and the pixel
% after the pixel considered, divided by this distance on the aligned
% chromatogram (id est 2 pixels). The correction was computed in both
% dimensions and then applied, this for the whole chromatogram.
% 
% This function requires the following ARGUMENTS:
% - "Ref" contains the reference chromatogram on which the target
% chromatogram will be aligned. This should be a two-dimensional matrix
% with 2nd D column as rows and 1st D column as columns. (each row is for a
% different 2nd D value, each column for a different 1st D value)
% Please note that
% it is important that this order is respected as the code allows
% wrap-around between the top of a column and the bottom of the next one,
% whenever a wrap-around between the beginning and end of a line would
% have no meaning.
% - "Target" contains the target chromatogram you want to align on the
% reference chromatogram. This should also be a two-dimensional matrix, of
% the same size than Ref.
% - "Peaks_Ref" contains the position of the alignment points in Ref.
% it should be a vector with 2 columns and as many lines as you want, each
% line giving the position of a point (please note that these points, as
% you have to find them both on Ref and on Target, should be maximum pixels
% of chosen peaks). Please note also that if Ref and Target are of size
% [m,n], Peaks_Ref(:,1) are values between 1 and n). 
% The positions should be expressed as
% "pixels" (a pixel is the distance between two measurements, the first
% measurement is at the position one pixel in both dimensions ([1,1]))
% - "Peaks_Target" contains the position of the alignment points in Target
% (i.e. in the target chromatogram).
% It should be a vector of the same size then Peaks_Ref, with alignment
% point positions indicated in the same order. You may use the function
% Find_Peaks_Other as an help to determine Peaks_Target from Peaks_Ref.
% - Please note that, if you want to use the version of the code described
% in the article, you will have to give to the code the widths of a
% "typical peak" in both dimensions (but this is an optional input, if you
% don't give this information, distances will just be given in units of
% [pixel]) (the use of the article version is recommended, whenever this
% parameter has only a marginal effect on the alignment quality).
% 
% This function accepts the following OPTIONAL ARGUMENTS
% (it was decided to use Matlab usual practice, thus if you want to give to
% the function a value for an optional argument, you have to proceed like
% this: Align2DChrom(Ref,Target,Peaks_Ref,Peaks_Target,'optional argument
% name','optional argument value'); (where 'optional argument name' is a
% string and 'optional argument value' is any kind of Matlab variable,
% depending upon the option):
% - "Peak_widths" is a vector (of the type
% [1st D typical peak width, 2nd D typical peak width] 
% containing two values, which gives an estimate (as a number of pixels) of
% the "typical peak width", the first number being for the first dimension,
% the second one for the second dimension. This can be roughly
% determined based on visual inspection of one or some "typical" peaks.
% Please note that if this option is not given to the code, distances will
% be considered in units of [pixel] and not in units of
% [typical peak width], what will change the result of the
% interpolation/extrapolation of the second dimension displacement.
% - "DisplacementInterpMeth" is a string which value is 
% 'natural-neighbor' to indicate which interpolation/extrapolation method
% has to be used for the displacement interpolation in second dimension 
% ('natural-neighbour' (Great Britain English) is also
% accepted). 'natural-neighbor' is used by default if nothing is indicated
% by the user about this point.
% (other option value is outdated).
% - "PowerFactor" is a number indicating the power of the inverse of the
% distance for the weighted average. A value of 2 seemed to work better
% than a value of 1... (in fact 1.8 was maybe the best in our case
% But this should depend upon the number of alignment
% points, of the particular deformations found in the chromatograms to
% align, etc.). If no value is
% given by the user for PowerFactor, a default value of 2 will be used.
% - "SaveVar" is another optional variable: it should be a string
% corresponding to the name of the file in which Aligned should be saved
% (for example if its value is 'example', Aligned would be saved in the
% file example.mat, as a 2D matrix).
% - "InterPixelInterpMeth" specifies the inter-pixel interpolation method
% (intensity value interpolation method)
% (when the displacement is applied, the value at a certain pixel in the
% aligned chromatogram should usually come from between four pixels of the
% target chromatogram; the interpolation specified in this variable is used
% for this). The values could be:
% > "spline", the Matlab implementation of the cubic spline interpolation
% (using the function interp2). // it is the least recommended for peak
% volumes conservation.
% > "linear", the matlab implementation of the bilinear interpolation
% (using the function interp2). // Gave less good
% results for peak volumes/areas computed by GC Image than the default
% value.
% > "cubic" (default value), the matlab implementation of bicubic
% interpolation (using the function interp2). // Is the best for minimal
% peak volumes alterations as computed with GC Image. The
% "linear" option (bilinear interpolation) might be tried if peaks volumes
% are to be computed in a different way, and comparisons done.
% 
% The OUTPUTS of the function are:
% - "Aligned", a 2D chromatogram of the same size than Target which is a
% transformed version of Target which should be better aligned to Ref. If
% this would not be enough the case, you are adviced to verify that your
% alignment points really correspond between Peaks_Ref and Peaks_Target, to
% try other points or more points (as a guideline, try to have an alignment
% point in each region of the chromatogram showing a different shifting
% trend, and try also to have most of the interesting part of your
% chromatogram within the convex hull of the alignment points).
% Please note that you can use the syntax:
% Aligned=Align2DChrom(Ref,Target,Peaks_Ref,Peaks_Target,varargin);
% if you are only interested by this output.
% - "Displacement" is a three dimensional matrix of size [m,n,2].
% [:,:,1] contains the estimated displacerment in the 2nd D for each pixel
% of the reference chromatogram (which is of size [m,n]), while [:,:,2]
% contains the estimated displacement value in the 1st D for each pixel of
% the reference chromatogram.
% - "Deform_output" contains the deformation correction factor which was
% applied to each pixel (this is a matrix of the same size than
% Displacement, containng the correction due to stretching/compression
% independently for both dimensions. The total correction that is applied
% corresponds to:
% (Deform_output(:,:,1).*Deform_output)/4
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
% Switzerland, Environmental Chemistry Modeling Laboratory, 2014. 
% Authors : Jonas Gros and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%  - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 -  %

% First, if ever the user used some options, give their values to the
% corresponding variables:
for i=1:2:length(varargin)
    if strcmpi(varargin{i},'DisplacementInterpMeth')
        DisplacementInterpMeth=varargin{i+1};
    elseif strcmpi(varargin{i},'PowerFactor')
        PowerFactor=varargin{i+1};
    elseif strcmpi(varargin{i},'SaveVar')
        SaveVar=varargin{i+1};
    elseif strcmpi(varargin{i},'Peak_widths')
        Peak_widths=varargin{i+1};
        elseif strcmpi(varargin{i},'InterPixelInterpMeth')
        InterPixelInterpMeth=varargin{i+1};
    else
        error(['The option ''' varargin{i} ''' is unknown'])
    end
end

% If DisplacementInterpMeth was not given, set it to 'natural-neighbor':
if ~exist('DisplacementInterpMeth','var')
     DisplacementInterpMeth='natural-neighbor';
end

% If PowerFactor was not given, set its value to a default value of 2:
if ~exist('PowerFactor','var')
    PowerFactor=2;
end

% If Peak_widths was not given, it will be considered that no peak width
% correction has to be added (i.e. that distances will be considered in
% units of pixels), this is achieved by giving the same value to the two
% elements of Peak_widths:
if ~exist('Peak_widths','var')
     Peak_widths=[1,1];
end

% If InterPixelInterpMeth was not given, set it to cubic:
if ~exist('InterPixelInterpMeth','var')
InterPixelInterpMeth='cubic';
end

% The code issues an error if Ref and Target do not have the same size.
% In case they do not have the same size, let's here complete the smallest
% one with zeros to get the same size:
[Ref,Target] = equalize_size(Ref,Target);
	

% normalized distance:
PeakWidth2ndD=Peak_widths(2);
PeakWidth1stD=Peak_widths(1);

% This is the displacement (between the Reference and Target) of the
% alignment points:
Peaks_displacement=-(Peaks_Ref-Peaks_Target);
% This will contain the aligned chromatogram:
Aligned=zeros(size(Target,1),size(Target,2));
% This will contain the displacement (between Ref and Target) for
% each point of the chromatogram (2 values for each point, one for each 
% dimension):
Displacement=zeros(size(Target,1),size(Target,2),2);

if or(strcmpi(DisplacementInterpMeth,'natural-neighbor'),...
        strcmpi(DisplacementInterpMeth,'natural-neighbour'))
    %     "Displacement2" will contain the displacement for the four added
    %     points (placed on the corners of a rectangle ~10% bigger than the
    %     chromatogram, placed around it):
    Displacement2=zeros(2,2,2);

    for w=-floor(0.05*size(Aligned,1)+1):...
            size(Aligned,1)+1+floor(0.05*size(Aligned,1)+1)+...
            ceil(0.05*size(Aligned,1)+1):...
            size(Aligned,1)+1+ceil(0.05*size(Aligned,1)+1)
        for x=-floor(0.05*size(Aligned,2)+1):...
                size(Aligned,2)+1+floor(0.05*size(Aligned,2)+1)+...
                ceil(0.05*size(Aligned,2)+1):...
                size(Aligned,2)+1+ceil(0.05*size(Aligned,2)+1)
    %             Distance of the pixel to the alignment points, as a vector:
            Distance_vec=ones(size(Peaks_Ref,1),1)*[w x]-fliplr(Peaks_Ref);   

            Distance=sqrt((Distance_vec(:,1).^2)+((Distance_vec(:,2)*PeakWidth2ndD/PeakWidth1stD).^2));

    % Displacement of the pixel (it's a weighted mean of the displacement of
    % each alignment point to a power (defined by PowerFactor) of the
    % inverse of the distance to the pixel considered):
            Displacement2(1+(w+floor(0.05*size(Aligned,1)+1))/...
            (size(Aligned,1)+1+floor(0.05*size(Aligned,1)+1)+...
            ceil(0.05*size(Aligned,1)+1)),1+(x+floor(0.05*size(Aligned,2)+1))/...
            (size(Aligned,2)+1+floor(0.05*size(Aligned,2)+1)+...
            ceil(0.05*size(Aligned,2)+1)),:)=...
            sum(fliplr(Peaks_displacement).*(1./(Distance.^PowerFactor)*...
            ones(1,2)))/sum(1./(Distance.^PowerFactor));

        end
    end
    % Add a peak at each corner (on the corners of a rectangle 10% bigger than
    % the chromatogram placed around it):
    Peaks_Ref=[Peaks_Ref;[0 0]+[-floor(0.05*size(Aligned,2)+1) -floor(0.05*size(Aligned,1)+1)];...
        [0 size(Aligned,1)+1]+[-floor(0.05*size(Aligned,2)+1) ceil(0.05*size(Aligned,1)+1)];...
        [size(Aligned,2)+1 0]+[ceil(0.05*size(Aligned,2)+1) -floor(0.05*size(Aligned,1)+1)];...
        [size(Aligned,2)+1 size(Aligned,1)+1]+[ceil(0.05*size(Aligned,2)+1) ceil(0.05*size(Aligned,1)+1)]];
    % Add the corresponding extrapolated displacement values:
    Peaks_displacement=[Peaks_displacement;...
        [Displacement2(1,1,2) Displacement2(1,1,1)];...
        [Displacement2(2,1,2) Displacement2(2,1,1)];...
        [Displacement2(1,2,2) Displacement2(1,2,1)];...
        [Displacement2(2,2,2) Displacement2(2,2,1)]];

    % Do the natural-neighbor interpolation for the normalized distances:
     Fdist1=TriScatteredInterp(Peaks_Ref(:,2),Peaks_Ref(:,1)*PeakWidth2ndD/PeakWidth1stD,...
        Peaks_displacement(:,2),'natural');

end

% Linear interp within the convex hull of the alignment points:
Hep=[Peaks_Ref(1:end-4,1),Peaks_displacement(1:end-4,1)]; % The four last 
% values are the added alignment points at the corner of the rectangle 10%
% bigger than the chromatogram.
% The things with the variables below with strange names is there to avoid
% bugs when two alignment points would have the same position in 1st D.
% Then a mean displacement is used.
Hep1=Hep(:,1);
Hep2=Hep(:,2);
Herp=[];
for hHh=1:size(Hep,1)
    Herp(hHh,1)=Hep1(hHh);
    Herp(hHh,2)= sum(Hep2(Hep1==Hep1(hHh)))/sum(Hep1==Hep1(hHh));
end
Hap=[];
k=1;
puet=1;
for kui=1:size(Herp,1)
    for zui=1:size(Hap,1)
        if Herp(kui,1)==Hap(zui,1)
            puet=0;
            Peaks_displacement(Herp(:,1)==Hap(zui,1),1)=Herp(zui,2);
        end
    end
    if puet
        Hap(k,:)=Herp(kui,:);
        k=k+1;
    else
        puet=1;
    end
end
Hum=interp1(Hap(:,1),Hap(:,2),1:(size(Target,2)),'linear','extrap');
Pks=Peaks_Ref(1:end-4,1);
Displ1D=Peaks_displacement(1:end-4,1);
% "Special" linear extrap, will be used outside of the convex hull of the
% alignment points:
Hum2=interp1([mean(Pks(Pks==min(Peaks_Ref(1:end-4,1)))),mean(Pks(Pks==max(Peaks_Ref(1:end-4,1))))],...
[mean(Displ1D(Pks==min(Peaks_Ref(1:end-4,1)))),mean(Displ1D(Pks==max(Peaks_Ref(1:end-4,1))))],...
1:(size(Target,2)),'linear','extrap');

for w=1:size(Aligned,1)
    for x=1:size(Aligned,2)
%         The "if" is here to avoid re-defining the value of the peaks used
%         for alignment.
        if Aligned(w,x)==0

%     Interpolation result is a kind of function which has to be evaluated
%     in [w,x]:
% And in 1st D the linear interpolation/extrapolation is used:
                if and(x>=min(Peaks_Ref(1:end-4,1)),x<=max(Peaks_Ref(1:end-4,1)))
                    Displacement(w,x,:)=[Fdist1(w,x*PeakWidth2ndD/PeakWidth1stD) ; Hum(x)];
                else
                    Displacement(w,x,:)=[Fdist1(w,x*PeakWidth2ndD/PeakWidth1stD) ; Hum2(x)]; 
                end
        end
    end
end
X=ones(size(Ref,1)*2,1)*(1:size(Ref,2));
Y=(((1-round(size(Ref,1)/2)):(size(Ref,1)+...
    (size(Ref,1)-round(size(Ref,1)/2))))')*ones(1,size(Ref,2));
Z=zeros(size(Target,1)*2,size(Target,2));
Z(1:floor(size(Target,1)/2),2:end)=Target(floor(size(Target,1)/2+1:size(Target,1)),1:end-1);
    Z(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)=Target;
Z((floor(size(Target,1)/2)+size(Target,1)+1):end,1:end-1)=...
        Target(1:ceil(size(Target,1)/2),2:end);
% % % if strcmpi(InterPixelInterpMeth,'spline')
% % %     Aligned=interp2(X,Y,Z,...
% % %         X(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
% % %         Displacement(:,:,2),...
% % %         Y(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
% % %         Displacement(:,:,1),'spline');
% % % elseif strcmpi(InterPixelInterpMeth,'cubic')
% % %     Aligned=interp2(X,Y,Z,...
% % %         X(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
% % %         Displacement(:,:,2),...
% % %         Y(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
% % %         Displacement(:,:,1),'cubic');
% % % elseif strcmpi(InterPixelInterpMeth,'linear')
% % %     Aligned=interp2(X,Y,Z,...
% % %         X(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
% % %         Displacement(:,:,2),...
% % %         Y(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
% % %         Displacement(:,:,1),'linear');
% % % end
% The function interp2 in matlab suffers an error if using the 2012a matlab
% version ( https://www.mathworks.com/support/bugreports/886257 ).
% As the user might be using the version of interp2 that suffers the error,
% we will use the work-around proposed in the bug report which link is
% given above:
if strcmpi(InterPixelInterpMeth,'spline')
    Xq = X(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
        Displacement(:,:,2);
    Yq = Y(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
        Displacement(:,:,1);
    Aligned=interp2(X,Y,Z,...
        Xq(:),...
        Yq(:),'spline');
    Aligned = reshape(Aligned,size(Xq));
elseif strcmpi(InterPixelInterpMeth,'cubic')
    Xq = X(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
        Displacement(:,:,2);
    Yq = Y(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
        Displacement(:,:,1);
    Aligned=interp2(X,Y,Z,...
        Xq(:),...
        Yq(:),'cubic');
    Aligned = reshape(Aligned,size(Xq));
elseif strcmpi(InterPixelInterpMeth,'linear')
    Xq = X(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
        Displacement(:,:,2);
    Yq = Y(round(size(Target,1)/2)+1:(round(size(Target,1)/2)+size(Target,1)),:)+...
        Displacement(:,:,1);
    Aligned=interp2(X,Y,Z,...
        Xq(:),...
        Yq(:),'linear');
    Aligned = reshape(Aligned,size(Xq));
end
% If values should come from outside of the chromatogram ("before it or
% after it"), set the value to 0.
Aligned(isnan(Aligned))=0;

% Add the alignment peaks diplacements into the variable Displacement
% (except the four points added to allow extrapolation with Sibson's
% natural-neighbor interpolation technique):
    for k=1:size(Peaks_displacement,1)-4
        Displacement(Peaks_Ref(k,2),Peaks_Ref(k,1),:)=fliplr(Peaks_displacement(k,:));
    end

% These two variables will contain the deformation in each dimension:
Deform1=zeros(size(Aligned,1),size(Aligned,2));
Deform2=zeros(size(Aligned,1),size(Aligned,2));

% Create a variable containing the interpolated/extrapolated displacements
% not only for the chromatogram, but also for the pixels just adjacent to
% it (but first, just initialise it):
Displacement_Extended=zeros(size(Aligned,1)+2,size(Aligned,2)+2,2);
% Introduce the values already computed:
Displacement_Extended(2:end-1,2:end-1,:)=Displacement(:,:,:);
% Do the interpolation for the normalized distances again, with a 1 pixel
% move forward
 Fdist1bis=TriScatteredInterp(Peaks_Ref(:,2)+1,(Peaks_Ref(:,1)+1)*PeakWidth2ndD/PeakWidth1stD,...
    Peaks_displacement(:,2),'natural');
Humbis=interp1(Hap(:,1),Hap(:,2),0:(size(Target,2)+1),'linear','extrap');
Hum2bis=interp1([mean(Pks(Pks==min(Peaks_Ref(1:end-4,1)))),mean(Pks(Pks==max(Peaks_Ref(1:end-4,1))))],...
[mean(Displ1D(Pks==min(Peaks_Ref(1:end-4,1)))),mean(Displ1D(Pks==max(Peaks_Ref(1:end-4,1))))],...
0:(size(Target,2)+1),'linear','extrap');
% "Special" linear extrap, will be used outside of the convex hull of the
% alignment points:
Hum2=interp1([mean(Pks(Pks==min(Peaks_Ref(1:end-4,1)))),mean(Pks(Pks==max(Peaks_Ref(1:end-4,1))))],...
[mean(Displ1D(Pks==min(Peaks_Ref(1:end-4,1)))),mean(Displ1D(Pks==max(Peaks_Ref(1:end-4,1))))],...
0:(size(Target,2)+1),'linear','extrap');
% Compute and add the values on the edges:
for w=[1 size(Aligned,1)+2]
    for x=1:(size(Aligned,2)+2)

            % %     Interpolation result is a kind of function which as to be evaluated
            % %     in [w,x], and the linear interpolation/extrapolation is used for
            % the 1st D:
            if and(x>=min(Peaks_Ref(1:end-4,1)),x<=max(Peaks_Ref(1:end-4,1)))
                Displacement_Extended(w,x,:)=[Fdist1bis(w,(x)*PeakWidth2ndD/PeakWidth1stD) ; Hum(x)];
            else
                Displacement_Extended(w,x,:)=[Fdist1bis(w,(x)*PeakWidth2ndD/PeakWidth1stD) ; Hum2(x)]; 
            end

    end
end
for w=2:(size(Aligned,1)+1)
    for x=[1 size(Aligned,2)+2]

            %     Interpolation result is a kind of function which as to be evaluated
            %     in [w,x], and linear interpolation/extrapolation for 1st D
            % displacement:
            if and(x>=min(Peaks_Ref(1:end-4,1)),x<=max(Peaks_Ref(1:end-4,1)))
                Displacement_Extended(w,x,:)=[Fdist1bis(w,(x)*PeakWidth2ndD/PeakWidth1stD) ; Hum(x)];
            else
                Displacement_Extended(w,x,:)=[Fdist1bis(w,(x)*PeakWidth2ndD/PeakWidth1stD) ; Hum2(x)]; 
            end

    end
end

% The correction is computed for each pixel:
for w=1:size(Aligned,1)
    for x=1:size(Aligned,2)
        Deform1(w,x)=2+(-Displacement_Extended((w-1)+1,x+1,1)+...
        Displacement_Extended((w+1)+1,x+1,1));
        Deform2(w,x)=2+(-Displacement_Extended(w+1,(x-1)+1,2)+...
        Displacement_Extended(w+1,(x+1)+1,2));
%     if the value is bigger than 2, there is compression, otherwise
%     enlargement.
    end
end

if or(any(any(Deform1<0)),any(any(Deform2<0)))
    warning(['Some displacement predictions lead to inversion of pixel order ',...
        '(as indicated by negative values in Deform_output, here set to zero) ',...
        'please change your set of alignment points!'])
    figure; plotChromato(double(or(Deform1<0,Deform2<0))); caxis([0,1]);
    title('\bfProblematic pixels in red (set to 0 in Aligned chromatogram)')
    Deform1(Deform1<0) = 0;
    Deform2(Deform2<0) = 0;
end

% The correction is applied:
Aligned=Aligned.*(Deform1.*Deform2)/4;

% If the user wants to save the results:
if exist('SaveVar','var')
    save([SaveVar '.mat'],'Aligned','-double');
end

% Give the deformation correction for each pixel as an ouput of the
% function:
SzE=[size(Deform1),2];
Deform_output=zeros(SzE);
Deform_output(:,:,1)=Deform1;
Deform_output(:,:,2)=Deform2;

%%%%%%%%%

% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
% Switzerland, Environmental Chemistry Modeling Laboratory, 2014. 
% Authors : Jonas Gros and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
