(* ::Package:: *)

(*

A package for orientation data analysis in structural geology

William C. Haneberg
bill@haneberg.com
www.haneberg.com

version 0.9

June 16, 2011

copyright \[CapitalAHat]\[Copyright]2011 William C. Haneberg

*)

BeginPackage["Stereonet`"]

       
OuterCircle::usage="OuterCircle[ ] is a Stereonet graphics primitive that draws the outline and ticks, and labels them. You likely won't have much use for it."

CartesianToCompass::usage="CartesianToCompass[\[Theta]] converts angles measured counterclockwise from the +x axis (Cartesian or mathematical convention) to an azimuth measured clockwise from North (compass convention)"

DiplineToPole::usage="DiplineToPole[{\[Delta], \[Theta]}], where \[Delta] is the dip angle and \[Theta] is the dip direction, converts the dip vector of a plane to the lower hemisphere pole to that plane."

PoleToDipline::usage="PoleToDipline[{\[Delta], \[Theta]}] converts a lower hemisphere pole to a plane with plunge \[Delta] and azimuth \[Theta] to the dip vector for that plane."

DiplineToCosines::usage="DiplineToCosines[{\[Delta], \[Theta]}], where \[Delta] is the dip angle and \[Theta] is the dip direction, converts a dip vector to its direction cosines. It can be used with any vector orientation measurement, for example the plunge and azimuth of a lineation."

CosinesToDipline::usage="CosinesToDipline[{x, y, z}] converts a set of direction cosines to the the dip angle and dip vector (or plunge and azimuth)."

DiplineToStrike::usage="DiplineToStrike[{\[Delta], \[Theta]}], where \[Delta] is the dip angle and \[Theta] is the dip direction, converts the dip vector to strike and dip using the right-hand rule."

StrikeToDipline::usage="StrikeToDipline[{\[Theta], \[Delta]}], where \[Delta] is the dip angle and \[Theta] is the strike, to a dip vector."

DiplineToNormal::usage="DiplineToNormal[{\[Delta], \[Theta]}], where \[Delta] is the dip angle and \[Theta] is the dip direction, converts the dip vector to the upward-directed normal to the plane."

ApparentDip::usage="ApparentDip[{\[Delta], \[Alpha]}], where \[Delta] is the true dip angle of the plane, calculates the apparent dip along a line \[Alpha] degrees from the true dip direction."

MeanVector::usage="MeanVector[data] calculates the mean vector from a list of dip vectors or lineations."

AngularVariance::usage="AngularVariance[data] calculates the angular variance of a list of dip vectors or lineations."

AngularStdDev::usage="AngularStdDev[data] calculates the angular standard deviation of a list of dip vectors or lineations."

FisherK::usage="FisherK[data] calculates Fisher's K for a list of dip vectors or lineations."

AngularConfidenceRadius::usage="AngularConfidenceRadius[data, \[Alpha]] calculates the angular confidence radius about the mean for a list of dip vectors or lineations for a specified level of significance. If no level of significance is given, the default value of \[Alpha] = 0.05 is used."

CroninBeta::usage="CroninBeta[data] calculates Cronin's \[Beta] statistic quantifying uncertainty in the strike angle at the \[Alpha] = 0.05 level of significance."

EigenFabricAnalysis::usage="EigenFabricAnalysis[data], where data is a list of dip vectors or lineations, returns 1) the mean vector, 2) the three eigenvalues, and 3) the three eigenvectors for the data."

HemisphericalVectorAngle::usage="HemisphericalVectorAngle[u, v] calculates the smallest angle between vectors u and v as measured within the lower hemisphere."

OutlierRemoval::usage="OutlierRemoval[data, \[Psi]] uses the Mahtab and Yegulalp (1982) method as described in Priest (1985), where \[Psi] is an angular tolerance measured in degrees. Using large values of \[Psi] relaxes the filtering and removes fewer outliers, whereas using small values of \[Psi] tightens the filtering and removes more outliers."

EqualAreaPoint::usage="EqualAreaPoint[{\[Delta], \[Theta]}, r] plots a point of radius r (relative to the total plot size) at dip angle \[Delta] and dip direction \[Theta] on a lower hemisphere equal area projection."

ConeSphereIntersection::usage="ConeSphereIntersection[\[Alpha]] generates a set of {x, y, z} points defining the intersection of a unit sphere with a vertical cone of angle \[Alpha]. The points can then be rotated into any position on the unit sphere."

CirclePts::usage="CirclePts[{\[Delta], \[Theta]}, \[Alpha]] calculates a list of points defining a circle of angular radius \[Alpha] and centered at {\[Delta], \[Theta]} on a unit lower hemisphere. It doesn't plot anything; rather, it simply calculates a list of points to be used by other functions. Try EqualAreaCirclePlot if you want to plot the equal area projection of a circle on the lower hemisphere."

EqualAreaCirclePlot::usage="EqualAreaCirclePlot[{\[Delta], \[Theta]}, \[Alpha]] draws the equal area projection of a circle of angular radius \[Alpha] and centered at {\[Delta], \[Theta]} on a unit lower hemisphere. It can accept most graphics style options."

ListEqualAreaPointPlot::usage="ListEqualAreaPointPlot[data, pointsize, options] creates a lower hemisphere equal area projection of vector data, which must be a list of dip angles and dip directions (or plunges and azimuths for lineations). It is not necessary to specify point size unless graphics options are used, in which case the point size must be given. Graphics options can be used to specify point colors, sizes, edge styles, etc."

ListEqualAreaPointColorPlot::usage="ListEqualAreaPointColorPlot[data, ColorRef, pointsize, options] creates a lower hemisphere equal area projection of vector data, which must be a list of dip angles and dip directions (or plunges and azimuths for lineations). It is not necessary to specify point size unless graphics options are used, in which case the point size must be given. Graphics options can be used to specify point colors, sizes, edge styles, etc."

ListEqualAreaArcPlot::usage="ListEqualAreaArcPlot[data, options] creates a lower hemisphere equal area projection of great circles (arcs) defined by the vector data, which must be a list of dip angles and dip directions (or plunges and azimuths for lineations). It is not necessary to specify point size unless graphics options are used, in which case the point size must be given. Graphics options can be used to specify line colors, thickness, etc."

LSAD::usage="LSAD[{\[Delta], \[Theta]}, data] calculates the Linear Sampling Angular Deviation as defined in Haneberg (2009, Envir. Engrg. Geosci. 15, 107-113), which can be used to optimize drilling directions in fractured rocks. The plunge and azimuth of a hypothetical borehole are given by {\[Delta], \[Theta]} and the fracture sets (typically mean values) are given by a list of dip vectors."

EqualAreaLSADPlot::usage="EqualAreaLSADPlot[data] creates a lower hemisphere equal area projection of contours of the LSAD (Linear Sampling Angular Deviation) function, given a list of fracture set dip vectors (typically mean orientations)."

ListEqualAreaContourPlot::usage="ListEqualAreaContourPlot[data, area, interval, options] creates a contoured equal area plot of vector data. Data are given as a list of dip vectors or the plunges and azimuths of lineations, the area is the size of the counting circle used (for example, 1% is the traditional but largely statistically unjustifiable default), and interval is the contour interval (in percent)."

ListKambPlot::usage="ListKambPlot[data, interval, options] creates a contoured equal area plot of vector data using Kamb's method, where data is a list of dip vectors or plunges/azimuths and interval is the contour interval in standard deviations."

FindDiplineClusters::usage="FindDiplineClusters[data, k] uses an angular adaptation of the Mathematica FindClusters function to objectively delineate clusters such as joint sets. The distance metric is the minimum angle between each pair of vectors, taking into account the fact that all are mapped onto the lower hemisphere. If k is not specified, Mathematica will attempt to identify the optimal number of clusters."

EqualAreaClusterPlot::usage="EqualAreaClusterPlot[data] takes the output of FindDiplineClusters and creates an equal area point plot with a different color assigned to each cluster."

FisherRandomDeviate::usage="FisherRandomDeviate[{\[Delta], \[Theta]}, \[Kappa]] generates a single random deviate from a Fisher distribution centered around {\[Delta], \[Theta]} and having the dispersion statistic \[Kappa]."

EqualAreaMesh::usage="EqualAreaMesh[ ] creates a polar mesh that can be combined with other equal area plots using Show. It takes no arguments."

ListRosePlot::usage="ListRosePlot[data_,\[CapitalDelta]\[Theta]_,\[CapitalDelta]r_,shadelevel_: GrayLevel[0.0] ] create a rose histgram for given strike"

ListStereoPointPlot;
ListStereoArcPlot;
CumFreqPlot;
CumFreqs;
KSOneList;
KSOneListPlot;
KSProb;
KSTwoList;
KSTwoListPlot;
LinePoint;
ListBoxWhiskerPlot;
ListStemPlot;
ListTernaryPlot;
RGBViewer;
SlopeAngle;
SlopeCurvature;
StereographicPlaneArc;
ThinPlateGrid;
Rainbow;
RainbowReverse;
BrownGreenCream;
BrownGreenWhite;
GreenYellowRed;
GreenWhiteRed;
RedWhiteGreen;
RedYellowGreen;
RedWhiteBlue;
BlueWhiteRed;

EqualAreaPlaneArc;
EqualAreaAngularCoords;
EqualAreaXYCoords;
CircleOverlapArea;
EucDist;
StereoXYCoords;
StereoLinePoint;
StereoPlaneArc;

EqualAreaMarklandPlot::usage="ListEqualAreaMarklandPlot[discos, {\[Delta]s, \[Theta]s}, \[Phi]] draws a Markland plot for kinematic rock slope stability analysis of wedge and plane failures. The variable discos is a list of discontinuity dip vectors, {\[Delta]s, \[Theta]s} is the dip vector of the rock slope face (assumed to be planar), and \[Phi]] is the angle of internal friction of the discontinuities."

EqualAreaFrictionCircle;


Begin["`Private`"]

colorlist = {Red, Green, Blue, Purple, Orange, Cyan, Magenta, Pink, 
   Gray};

outercircle = {GrayLevel[0.], Dashing[{0.}], Thickness[0.001], 
   Opacity[1.],
   Circle[{0., 0.}, 1.],
   Table[
    Line[{{0.98 Cos[\[Theta]], 0.98 Sin[\[Theta]]}, {Cos[\[Theta]], 
       Sin[\[Theta]]}}], {\[Theta], 0, 2 \[Pi], \[Pi]/18}], 
   Text["0", {0, 1.05}
    ],
   Text["90", {1.05, 0}, {-1, 0}],
   Text["180", {0, -1.05}],
   Text["270", {-1.05, 0}, {1, 0}],
   Options[ListStereoPlot]};

OuterCircle[] := Graphics[outercircle];

CartesianToCompass[cart_] := Module[{comp},
  comp = -cart + 90.;
  If[comp < 0., comp = comp + 360., comp];
  comp
  ]

DiplineToPole[{dip_, ddr_}] := Module[{\[Theta], \[Delta]},
  \[Delta] = 90. - dip;
  If[ddr <= 180., \[Theta] = ddr + 180., \[Theta] = ddr - 180.];
  {\[Delta], \[Theta]}
  ]

PoleToDipline[{dip_, ddr_}] := Module[{\[Theta], \[Delta]},
  \[Delta] = 90. - dip;
  If[ddr <= 180., \[Theta] = ddr + 180., \[Theta] = ddr - 180.];
  {\[Delta], \[Theta]}
  ]

DiplineToCosines[{dip_, ddr_}] := Module[{\[Theta], \[Delta]},
  \[Delta] = dip Degree;
  \[Theta] = ddr Degree;
  Chop[N[{Sin[\[Theta]] Cos[\[Delta]], 
     Cos[\[Theta]] Cos[\[Delta]], -Sin[\[Delta]]}]]
  ]

CosinesToDipline[{x_, y_, z_}] := Module[{\[Delta], \[Theta]},
  \[Delta] = -ArcSin[z]/Degree;
  \[Theta] = ArcTan[y, x]/Degree;
  If[\[Theta] < 0, \[Theta] = \[Theta] + 360];
  Return[N[{\[Delta], \[Theta]}]]
  ]

DiplineToStrike[{dip_, ddn_}] := 
 If[ddn > 90., {ddn - 90., dip}, {ddn + 270, dip}]

StrikeToDipline[{strike_, dip_}] :=
 
 If[strike < 270., {dip, strike + 90.}, {dip, strike - 270.}]

ApparentDip[{\[Delta]_, \[Alpha]_}] :=
 (* \[Delta] is dip angle, \
\[Alpha] is angle between true and apparent dip directions *)
 
 Chop[ArcTan[Tan[(\[Delta] - 0.0000001) Degree] Cos[\[Alpha] Degree]]/
   Degree]

MeanVector[indata_List] := Module[{data, cosines, len, resultant},
  data = indata;
  len = Length[data];
  cosines = Map[DiplineToCosines, data];
  resultant = Normalize[Total[cosines]];
  CosinesToDipline[resultant]
  ]

MeanVectorAlternate[indata_List] := 
 Module[{data, cosines, len, resultant, temp, maxangle},
  
  data = indata;
  len = Length[data];
  
  maxangle = 
   Max[Table[
     VectorAngle[DiplineToCosines[data[[i]]], 
      DiplineToCosines[data[[j]]]], {i, len}, {j, len}]];
  
  If[maxangle > \[Pi]/2.,
   Do[
    If[data[[i, 2]] > 180. , 
     data[[i]] = {-data[[i, 1]], data[[i, 2]] - 180.}], {i, len}]
   ];
  
  cosines = Map[DiplineToCosines, data];
  resultant = Normalize[Total[cosines]];
  temp = CosinesToDipline[resultant];
  
  If[temp[[1]] > 0, temp, {-temp[[1]], temp[[2]] + 180.}]
  ]

AngularVariance[data_List] := Module[{len, meanvec, cosines, angles},
  len = Length[data];
  meanvec = DiplineToCosines[MeanVector[data]];
  cosines = Map[DiplineToCosines, data];
  1./((len-1) Degree^2) Sum[VectorAngle[cosines[[i]],meanvec]^2,{i,1,len}] 
  ]

AngularStdDev[data_List] := Sqrt[AngularVariance[data]]

FisherK[data_List] := Module[{len, c, R},
  len = Length[data];
  c = Map[DiplineToCosines, data];
  R = Sqrt[(\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 1]\)]\)\))^2 + (\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 2]\)]\)\))^2 + (\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 3]\)]\)\))^2];
  N[(len - 1)/(len - R)]
  
  ]

AngularConfidenceRadius[data_List, p_: 0.05] := Module[{},
  len = Length[data];
  c = Map[DiplineToCosines, data];
  R = Sqrt[(\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 1]\)]\)\))^2 + (\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 2]\)]\)\))^2 + (\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 3]\)]\)\))^2];
  ArcCos[1. - ((len - R)/R) ((1./p)^(1./(len - 1)) - 1.)]/Degree
  ]

CroninBeta[data_List] := Module[{r, \[Delta], \[Theta]},
  r = AngularConfidenceRadius[data, 0.05]  \[Degree];
  \[Delta] = MeanVector[data][[1]]   \[Degree];
  \[Theta] = ArcSin[(Sin[r] Sin[\[Delta]])/(Cos[r] Cos[\[Delta]])];
  ((Sin[r] Cos[\[Theta]])/(
    Cos[r] Cos[\[Delta]] - Sin[r] Sin[\[Theta]] Sin[\[Delta]]))/
   Degree
  
  ]

EigenFabricAnalysis[data_List] := Module[{len, c, m},
  len = Length[data];
  c = Map[DiplineToCosines, data];
  m = ({
     {\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]
\*SuperscriptBox[\(c[\([i, 1]\)]\), \(2\)]\), \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 1]\)] c[\([i, 2]\)]\)\), \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 1]\)] c[\([i, 3]\)]\)\)},
     {\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 2]\)] c[\([i, 1]\)]\)\), \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]
\*SuperscriptBox[\(c[\([i, 2]\)]\), \(2\)]\), \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 2]\)] c[\([i, 3]\)]\)\)},
     {\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 3]\)] c[\([i, 1]\)]\)\), \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(c[\([i, 3]\)] c[\([i, 2]\)]\)\), \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]
\*SuperscriptBox[\(c[\([i, 3]\)]\), \(2\)]\)}
    });
  vals = Eigenvalues[m];
  vecs = Eigenvectors[m];
  
  Return[{Round[CosinesToDipline[-vecs[[1]]]], vals, vecs}]
  ]

HemisphericalVectorAngle[u_, v_] :=
 Module[{\[Alpha], \[Beta]},
  \[Alpha] = VectorAngle[u, v];
  \[Beta] = \[Pi] - VectorAngle[u, v];
  Min[{\[Alpha], \[Beta]}]
  ]

OutlierRemoval[data_List, \[Psi]_] := 
 Module[{len, c, j, t = 0, ptlist, outpts},
  (* Mahtab and Yegulalp (1982) method as described in Priest (1985) *)
\
  
  len = Length[data];
  c = 1. - Cos[\[Psi] Degree];
  outpts = Table[Null, {0}];
  
  (* Calculate number of points within radius psi of each point that \
must be exceeded in order for the point to be "dense" at a confidence \
level of 1-p. *)
  
  While[
   (1 - \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 0\), \(t\)]
\*FractionBox[\(Exp[\(-len\)\ c]\ 
\*SuperscriptBox[\((len\ c)\), \(j\)]\), \(j!\)]\)) > c,
   t++
   ];
  
  (* Loop through data and discard any points that do not satisfy the \
density criterion. Append the rest to an output list.*)
  
  Do[
   Module[{temp, npts},
    temp = 
     Table[HemisphericalVectorAngle[DiplineToCosines[data[[i]]], 
       DiplineToCosines[data[[j]]]], {j, len}];
    npts = Length[Select[temp, # < \[Psi] Degree &]];
    If[npts > t, AppendTo[outpts, data[[i]]]]
    ]
   , {i, len}
   ];
  outpts
  ]

DiplineToNormal[{dip_, ddn_}] := 
 Module[{\[Theta], \[Delta], x, y, z},
  \[Theta] = ddn Degree;
  \[Delta] = (-dip + 90.) Degree;
  x = Cos[\[Delta]] Sin[\[Theta]];
  y = Cos[\[Delta]] Cos[\[Theta]];
  z = Sin[\[Delta]];
  Chop[N[{x, y, z}]]
  
  ]

EqualAreaLinePoint[{\[Delta]_, \[Theta]_}, pointrad_] := 
 Module[{az, radius, dx, dy, graphicslist},
  az = \[Pi]/2. - \[Theta];
  If[az < 0., az = 2 \[Pi] + az];
  radius = Sqrt[2.] Sin[\[Pi]/4. - \[Delta]/2. ];
  dx = radius*Cos[az];
  dy = radius*Sin[az];
  Tooltip[
   Disk[{dx, dy}, pointrad], {Round[\[Delta]/Degree] Degree, 
    Round[\[Delta]/Degree] Degree}]
  ]

ConeSphereIntersection[\[Alpha]_] := 
 Module[{xmin, xmax, \[CapitalDelta], c, , t1, t2, \[CapitalDelta]x, 
   z, temp},
  (*
  Generates a set of {x,y,
  z} points defining the intersection of a unit sphere with a \
vertical cone of angle \[Alpha]. 
  The points can then be rotated into any position on the unit \
sphere.
  *)
  xmin = -1/(Sqrt[1 + Cot[\[Alpha]/2]^2]);
  xmax = (-xmin);
  \[CapitalDelta]x = (xmax - xmin)/200.;
  c = Tan[\[Alpha]/2.];
  z = Sqrt[xmin^2/c^2];
  temp = Sqrt[1 - x^2 - x^2 Cot[\[Alpha]/2]^2]/Sqrt[
   1 + Cot[\[Alpha]/2]^2];
  t1 = Table[{x, temp, -z}, {x, xmin + \[CapitalDelta]x, 
     xmax - \[CapitalDelta]x, \[CapitalDelta]x}];
  t2 = Table[{x, -temp, -z}, {x, xmax - \[CapitalDelta]x, 
     xmin + \[CapitalDelta]x, -\[CapitalDelta]x}];
  PrependTo[t1, t2[[Length[t2]]]];
  N[Chop[Join[t1, t2]]]
  ]

CirclePts[{dip_, ddn_}, \[Alpha]_] := Module[{pts, r},
  Off[RotationTransform::"spln"];
  pts = ConeSphereIntersection[2 \[Alpha] Degree];
  r = RotationTransform[{{0., 0., -1.}, 
     DiplineToCosines[{dip, ddn}]}, {0., 0., 0.}];
  On[RotationTransform::"spln"];
  pts = Map[CosinesToDipline, Map[r, pts]]
  ]

EqualAreaCirclePlot[{dip_, ddn_}, r_, opts___] :=
 
 Module[{pts, lowerpts, upperpts, graphicslist, styles, circlepts, 
   circlemask, \[Theta], startingpt},
  (*
  Draw a circular mask to use if the plotted circles fall \
outside of the unit circle used for the equal area net.
  *)
  
  
  circlepts = 
   Table[{Cos[\[Theta]], Sin[\[Theta]]}, {\[Theta], 0.01, 2. Pi, 
     Pi/100.}];
  circlepts = 
   AppendTo[
    circlepts, {{1.0, -10^-6}, {1.3, -10^-6}, {1.3, -1.3}, {-1.3, \
-1.3}, {-1.3, 1.3}, {1.3, 1.3}, {1.3, 10^-6}, {1.0, 10^-6}}];
  
  circlepts = Partition[Flatten[circlepts], 2];
  
  circlemask = {GrayLevel[1], Polygon[circlepts]};
  
  (*
  Calculate the points that consitute the circle to be projected \
onto the unit sphere.
  *)
  
  lowerpts = Map[EqualAreaXYCoords, CirclePts[{dip  , ddn  }, r  ] ];
  
  upperpts = Select[lowerpts, Sqrt[#[[1]]^2 + #[[2]]^2] > 1 &];
  upperpts = Map[EqualAreaAngularCoords, upperpts];
   
  upperpts = 
   Table[{-upperpts[[i, 1]], Mod[upperpts[[i, 2]] + 180., 360.]}, {i, 
     Length[upperpts]}];
  
  upperpts = Map[EqualAreaXYCoords, upperpts];
  
  
  Off[ReplaceAll::argr];
  styles = Style /. opts;
  On[ReplaceAll::argr];
  
  graphicslist = {Black, EdgeForm[Black]};
  
  If[
   Length[styles] > 1, AppendTo[graphicslist, styles]
   ];
  
  AppendTo[graphicslist, Line[lowerpts]];
  
  AppendTo[graphicslist, {Line[upperpts]}];
  
  AppendTo[
   graphicslist, {White, EdgeForm[], Opacity[1], outercircle, 
    circlemask}
   ];
  
  
  Show[Graphics[Flatten[graphicslist]], 
   PlotRange -> {{-1.1, 1.1}, {-1.1, 1.1}}, 
   PlotRangeClipping -> True]
  ]



Options[CumFreqPlot] = {Axes->True,Frame->False};
CumFreqPlot[indata_,minval_,maxval_,opts___]:=Block[{len,dz,temp,temp2=indata},
len=Length[temp2];	
dz=1./(len);	
data=CumFreqs[temp2];
temp=Table[{Line[{{data[[i,1]],data[[i,2]]-dz},data[[i]]}],Line[{data[[i]],{data[[i+1,1]],data[[i,2]]}}]},{i,len-1}];
Show[Graphics[{Line[{{minval,0},{data[[1,1]],0}}],temp,Line[{{data[[len,1]],data[[len,2]]-dz},data[[len]]}],Line[{data[[len]],{maxval,1}}]}],opts,Options[CumFreqPlot],PlotRange->{{minval,maxval},{-0.05,1.05}}]]


CumFreqs[indata_List]:=Block[{len,sorted,cumfreqlist,i,data=indata},	
	data;
	len=Length[data];	
		cumfreqlist=Table[{0,0},{len}];	
	sorted=Sort[data];	
	Do[
		Block[{},
			cumfreqlist[[i,1]]=sorted[[i]];		
			cumfreqlist[[i,2]]=i/len;
		]	,{i,len}	
	];
N[cumfreqlist]
]


KSOneList[data_] := Block[{temp,temp2,temp3,datalen,i,meanval,sdev,cumfreqs,minval,maxval,obs,calc,dz,KS},
datalen=Length[data];	
meanval=Mean[data];	
sdev=StandardDeviation[data];	
cumfreqs=CumFreqs[data];	
obs=Table[cumfreqs[[i,2]],{i,datalen}];	
calc=Table[CDF[NormalDistribution[meanval,sdev],cumfreqs[[i,1]]],{i,datalen}];
	KS=Max[Table[Abs[obs[[i]]-calc[[i]]],{i,Length[obs]}]];	
	Return[KS]
]


Options[KSOneListPlot]={Frame->False,Axes->True};
KSOneListPlot[data_,minplot_,maxplot_,opts___] := Block[{temp,temp2,temp3,datalen,i,meanval,sdev,cumfreqs,minval,maxval,obs,calc,dz,KS},	
	datalen=Length[data];		
	meanval=Mean[data];	
	sdev=StandardDeviation[data];	
	cumfreqs=CumFreqs[data];	
	obs=Table[cumfreqs[[i,2]],{i,datalen}];	
	calc=Table[CDF[NormalDistribution[meanval,sdev],cumfreqs[[i,1]]],{i,datalen}];
	dz=1./datalen;	temp2=Table[{Line[{{cumfreqs[[i,1]],cumfreqs[[i,2]]-dz},cumfreqs[[i]]}],		Line[{cumfreqs[[i]],{cumfreqs[[i+1,1]],cumfreqs[[i,2]]}}]},{i,datalen-1}];	temp3=Graphics[{Line[{{minplot,0},{cumfreqs[[1,1]],0}}], 		temp2,		Line[{{cumfreqs[[datalen,1]],cumfreqs[[datalen,2]]-dz},cumfreqs[[datalen]]}],		Line[{cumfreqs[[datalen]],{maxplot,1}}]}];	calcplot=Plot[CDF[NormalDistribution[meanval,sdev],x],{x,minplot,maxplot},PlotStyle->{Dashing[{0.02,0.02}]},DisplayFunction->Identity,PlotRange->{{minplot,maxplot},{0,1}}];	Show[temp3,calcplot,DisplayFunction->$DisplayFunction,PlotRange->{{minplot,maxplot},{-0.05,1.05}},opts,Options[KSOneListPlot]]
]


KSProb[KS_,len_]:= Module[{QKS,\[Lambda]},
\[Lambda] =( Sqrt[len]+0.12 + 0.11/Sqrt[len])*KS;
Return[2. \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 1\), \(100\)]\((
\*SuperscriptBox[\((\(-1\))\), \(j - 1\)] Exp[\(-2\) 
\*SuperscriptBox[\(j\), \(2\)]\ 
\*SuperscriptBox[\(\[Lambda]\), \(2\)]])\)\)]
]


KSTwoList[data1__List,data2__List,npts_] :=

Block[{freqs1,freqs2,\[CapitalDelta]x,nfreqs,nvals,temp,minval,maxval},

minval=Min[Append[data1,data2]];
maxval=Max[Append[data1,data2]];

freqs1=CumFreqs[data1];
freqs2= CumFreqs[data2];

freqs1=Append[Prepend[freqs1,{minval,0.}],{maxval,1.}];
freqs2=Append[Prepend[freqs2,{minval,0.}],{maxval,1.}];

nfreqs=Length[freqs1]-1;
nvals = npts+1;

\[CapitalDelta]x = (maxval-minval)/npts;
xvals=Table[minval+(i-1)*\[CapitalDelta]x ,{i,nvals}];

yvals1=Table[0.,{nvals}];
yvals2=Table[0.,{nvals}];

Do[
Do[If[xvals[[i]] >= freqs1[[j,1]] && xvals[[i]] < freqs1[[j+1,1]],yvals1[[i]] = freqs1[[j,2]]],{j,1,nfreqs}
],{i,nvals}
];

Do[
Do[If[xvals[[i]] >= freqs2[[j,1]] && xvals[[i]] < freqs2[[j+1,1]],yvals2[[i]] = freqs2[[j,2]]],{j,1,nfreqs}
],{i,nvals}
];

KS=Max[Abs[yvals1-yvals2]];
Return[KS]
]


Options[KSTwoListPlot]={Frame->False,Axes->True};
KSTwoListPlot[data1_,data2_, minplot_,maxplot_,opts___] := Block[{temp,temp2,temp3,temp4,temp5,datalen,i,cumfreqs1, cumfreqs2,minval,maxval,obs1, obs2,dz,KS},	
	datalen=Length[data1];			
	cumfreqs1=CumFreqs[data1];
cumfreqs2=CumFreqs[data2];	
	obs1=Table[cumfreqs1[[i,2]],{i,datalen}];	
	obs2=Table[cumfreqs2[[i,2]],{i,datalen}];
dz=1./datalen;
	
temp2=Table[{Line[{{cumfreqs1[[i,1]],cumfreqs1[[i,2]]-dz},cumfreqs1[[i]]}],		Line[{cumfreqs1[[i]],{cumfreqs1[[i+1,1]],cumfreqs1[[i,2]]}}]},{i,datalen-1}];	

temp3=Graphics[{Line[{{minplot,0},{cumfreqs1[[1,1]],0}}], 		temp2,		Line[{{cumfreqs1[[datalen,1]],cumfreqs1[[datalen,2]]-dz},cumfreqs1[[datalen]]}],		Line[{cumfreqs1[[datalen]],{maxplot,1}}]}];

temp4=Table[{Line[{{cumfreqs2[[i,1]],cumfreqs2[[i,2]]-dz},cumfreqs2[[i]]}],		Line[{cumfreqs2[[i]],{cumfreqs2[[i+1,1]],cumfreqs2[[i,2]]}}]},{i,datalen-1}];	

temp5=Graphics[{Dashing[{0.01}],Line[{{minplot,0},{cumfreqs2[[1,1]],0}}], 		temp4,		Line[{{cumfreqs2[[datalen,1]],cumfreqs2[[datalen,2]]-dz},cumfreqs2[[datalen]]}],		Line[{cumfreqs2[[datalen]],{maxplot,1}}]}];

Show[temp3,temp5, DisplayFunction->$DisplayFunction,PlotRange->{{minplot,maxplot},{-0.05,1.05}},opts,Options[KSTwoListPlot]]
]


LinePoint[plunge_,azimuth_,pointrad_,pointflag_]:= 
Module[{radius,dx,dy},
az = \[Pi]/2. - azimuth;
radius=Tan[\[Pi]/4.-plunge/2. ];
dx = radius*Cos[az];
dy = radius*Sin[az];
If[pointflag=="Open"|| pointflag=="open",Circle[{dx,dy},pointrad],Disk[{dx,dy},pointrad]]
]


Options[ListBoxWhiskerPlot]={AspectRatio->1,Axes->None,Frame->True};
ListBoxWhiskerPlot[data_,del_:0.1,opts___]:=Block[{i},results=Table[{},{i,Length[data]}];
results=Graphics[{GrayLevel[0],Table[{Line[{{i,data[[i,4]]},{i,data[[i,5]]}}],Line[{{i,data[[i,1]]},{i,data[[i,2]]}}],Line[{{i-del/2,data[[i,1]]},{i+del/2,data[[i,1]]}}],Line[{{i-del/2,data[[i,5]]},{i+del/2,data[[i,5]]}}],Line[{{i-del,data[[i,2]]},{i+del,data[[i,2]]},{i+del,data[[i,4]]},{i-del,data[[i,4]]},{i-del,data[[i,2]]}}],Line[{{i-del,data[[i,3]]},{i+del,data[[i,3]]}}]},{i,1,Length[data]}]}];
Show[results,PlotRange->{{0,Length[data]+1},All},FrameTicks->{Table[i,{i,1,Length[data]}],Automatic},opts,Options[ListBoxWhiskerPlot]]]


Options[ListStemPlot] = {Axes->True,Frame->False};
ListStemPlot[indata_,ballsize_,opts___]:=
Block[{stems,balls,len},
len=Length[indata];
If[Length[Dimensions[indata]] == 1,
Block[{},
stems=Table[Line[{{i,0},{i,indata[[i]]}}],{i,len}];
balls=Table[Point[{i,indata[[i]]}],{i,len}];
Show[Graphics[{stems,PointSize[ballsize],balls}],opts,Options[ListStemPlot]];
],
Block[{},stems=Table[Line[{{indata[[i,1]],0},{indata[[i,1]],indata[[i,2]]}}],{i,len}];
balls=Table[Point[{indata[[i,1]],indata[[i,2]]}],{i,len}];
Show[Graphics[{stems,PointSize[ballsize],balls}],opts,Options[ListStemPlot]];
]
]
]


ListTernaryPlot[data_,labels_,pointsize_:0.02,pointshade_: GrayLevel[0.]]:=
Module[{tan60,tan30,cos60,sin60,\[CapitalDelta]x,triangle,grid1,grid2,grid3,baseticks,i,datalen},
tan60 = Tan[60. Degree];
tan30 = Tan[30. Degree];
cos60 = Cos[60. Degree];
sin60 = Sin[60. Degree];
\[CapitalDelta]x = 1./tan60;
datalen=Length[data];
labelA=ToString[labels[[1]]];
labelB = ToString[labels[[2]]];
labelC = ToString[labels[[3]]];
triangle=
Line[{{0.,1.},{\[CapitalDelta]x,0.},{-\[CapitalDelta]x ,0.},{0,1}}];
baseticks=Table[-\[CapitalDelta]x +i (2 \[CapitalDelta]x)/10. ,{i,1,9}];
gridlen=Table[2 (1-j) tan30,{j,0.1,0.9,0.1}];
grid1=Table[Line[{{-gridlen[[i]]/2.,i/10.},{gridlen[[i]]/2.,i/10.}}],{i,1,9}];
grid2=Table[Line[{{baseticks[[i]],0.},{baseticks[[i]]+gridlen[[i]]cos60,gridlen[[i]] sin60}}],{i,1,9}];
grid3=Table[Line[{{baseticks[[10-i]],0.},{baseticks[[10-i]]-gridlen[[i]]cos60,gridlen[[i]] sin60}}],{i,1,9}];
datapts=Table[Disk[{((data[[i,1]]-1.)cos60  + data[[i,2]])/(cos60*tan60),data[[i,1]]},pointsize],{i,datalen}];
Show[
Graphics[
{triangle,Dashing[{0.007}],grid1,grid2,grid3,Dashing[{1}],pointshade,datapts,GrayLevel[0],Text[labelA,{0.,1.05}],Text[labelB,{0.6,0.},{-1,0}],Text[labelC,{-0.6,0.},{1,0}]}
],AspectRatio->0.75,Axes->False,PlotRange->{{-0.8,0.8},{-0.05,1.1}}
]
]


RGBViewer[r_,g_,b_]:= Show[Graphics[{RGBColor[r,g,b],Rectangle[{0,0},{1,1}]}]]


SlopeAngle[indata_,spacing_]:=
Block[
{r,c,nr,nc,results,bvals,data,i,temp,c1},
bvals = Table[0,{i,4}];
data=indata;
nr=Length[data];
nc=Length[data[[1]]];
dx2=2*spacing;
temp=Table[0,{nr+2},{nc+2}];



bvals[[1]] = 2.*data[[nr]]-data[[nr-1]];


bvals[[2]] = Table[2.*data[[r,nc]]-data[[r,nc-1]],{r,1,nr}];

bvals[[3]] = 2.*data[[1]]-data[[2]];

bvals[[4]] = Table[2.*data[[r,2]]-data[[r,3]],{r,1,nr}];

Do[temp[[1,c+1]]=bvals[[3,c]],{c,1,nc}];
Do[temp[[nr+2,c+1]]=bvals[[1,c]],{c,1,nc}];
Do[temp[[r+1,nc+2]]=bvals[[2,r]],{r,1,nr}];
Do[temp[[r+1,1]]=bvals[[4,r]],{r,1,nr}];
Do[temp[[r+1,c+1]] = data[[r,c]],{r,1,nr},{c,1,nc}];

Return[
Table[ArcTan[Sqrt[(temp[[r+1,c]]-temp[[r-1,c]])^2 + (temp[[r,c+1]]-temp[[r,c-1]])^2]/(dx2)]/Degree,
{r,2,nr+1},{c,2,nc+1}]
]; 
Clear[temp]
]


SlopeCurvature[indata_,spacing_]:=
Block[
{r,c,nr,nc,results,bvals,data,i,temp,c1},
bvals = Table[0,{i,4}];
data=indata;
nr=Length[data];
nc=Length[data[[1]]];
dx=spacing;
temp=Table[0,{nr+2},{nc+2}];
c1=90./dx/Pi;

bvals[[1]] = 2.*data[[nr]]-data[[nr-1]];

bvals[[2]] = Table[2.*data[[r,nc]]-data[[r,nc-1]],{r,1,nr}];

bvals[[3]] = 2.*data[[1]]-data[[2]];

bvals[[4]] = Table[2.*data[[r,2]]-data[[r,3]],{r,1,nr}];


Do[temp[[1,c+1]]=bvals[[3,c]],{c,1,nc}];
Do[temp[[nr+2,c+1]]=bvals[[1,c]],{c,1,nc}];
Do[temp[[r+1,nc+2]]=bvals[[2,r]],{r,1,nr}];
Do[temp[[r+1,1]]=bvals[[4,r]],{r,1,nr}];
Do[temp[[r+1,c+1]] = data[[r,c]],{r,1,nr},{c,1,nc}];

Return[
Table[(temp[[r,c+1]]+temp[[r,c-1]]+temp[[r-1,c]]+temp[[r+1,c]]-4.*temp[[r,c]])/dx/dx,
{r,2,nr+1},{c,2,nc+1}]
]; 
Clear[temp]
]


StereographicPlaneArc[strike_,dip_]:= 
Module[{r,dx,dy,dipazimuth},
dx = Tan[dip]*Cos[-\[Pi] -strike];
dy = Tan[dip]*Sin[-\[Pi] -strike];
dipazimuth = -strike ;
radius=Tan[dip ] + Tan[\[Pi]/4.-dip/2 ];
Circle[{dx,dy},radius,{dipazimuth - (\[Pi]/2.-dip),dipazimuth + (\[Pi]/2.-dip)}]
]


ThinPlateGrid[indata_,xvals_,yvals_,dx_,tolconst_]:=
Block[{len,xmin,xmax,nx,ymin,ymax,ny,temp,fixedpts,fixedptlocs},
xmin=xvals[[1]];
xmax=xvals[[2]];

ymin=yvals[[1]];
ymax=yvals[[2]];

len=Length[indata];
fixedpts = indata;
fixedptlocs = Table[{0,0},{len}];

nc =Round[ (xmax-xmin)/dx]+5;
nr =Round[ (ymax-ymin)/dx]+5;

new=Table[0.,{nr},{nc}];

Do[
Block[{},
fixedpts[[k,1]]=Round[(indata[[k,1]]-xmin)/dx+3];
fixedpts[[k,2]]=Round[(indata[[k,2]]-ymin)/dx+3];
fixedptlocs[[k]] = {fixedpts[[k,1]],fixedpts[[k,2]]};


new[[fixedpts[[k,2]],fixedpts[[k,1]]]] = indata[[k,3]]
],{k,len}
];

maxerr=10^6;
temp=Table[fixedpts[[i,3]],{i,len}];
toler=tolconst (Max[temp]-Min[temp]);
ct = 0;

new[[2]] = new[[3]];
new[[1]] = new[[3]];

new[[nr-1]] = new[[nr-2]];
new[[nr]] = new[[nr-2]];

Do[
Block[{},
new[[r,2]]=new[[r,3]];
new[[r,1]]= new[[r,3]];
new[[r,nc-1]]=new[[r,nc-2]];
new[[r,nc]]=new[[r,nc-2]]],{r,nr}
];

While[maxerr>toler,
Block[{k},
old=new;
maxerr = 0.;
Do[
Block[{z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12},
z1=old[[r,c+1]];
z2= old[[r+1,c]];
z3=new[[r,c-1]];
z4=new[[r-1,c]];
z5=old[[r,c+2]];
z6=old[[r+1,c+1]];
z7=old[[r+2,c]];
z8=old[[r+1,c-1]];
z9=new[[r,c-2]];
z10 = new[[r-1,c-1]];
z11 = new[[r-2,c]];
z12 = new[[r-1,c+1]];

If[Length[{{c,r}} \[Intersection] fixedptlocs] == 0,
new[[r,c]]= 1/20. (8(z1+z2+z3+z4) - 2 (z6+z8+z10+z12) -1 (z5 + z7+ z9 + z11))
]
],{r,3,nr-2},{c,3,nc-2}
];

Do[
new[[fixedpts[[i,2]],fixedpts[[i,1]]]] = fixedpts[[i,3]],{i,len}
];

new[[2]] = new[[3]];
new[[1]] = new[[3]];

new[[nr-1]] = new[[nr-2]];
new[[nr]] = new[[nr-2]];

Do[
Block[{},
new[[r,2]]=new[[r,3]];
new[[r,1]]= new[[r,3]];
new[[r,nc-1]]=new[[r,nc-2]];
new[[r,nc]]=new[[r,nc-2]]],{r,nr}
];

maxerr=Max[Abs[new-old]];
ct++;
]
];
Return[Table[new[[r,c]],{r,3,nr-2},{c,3,nc-2}]]
]


Rainbow[z_]:=Module[{RedCurve,GreenCurve,BlueCurve},
RedCurve=Interpolation[{{0.,1.},{0.17,1.},{0.34,1.},{0.50, 0.},{0.67,0.},{0.85,0.03},{1.,0.56}},InterpolationOrder->1];
GreenCurve=Interpolation[{{0.,0.},{0.17,0.05},{0.34,1.},{0.50, 1.},{0.67,0.},{0.85,0.18},{1.,0.37}},InterpolationOrder->1];
BlueCurve=Interpolation[{{0.,0.},{0.17,0.0},{0.34,0.},{0.50, 0.},{0.67,1.},{0.85,0.33},{1.,0.60}},InterpolationOrder->1];
RGBColor[RedCurve[z],GreenCurve[z],BlueCurve[z]]
]


RainbowReverse[z_]:=Module[{RedCurve,GreenCurve,BlueCurve},
RedCurve=Interpolation[{{1.,1.},{0.85,1.},{0.67,1.},{0.50, 0.},{0.34,0.},{0.17,0.03},{0.,0.56}},InterpolationOrder->1];
GreenCurve=Interpolation[{{1.,0.},{0.85,0.05},{0.67,1.},{0.50, 1.},{0.34,0.},{0.17,0.18},{0.,0.37}},InterpolationOrder->1];
BlueCurve=Interpolation[{{1.,0.},{0.85,0.0},{0.67,0.},{0.50, 0.},{0.34,1.},{0.17,0.33},{0.,0.60}},InterpolationOrder->1];
RGBColor[RedCurve[z],GreenCurve[z],BlueCurve[z]]
]


BrownGreenCream[z_]:=Module[{RedCurve,GreenCurve,BlueCurve},
RedCurve=Interpolation[{{0.,0.63},{0.5,0.04},{1.,1.}},InterpolationOrder->1];
GreenCurve=Interpolation[{{0.,0.32},{0.5,0.79},{1.,0.97}},InterpolationOrder->1];
BlueCurve=Interpolation[{{0.,0.18},{0.5,0.17},{1.,0.86}},InterpolationOrder->1];
RGBColor[RedCurve[z],GreenCurve[z],BlueCurve[z]]
];

BrownGreenWhite[z_]:=Module[{RedCurve,GreenCurve,BlueCurve},
RedCurve=Interpolation[{{0.,0.63},{0.5,0.04},{0.9,1.},{1.,1.}},InterpolationOrder->1];
GreenCurve=Interpolation[{{0.,0.32},{0.5,0.79},{0.9,0.97},{1.,1.}},InterpolationOrder->1];
BlueCurve=Interpolation[{{0.,0.18},{0.5,0.17},{0.9,0.86},{1.,1.}},InterpolationOrder->1];
RGBColor[RedCurve[z],GreenCurve[z],BlueCurve[z]]
];

GreenYellowRed[z_]:=Module[{RedCurve,GreenCurve,BlueCurve},
RedCurve=Interpolation[{{0.,0.},{0.5,1.},{1.,1.}},InterpolationOrder->1];
GreenCurve=Interpolation[{{0.,1.},{0.5,1.},{1.,0.}},InterpolationOrder->1];
RGBColor[RedCurve[z],GreenCurve[z],0.]
];

GreenWhiteRed[z_]:=Module[{RedCurve,GreenCurve,BlueCurve},
RedCurve=Interpolation[{{0.,0.},{0.5,1.},{1.,1.}},InterpolationOrder->1];
GreenCurve=Interpolation[{{0.,1.},{0.5,1.},{1.,0.}},InterpolationOrder->1];
BlueCurve=Interpolation[{{0.,0.0},{0.5,1.},{1.,0.}},InterpolationOrder->1];
RGBColor[RedCurve[z],GreenCurve[z],BlueCurve[z]]
];

RedWhiteGreen[z_]:=Module[{RedCurve,GreenCurve,BlueCurve},
GreenCurve=Interpolation[{{0.,0.},{0.5,1.},{1.,1.}},InterpolationOrder->1];
RedCurve=Interpolation[{{0.,1.},{0.5,1.},{1.,0.}},InterpolationOrder->1];
BlueCurve=Interpolation[{{0.,0.0},{0.5,1.},{1.,0.}},InterpolationOrder->1];
RGBColor[RedCurve[z],GreenCurve[z],BlueCurve[z]]
];

RedYellowGreen[z_]:=Module[{RedCurve,GreenCurve,BlueCurve},
GreenCurve=Interpolation[{{0.,0.},{0.5,1.},{1.,1.}},InterpolationOrder->1];
RedCurve=Interpolation[{{0.,1.},{0.5,1.},{1.,0.}},InterpolationOrder->1];
RGBColor[RedCurve[z],GreenCurve[z],0.]
];

RedWhiteBlue[z_]:=Module[{RedCurve,GreenCurve,BlueCurve},
RedCurve=Interpolation[{{0.,1.},{0.5,1.},{1.,0.}},InterpolationOrder->1];
GreenCurve=Interpolation[{{0.,0.},{0.5,1.},{1.,0.}},InterpolationOrder->1];
BlueCurve=Interpolation[{{0.,0.0},{0.5,1.},{1.,1.}},InterpolationOrder->1];
RGBColor[RedCurve[z],GreenCurve[z],BlueCurve[z]]
];

BlueWhiteRed[z_]:=Module[{RedCurve,GreenCurve,BlueCurve},
BlueCurve=Interpolation[{{0.,1.},{0.5,1.},{1.,0.}},InterpolationOrder->1];
GreenCurve=Interpolation[{{0.,0.},{0.5,1.},{1.,0.}},InterpolationOrder->1];
RedCurve=Interpolation[{{0.,0.0},{0.5,1.},{1.,1.}},InterpolationOrder->1];
RGBColor[RedCurve[z],GreenCurve[z],BlueCurve[z]]
];


ListRosePlot[data_,\[CapitalDelta]\[Theta]_,\[CapitalDelta]r_,shadelevel_: GrayLevel[0.0] ]:=
Module[{r,binned,len,lenbinned,maxdat,maxbinned,flag},
len=Length[data];

If[Length[Select[data,(# > 90. && # < 270.)&]] == 0,flag="bi",flag="uni"];

binned=BinCounts[data,{0.,360.,N[\[CapitalDelta]\[Theta]/Degree]}];	
lenbinned=Length[binned];
maxdat=Max[data];
maxbinned=Max[binned];
maxrad=Ceiling[maxbinned/\[CapitalDelta]r]*\[CapitalDelta]r;

If[flag == "bi",
Show[
Graphics[{shadelevel,
Table[Disk[{0.,0.},binned[[i]],{  \[Pi]/2-i \[CapitalDelta]\[Theta] ,\[Pi]/2-(i-1) \[CapitalDelta]\[Theta]}],{i,lenbinned}],
Table[Disk[{0.,0.},binned[[i]],{ - (\[Pi]/2)-i \[CapitalDelta]\[Theta] ,-(\[Pi]/2)-(i-1) \[CapitalDelta]\[Theta] }],{i,lenbinned}],GrayLevel[0.0],
Table[Circle[{0.,0.},r],{r,\[CapitalDelta]r,maxrad,\[CapitalDelta]r}],
Line[{{-maxrad,0},{maxrad,0}}],
Line[{{0,-maxrad},{0,maxrad}}],Line[{{0.50 maxrad,0.86 maxrad},{-0.50 maxrad,-0.86 maxrad}}],
Line[{{-0.50 maxrad,0.86 maxrad},{0.50 maxrad,-0.86 maxrad}}],Line[{{0.86 maxrad,0.50 maxrad},{-0.86maxrad,-0.50 maxrad}}],
Line[{{-0.86 maxrad,0.50 maxrad},{0.86 maxrad,-0.50 maxrad}}],
Text["0",{0,1.05 maxrad}],
Text["90",{1.05 maxrad,0},{-1,0}],
Text["180",{0,-1.05 maxrad}],
Text["270",{-1.05 maxrad,0},{1,0}]}
], AspectRatio->1.,PlotRange->{{-1.2 maxrad,1.2 maxrad},{-1.2 maxrad,1.2 maxrad}}
],
Show[
Graphics[{shadelevel,
Table[Disk[{0.,0.},binned[[i]],{  \[Pi]/2-i \[CapitalDelta]\[Theta] ,\[Pi]/2-(i-1) \[CapitalDelta]\[Theta]}],{i,lenbinned}],GrayLevel[0.0],
Table[Circle[{0.,0.},r],{r,\[CapitalDelta]r,maxrad,\[CapitalDelta]r}],
Line[{{-maxrad,0},{maxrad,0}}],
Line[{{0,-maxrad},{0,maxrad}}],Line[{{0.50 maxrad,0.86 maxrad},{-0.50 maxrad,-0.86 maxrad}}],
Line[{{-0.50 maxrad,0.86 maxrad},{0.50 maxrad,-0.86 maxrad}}],Line[{{0.86 maxrad,0.50 maxrad},{-0.86maxrad,-0.50 maxrad}}],
Line[{{-0.86 maxrad,0.50 maxrad},{0.86 maxrad,-0.50 maxrad}}],
Text["0",{0,1.05 maxrad}],
Text["90",{1.05 maxrad,0},{-1,0}],
Text["180",{0,-1.05 maxrad}],
Text["270",{-1.05 maxrad,0},{1,0}]}
],AspectRatio->1.,PlotRange->{{-1.2 maxrad,1.2 maxrad},{-1.2 maxrad,1.2 maxrad}}
]]
]


ListEqualAreaPointPlot[indata_, pointrad_: 0.02, opts___] :=
 
 Module[{data, len, graphicslist, styles, fontsize},
  
  data = indata Degree;
  len = Length[data];
  
  Off[ReplaceAll::argr];
  styles = Style /. opts;
  On[ReplaceAll::argr];
  
  fontsize = 0.025;
  
  graphicslist = {};
  
  If[
   Length[styles] > 1, AppendTo[graphicslist, styles]
   ];
  
  AppendTo[graphicslist, {Table[
     EqualAreaLinePoint[data[[i]] , pointrad], {i, len}
     ],
    Black, Opacity[1.],
    Circle[{0., 0.}, 1.],
    Table[
     Line[{{0.98 Cos[\[Theta]], 0.98 Sin[\[Theta]]}, {Cos[\[Theta]], 
        Sin[\[Theta]]}}], {\[Theta], 0, 2 \[Pi], \[Pi]/18}
     ],
    Text[Style["0", FontSize -> Scaled[fontsize]], {0, 
      1 + 2 fontsize}],
    Text[Style["90", FontSize -> Scaled[fontsize]], {1 + 2 fontsize, 
      0}],
    Text[Style["180", 
      FontSize -> Scaled[fontsize]], {0, -1 - 2 fontsize}],
    Text[Style["270", FontSize -> Scaled[fontsize]], {-1 - 2 fontsize,
       0}],
    Options[ListStereoPlot]}
   ];
  
  Graphics[Flatten[graphicslist]]
  
  ]

ListEqualAreaPointColorPlot[indata_, colorref_, pointrad_: 0.02, opts___] :=
 
 Module[{data, len, graphicslist, styles, fontsize},
  
  data = indata Degree;
  len = Length[data];
  
  Off[ReplaceAll::argr];
  styles = Style /. opts;
  On[ReplaceAll::argr];
  
  fontsize = 0.025;
  
  graphicslist = {};
  
  If[
   Length[styles] > 1, AppendTo[graphicslist, styles]
   ];
  
  AppendTo[graphicslist, {
    Black, Opacity[1.],
    Circle[{0., 0.}, 1.],
    Table[
     Line[{{0.98 Cos[\[Theta]], 0.98 Sin[\[Theta]]}, {Cos[\[Theta]], 
        Sin[\[Theta]]}}], {\[Theta], 0, 2 \[Pi], \[Pi]/18}
     ],
    Text[Style["0", FontSize -> Scaled[fontsize]], {0, 
      1 + 2 fontsize}],
    Text[Style["90", FontSize -> Scaled[fontsize]], {1 + 2 fontsize, 
      0}],
    Text[Style["180", 
      FontSize -> Scaled[fontsize]], {0, -1 - 2 fontsize}],
    Text[Style["270", FontSize -> Scaled[fontsize]], {-1 - 2 fontsize,
       0}],
    Options[ListStereoPlot]}
   ];
{Graphics[Table[
az = \[Pi]/2. - data[[i]][[2]];
  If[az < 0., az = 2 \[Pi] + az];
  radius = Sqrt[2.] Sin[\[Pi]/4. - data[[i]][[1]]/2. ];
  dx = radius*Cos[az];
  dy = radius*Sin[az];
colorin=Rescale[colorref][[i]][[1]];
  {ColorData["Rainbow"][colorin],Tooltip[
   Disk[{dx, dy}, pointrad], {Round[data[[i]][[1]]/Degree] Degree, 
    Round[data[[i]][[1]]/Degree] Degree}]}, {i, len}]],
  Graphics[Flatten[graphicslist]]}
  
  ]

ListEqualAreaArcPlot[indata_, opts___] :=
 
 Module[{data, len, graphicslist, styles, fontsize},
  
  data = indata;
  
  len = Length[data];
  
  Off[ReplaceAll::argr];
  styles = Style /. opts;
  On[ReplaceAll::argr];
  
  fontsize = 0.025;
  
  graphicslist = {};
  
  If[
   Length[styles] >= 1, AppendTo[graphicslist, styles]
   ];
  
  AppendTo[graphicslist, {Table[
     EqualAreaPlaneArc[data[[i]]], {i, len}
     ],
    Thickness[Small],
    Black, Opacity[1.],
    Circle[{0., 0.}, 1.],
    Table[
     Line[{{0.98 Cos[\[Theta]], 0.98 Sin[\[Theta]]}, {Cos[\[Theta]], 
        Sin[\[Theta]]}}], {\[Theta], 0, 2 \[Pi], \[Pi]/18}
     ],
    Text[Style["0", FontSize -> Scaled[fontsize]], {0, 
      1 + 2 fontsize}],
    Text[Style["90", FontSize -> Scaled[fontsize]], {1 + 2 fontsize, 
      0}],
    Text[Style["180", 
      FontSize -> Scaled[fontsize]], {0, -1 - 2 fontsize}],
    Text[Style["270", FontSize -> Scaled[fontsize]], {-1 - 2 fontsize,
       0}],
    Options[ListStereoPlot]}
   ];
  
  Graphics[Flatten[graphicslist]]
  
  ]

LSAD[{\[Delta]_, \[Theta]_}, discos_] := 
 Module[{borehole, discosines, angles, angvar},
  len = Length[discos];
  discosines = Map[DiplineToNormal, discos];
  borehole = DiplineToCosines[{\[Delta], \[Theta]}];
  angles = 
   Table[Min[{VectorAngle[borehole, -discosines[[i]]]^2, 
      VectorAngle[borehole, discosines[[i]]]^2}], {i, len}];
  
  angvar = 1/(len Degree^2) \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(len\)]\(angles[\([i]\)]\)\);
  Sqrt[angvar]
  ]

EqualAreaLSADPlot[discos_List, opts___] := 
 Module[{cplot, ctable, circlepts, circlemask, otherstuff, 
   graphicslist, fontsize},
  
  ctable = Table[c, {c, 0, 90, 2}];
  
  cplot = ListContourPlot[
    Table[
     LSAD[EqualAreaAngularCoords[{x, y}], discos], {x, -1, 1, 
      0.03}, {y, -1, 1, 0.03}],
    DataRange -> {{-1, 1}, {-1, 1}}, Contours -> ctable, 
    Frame -> False, PlotRange -> All, DisplayFunction -> Identity, 
    Axes -> False, 
    CoordinatesToolOptions -> {"DisplayFunction" -> 
       Function[pt, Round[EqualAreaAngularCoords[pt]] Degree] }, 
    Options[ListStereoPlot], opts
    ];
  
  circlepts = 
   Table[{Cos[\[Theta]], Sin[\[Theta]]}, {\[Theta], 0.01, 
     2. Pi - 0.01, Pi/100.}];
  
  circlepts = 
   AppendTo[
    circlepts, {{1.0, -10^-6}, {1.3, -10^-6}, {1.3, -1.3}, {-1.3, \
-1.3}, {-1.3, 1.3}, {1.3, 1.3}, {1.3, 10^-6}, {1.0, 10^-6}}];
  
  circlepts = Partition[Flatten[circlepts], 2];
  
  circlemask = {GrayLevel[1], Polygon[circlepts]};
   graphicslist = {circlemask};
  
  (*
  Draw a circular mask to use if the plotted circles fall \
outside of the unit circle used for the equal area net.
  *)
  
  fontsize = 0.025;
  
  otherstuff = {White, EdgeForm[], Opacity[1],
    Black, Dashing[{0.}], Thickness[0.001], Opacity[1.],
    Circle[{0., 0.}, 1.],
    Table[
     Line[{{0.98 Cos[\[Theta]], 0.98 Sin[\[Theta]]}, {Cos[\[Theta]], 
        Sin[\[Theta]]}}], {\[Theta], 0, 2 \[Pi], \[Pi]/18}], 
    Text[Style["0", FontSize -> Scaled[fontsize]], {0, 
      1 + 2 fontsize}],
    Text[Style["90", FontSize -> Scaled[fontsize]], {1 + 2 fontsize, 
      0}],
    Text[Style["180", 
      FontSize -> Scaled[fontsize]], {0, -1 - 2 fontsize}],
    Text[Style["270", FontSize -> Scaled[fontsize]], {-1 - 2 fontsize,
       0}]};
  
  AppendTo[graphicslist, otherstuff];
  ;
  
  
  Show[cplot, Graphics[Flatten[graphicslist]]]
  
  ]

FindDiplineClusters[data_, nclusters___] :=
 Map[
  CosinesToDipline, 
  FindClusters[Map[DiplineToCosines, data], nclusters, 
   DistanceFunction -> HemisphericalVectorAngle], {2}
  ]

EqualAreaClusterPlot[clusteredpoints_List] := Show[
  Table[
   ListEqualAreaPointPlot[clusteredpoints[[i]], 0.02, 
    Style -> {EdgeForm[Black], Opacity[0.5], colorlist[[i]]}], {i, 
    Length[clusteredpoints]}
   ]
  ]


ListEqualAreaContourPlot[data_, circlearea_, contourinterval_, 
  opts___] := 
 Module[{npts, A, \[Sigma], \[Mu], r, gridpts, gridvals, xycoords, ct,
    outline, contourplot, circlepts, circlemask, \[CapitalDelta], 
   nrows, ctable, x, y, c, maskrange, fontsize, interpolatedvalues, 
   gridsize},
  
  npts = Length[data];
  
  r = Sqrt[circlearea];
  A = \[Pi]*r^2;
  
  (*
  gridsize=0.1;
  *)
  
  gridsize = r;
  
  \[CapitalDelta] = 1;
  gridpts = 
   Table[{x + gridsize/2, 
     y + gridsize/2}, {x, -\[CapitalDelta], \[CapitalDelta], 
     gridsize}, {y, -\[CapitalDelta], \[CapitalDelta], gridsize}];
  nrows = Length[gridpts];
  xycoords = Map[EqualAreaXYCoords, data];
  gridvals = Table[0, {i, nrows}, {j, nrows}];
  
  Do[
   If[EucDist[xycoords[[c]], gridpts[[i, j]]] <= r, gridvals[[j, i]] ++],
   {i, nrows}, {j, nrows}, {c, npts}
   ];
  
  (*
  Do[
  Module[{dist},
  dist = EucDist[{0.,0.},gridpts[[i,j]]] ;
  If[dist > (1. - r )&& dist < (1. + r ),gridvals[[i,j]] *=  A/
  CircleOverlapArea[{0.,0.},gridpts[[i,j]],1.,r];
  ]
  ],
  {i,nrows},{j,nrows}
  ];
  *)
  
  Do[
   gridvals[[i, j]] = Min[{gridvals[[i, j]]/npts, 1.}], {i, 
    nrows}, {j, nrows}
   ];
  
  ctable = Table[c, {c, contourinterval, 1., contourinterval}];

  
  If[Max[gridvals] > 1, 
   Print["WARNING!  MAX GRIDVAL = ", Max[gridvals]]];
  
  
  (*
  interpolatedvalues=ListInterpolation[gridvals,{{-1,1},{-1,
  1}}];
  
  contourplot=ContourPlot[interpolatedvalues[y,x],{x,-1,1},{y,-1,1},
  Contours->ctable,Frame->False,DisplayFunction->Identity,Axes->False,
  ContourLabels->Automatic,PlotRange->{contourinterval,1.+
  contourinterval},opts];
  *)
  
  contourplot = 
   ListContourPlot[ gridvals, DataRange -> {{-1, 1}, {-1, 1}}, 
    Contours -> ctable, Frame -> False, DisplayFunction -> Identity, 
    Axes -> False, ContourLabels -> Automatic, MaxPlotPoints -> nrows,
     PlotRange -> {0.99*contourinterval, 1. + contourinterval}, 
    opts];
  
  
  fontsize = 0.025;
  
  outline = Graphics[{
     Black,
     Circle[{0., 0.}, 1.],
     Table[
      Line[{{0.98 Cos[\[Theta]], 0.98 Sin[\[Theta]]}, {Cos[\[Theta]], 
         Sin[\[Theta]]}}], {\[Theta], 0, 2 \[Pi], \[Pi]/18}],
     
     Text[
      Style["0", FontSize -> Scaled[fontsize]], {0, 1 + 2 fontsize}],
     Text[
      Style["90", FontSize -> Scaled[fontsize]], {1 + 2 fontsize, 
       0}],
     Text[
      Style["180", 
       FontSize -> Scaled[fontsize]], {0, -1 - 2 fontsize}],
     Text[
      Style["270", FontSize -> Scaled[fontsize]], {-1 - 2 fontsize, 0}]
     
     }
    ];
  
  maskrange = 1.05;
  
  circlepts = 
   Table[{Cos[\[Theta]], Sin[\[Theta]]}, {\[Theta], 0.01, 
     2. Pi - 0.01, Pi/100.}];
  
  circlepts = 
   AppendTo[
    circlepts, {{1.0, -10^-6}, {maskrange, -10^-6}, {maskrange, \
-maskrange}, {-maskrange, -maskrange}, {-maskrange, 
      maskrange}, {maskrange, maskrange}, {maskrange, 10^-6}, {1.0, 
      10^-6}}];
  
  circlepts = Partition[Flatten[circlepts], 2];
  
  circlemask = Graphics[{GrayLevel[1], Polygon[circlepts]}];
  
  Show[circlemask, contourplot, circlemask, outline, 
   AspectRatio -> 1., DisplayFunction -> $DisplayFunction, 
   Axes -> False]
  ]



ListKambPlot[data_, contourinterval_, opts___] := 
 Module[{npts, A, \[Sigma], \[Mu], r, gridpts, gridvals, xycoords, ct,
    outline, contourplot, circlepts, circlemask, \[CapitalDelta], 
   nrows, ctable, x, y, c, maskrange, fontsize},
  
  npts = Length[data];
  
  r = 3./Sqrt[(npts + 9)];
  
  A = \[Pi]*r^2;
  \[Mu] = npts*A;
  \[Sigma] = Sqrt[\[Mu] (1. - A)];
  \[CapitalDelta] = r*Floor[1/r];
  gridpts = 
   Table[{x, y}, {x, -\[CapitalDelta], \[CapitalDelta], 
     r}, {y, -\[CapitalDelta], \[CapitalDelta], r}];
  nrows = Length[gridpts];
  
  xycoords = Map[EqualAreaXYCoords, data];
  
  gridvals = Table[0, {i, nrows}, {j, nrows}];
  
  Do[
   If[EucDist[xycoords[[c]], gridpts[[i, j]]] <= r, gridvals[[j, i]] ++],
   {i, nrows}, {j, nrows}, {c, npts}
   ];
  
  (*
  
  Optional section to "correct" for the number of poles in counting \
circles that intersect the edge of the net. This may, however,
  over estimate densities along the edge and should be used with \
caution!
  
  Do[
  Module[{dist},
  dist = EucDist[{0.,0.},gridpts[[i,j]]] ;
  If[dist > (1. - r )&& dist < (1. + r ),gridvals[[i,j]] *=  A/
  CircleOverlapArea[{0.,0.},gridpts[[i,j]],1.,r]
  ]
  ],
  {i,nrows},{j,nrows}
  ];
  
  *)
  
  ctable = 
   Table[c, {c, contourinterval*\[Sigma], Max[gridvals]*\[Sigma], 
     contourinterval*\[Sigma]}];
  
  
  (*
  interpolatedvalues=ListInterpolation[gridvals,{{-1,1},{-1,
  1}}];
  contourplot=ContourPlot[interpolatedvalues[y,x],{x,-1,1},{y,-1,1},
  Contours->ctable,ContourShading->False,Frame->False,PlotRange->All,
  DisplayFunction->Identity,PlotPoints->50,Axes->False];
  *)
  
  contourplot = 
   ListContourPlot[ gridvals, DataRange -> {{-1, 1}, {-1, 1}}, 
    Contours -> ctable, Frame -> False, DisplayFunction -> Identity, 
    Axes -> False, ContourLabels -> Automatic, MaxPlotPoints -> nrows,
     PlotRange -> {Min[ctable], Max[ctable]}, opts];
  
  fontsize = 0.025;
  
  outline = Graphics[{
     Black,
     Circle[{0., 0.}, 1.],
     Table[
      Line[{{0.98 Cos[\[Theta]], 0.98 Sin[\[Theta]]}, {Cos[\[Theta]], 
         Sin[\[Theta]]}}], {\[Theta], 0, 2 \[Pi], \[Pi]/18}],
     
     Text[
      Style["0", FontSize -> Scaled[fontsize]], {0, 1 + 2 fontsize}],
     Text[
      Style["90", FontSize -> Scaled[fontsize]], {1 + 2 fontsize, 
       0}],
     Text[
      Style["180", 
       FontSize -> Scaled[fontsize]], {0, -1 - 2 fontsize}],
     Text[
      Style["270", FontSize -> Scaled[fontsize]], {-1 - 2 fontsize, 0}]
     
     }
    ];
  
  maskrange = 1.05;
  
  circlepts = 
   Table[{Cos[\[Theta]], Sin[\[Theta]]}, {\[Theta], 0.01, 
     2. Pi - 0.01, Pi/100.}];
  
  circlepts = 
   AppendTo[
    circlepts, {{1.0, -10^-6}, {maskrange, -10^-6}, {maskrange, \
-maskrange}, {-maskrange, -maskrange}, {-maskrange, 
      maskrange}, {maskrange, maskrange}, {maskrange, 10^-6}, {1.0, 
      10^-6}}];
  
  circlepts = Partition[Flatten[circlepts], 2];
  
  circlemask = Graphics[{GrayLevel[1], Polygon[circlepts]}];
  
  Print["N  = ", npts];
  Print["A = ", A/\[Pi], "%"];
  Print["\[Mu]  = ", \[Mu]];
  Print["\[Sigma]  = ", \[Sigma]];
  Print["CI = ", contourinterval, " \[Sigma]"];
  
  Show[circlemask, contourplot, circlemask, outline, 
   AspectRatio -> 1., DisplayFunction -> $DisplayFunction, 
   Axes -> False]
  ]

ListKambPlotOriginal[data_, contourinterval_, opts___] := 
 Module[{npts, A, \[Sigma], \[Mu], r, gridpts, gridvals, xycoords, ct,
    outline, contourplot, circlepts, circlemask, \[CapitalDelta], 
   nrows, ctable, x, y, c, maskrange, fontsize},
  
  npts = Length[data];
  
  r = 3./Sqrt[\[Pi] (npts + 9)];
  
  A = \[Pi]*r^2;
  \[Mu] = npts*A;
  \[Sigma] = Sqrt[\[Mu] (1. - A)];
  
  Print["N  = ", npts];
  Print["A = ", A];
  Print["\[Mu]  = ", \[Mu]];
  Print["\[Sigma]  = ", \[Sigma]];
  Print["CI = ", contourinterval, " \[Sigma]"];
  
  \[CapitalDelta] = r*Floor[1/r];
  gridpts = 
   Table[{x, y}, {x, -\[CapitalDelta], \[CapitalDelta], 
     r}, {y, -\[CapitalDelta], \[CapitalDelta], r}];
  nrows = Length[gridpts];
  
  xycoords = Map[EqualAreaXYCoords, data];
  
  gridvals = Table[0, {i, nrows}, {j, nrows}];
  
  Do[
   If[EucDist[xycoords[[c]], gridpts[[i, j]]] <= r, gridvals[[j, i]] ++],
   {i, nrows}, {j, nrows}, {c, npts}
   ];
  
  Do[
   Module[{dist},
    dist = EucDist[{0., 0.}, gridpts[[i, j]]] ;
    If[dist > (1. - r ) && dist < (1. + r ), 
     gridvals[[i, j]] *=  
      A/CircleOverlapArea[{0., 0.}, gridpts[[i, j]], 1., r]
     ]
    ],
   {i, nrows}, {j, nrows}
   ];
  
  ctable = 
   Table[c, {c, contourinterval*\[Sigma], Max[gridvals], 
     contourinterval*\[Sigma]}];
  
  (*
  interpolatedvalues=ListInterpolation[gridvals,{{-1,1},{-1,
  1}}];
  contourplot=ContourPlot[interpolatedvalues[y,x],{x,-1,1},{y,-1,1},
  Contours->ctable,ContourShading->False,Frame->False,PlotRange->All,
  DisplayFunction->Identity,PlotPoints->50,Axes->False];
  *)
  
  contourplot = 
   ListContourPlot[gridvals, DataRange -> {{-1, 1}, {-1, 1}}, 
    Contours -> ctable, Frame -> False, DisplayFunction -> Identity, 
    Axes -> False, ContourLabels -> Automatic, opts];
  
  fontsize = 0.025;
  
  outline = Graphics[{
     Black,
     Circle[{0., 0.}, 1.],
     Table[
      Line[{{0.98 Cos[\[Theta]], 0.98 Sin[\[Theta]]}, {Cos[\[Theta]], 
         Sin[\[Theta]]}}], {\[Theta], 0, 2 \[Pi], \[Pi]/18}],
     
     Text[
      Style["0", FontSize -> Scaled[fontsize]], {0, 1 + 2 fontsize}],
     Text[
      Style["90", FontSize -> Scaled[fontsize]], {1 + 2 fontsize, 
       0}],
     Text[
      Style["180", 
       FontSize -> Scaled[fontsize]], {0, -1 - 2 fontsize}],
     Text[
      Style["270", FontSize -> Scaled[fontsize]], {-1 - 2 fontsize, 0}]
     
     }
    ];
  
  maskrange = 1.05;
  
  circlepts = 
   Table[{Cos[\[Theta]], Sin[\[Theta]]}, {\[Theta], 0.01, 
     2. Pi - 0.01, Pi/100.}];
  
  circlepts = 
   AppendTo[
    circlepts, {{1.0, -10^-6}, {maskrange, -10^-6}, {maskrange, \
-maskrange}, {-maskrange, -maskrange}, {-maskrange, 
      maskrange}, {maskrange, maskrange}, {maskrange, 10^-6}, {1.0, 
      10^-6}}];
  
  circlepts = Partition[Flatten[circlepts], 2];
  
  circlemask = Graphics[{GrayLevel[1], Polygon[circlepts]}];
  
  Show[circlemask, contourplot, circlemask, outline, 
   AspectRatio -> 1., DisplayFunction -> $DisplayFunction, 
   Axes -> False]
  ]

EqualAreaPlaneArc[{\[Delta]_, \[Theta]_}] := Module[{\[Alpha]},
  
  Tooltip[
   Rotate[
    Line[
     Table[
      EqualAreaXYCoords[{ApparentDip[{\[Delta] , \[Alpha] }], 
        90 + \[Alpha]}], {\[Alpha], -90., 90., 2.}
      ]
     ], (90 - \[Theta] ) Degree, {0, 0}],
   {Round[\[Delta]] Degree, Round[\[Theta]] Degree}]
  ]


Options[ListStereoPlot] = {AspectRatio -> 1., Axes -> False, 
   PlotRange -> {{-1.2, 1.2}, {-1.2, 1.2}}};

ATan[x_, y_] := If[x  != 0 && y != 0, ArcTan[x, y], \[Pi]/2]

EqualAreaAngularCoords[{x_, y_}] := Module[{\[Theta], \[Delta], r},
  (*
  Note: Returns results in degrees.
  *)
  r = Sqrt[x^2 + y^2];
  \[Delta] = 1/2 (\[Pi] - 4 ArcSin[r/Sqrt[2.]]);
  \[Theta] = \[Pi]/2. - ATan[x, y];
  If[\[Theta] < 0, \[Theta] = \[Theta] + 2. \[Pi]];
  {\[Delta], \[Theta]}/Degree
  ]

EqualAreaXYCoords[{\[Delta]_, \[Theta]_}] := Module[{c},
  (* Expects input in degrees *)
  
  c = DiplineToCosines[{\[Delta], \[Theta]}];
  Chop[{c[[1]]/Sqrt[2] Sqrt[2./(1. - c[[3]])], 
    c[[2]]/Sqrt[2] Sqrt[2./(1. - c[[3]])]}]
  ]

EqualAreaMesh[] :=
 Module[{circles, spokes, \[Delta], \[Theta]},
  circles = 
   Table[Circle[{0, 0}, 
     Sqrt[2.] Cos[(\[Pi]/2. + \[Delta] Degree)/2.]], {\[Delta], 10, 
     80, 10}];
  spokes = 
   Table[Line[{{0, 0}, EqualAreaXYCoords[{0., \[Theta]}]}], {\[Theta],
      0, 350, 10}];
  Graphics[{GrayLevel[0.5], spokes, circles}]
  ]

CircleOverlapArea[pt1_, pt2_, r1_, r2_] :=
 
 Module[{angle1, angle2, c, area},
  c = EucDist[pt1, pt2];
  angle1 = 2. ArcCos[(r1^2 + c^2 - r2^2)/(2. r1 c)];
  angle2 = 2. ArcCos[(r2^2 + c^2 - r1^2)/(2. r2 c)];
  Return[0.5 (r1^2 (angle1 - Sin[angle1]) + 
      r2^2 (angle2 - Sin[angle2]))]
  ]


EucDist[pt1_List, pt2_List] :=
  
  N[Sqrt[(pt1[[1]] - pt2[[1]])^2 + (pt1[[2]] - pt2[[2]])^2]];

StereoXYCoords[{\[Delta]_, \[Theta]_}] := 
 
 (* Expects input in degrees *)
 
 Module[{r},
  r = Tan[\[Pi]/4. - (\[Delta] Degree )/2. ];
  Chop[{r*Sin[\[Theta] Degree ], r*Cos[\[Theta] Degree]}]
  ]

StereoLinePoint[{\[Delta]_, \[Theta]_}, pointrad_] := 
 Module[{az, r, dx, dy, graphicslist},
  Tooltip[
   Disk[StereoXYCoords[{\[Delta], \[Theta]}], 
    pointrad], {Round[\[Delta]/Degree] Degree, 
    Round[\[Theta]/Degree] Degree}]
  ]

StereoPlaneArc[{\[Delta]_, \[Theta]_}] := Module[{\[Alpha]},
  Tooltip[
   Rotate[
    Line[
     Table[
      StereoXYCoords[{ApparentDip[{\[Delta] , \[Alpha] }], 
        90 + \[Alpha]}], {\[Alpha], -90., 90., 2.}
      ]
     ], (90 - \[Theta] ) Degree, {0, 0}],
   {Round[\[Delta]] Degree, Round[\[Theta]] Degree}]
  ]

ListStereoPointPlot[indata_, pointrad_: 0.02, opts___] :=
 
 Module[{data, len, graphicslist, styles, fontsize},
  
  data = indata ;
  len = Length[data];
  
  Off[ReplaceAll::argr];
  styles = Style /. opts;
  On[ReplaceAll::argr];
  
  fontsize = 0.025;
  
  graphicslist = {};
  
  If[
   Length[styles] > 1, AppendTo[graphicslist, styles]
   ];
  
  AppendTo[graphicslist, {Table[
     StereoLinePoint[data[[i]] , pointrad], {i, len}
     ],
    Black, Opacity[1.],
    Circle[{0., 0.}, 1.],
    Table[
     Line[{{0.98 Cos[\[Theta]], 0.98 Sin[\[Theta]]}, {Cos[\[Theta]], 
        Sin[\[Theta]]}}], {\[Theta], 0, 2 \[Pi], \[Pi]/18}
     ],
    Text[Style["0", FontSize -> Scaled[fontsize]], {0, 
      1 + 2 fontsize}],
    Text[Style["90", FontSize -> Scaled[fontsize]], {1 + 2 fontsize, 
      0}],
    Text[Style["180", 
      FontSize -> Scaled[fontsize]], {0, -1 - 2 fontsize}],
    Text[Style["270", FontSize -> Scaled[fontsize]], {-1 - 2 fontsize,
       0}],
    Options[ListStereoPlot]}
   ];
  
  Graphics[Flatten[graphicslist]]
  ]

ListStereoArcPlot[indata_, opts___] :=
 
 Module[{data, len, graphicslist, styles, fontsize},
  
  data = indata;
  
  len = Length[data];
  
  Off[ReplaceAll::argr];
  styles = Style /. opts;
  On[ReplaceAll::argr];
  
  fontsize = 0.025;
  
  graphicslist = {};
  
  If[
   Length[styles] >= 1, AppendTo[graphicslist, styles]
   ];
  
  AppendTo[graphicslist, {Table[
     StereoPlaneArc[data[[i]]], {i, len}
     ],
    Thickness[Small],
    Black, Opacity[1.],
    Circle[{0., 0.}, 1.],
    Table[
     Line[{{0.98 Cos[\[Theta]], 0.98 Sin[\[Theta]]}, {Cos[\[Theta]], 
        Sin[\[Theta]]}}], {\[Theta], 0, 2 \[Pi], \[Pi]/18}
     ],
    Text[Style["0", FontSize -> Scaled[fontsize]], {0, 
      1 + 2 fontsize}],
    Text[Style["90", FontSize -> Scaled[fontsize]], {1 + 2 fontsize, 
      0}],
    Text[Style["180", 
      FontSize -> Scaled[fontsize]], {0, -1 - 2 fontsize}],
    Text[Style["270", FontSize -> Scaled[fontsize]], {-1 - 2 fontsize,
       0}],
    Options[ListStereoPlot]}
   ];
  
  Graphics[Flatten[graphicslist]]
  
  ]

FisherRandomDeviate[{dip_, ddn_}, \[Kappa]_] :=
 
 Module[{\[CapitalPhi], \[CapitalTheta], \[Lambda], R1, R2, vertical, 
   target, rotationmatrix},
  (*
   \[CapitalPhi] is dip direction and \[CapitalTheta] is dip \
angle somewhat
  following the convention of Fisher, Lewis, and
  Embelton (1987) 
  *)
  
  R1 = Random[];
  R2 = Random[];
  \[Lambda] = Exp[-2 \[Kappa]];
  
  (*
  First, calculate random deviates for a vertical dip vector
  *)
  \[CapitalTheta] = (\[Pi]/2. - 
      2 ArcSin[
        Sqrt[-Log[R1 (1 - \[Lambda]) + \[Lambda]]/(2 \[Kappa])]])/
    Degree;
  \[CapitalPhi] = (2. \[Pi] R2)/Degree;
  
  vertical = {0, 0, -1};
  target = DiplineToCosines[{dip, ddn}];
  rotationmatrix = RotationMatrix[{vertical, target}];
  
  CosinesToDipline[
   rotationmatrix.DiplineToCosines[{\[CapitalTheta], \[CapitalPhi]}]
   ]
  ]

EqualAreaFrictionCircle[\[Phi]_, opts___] := 
 Module[{styles, graphicslist, radius},
  Off[ReplaceAll::argr];
  styles = Style /. opts;
  On[ReplaceAll::argr];
  
  graphicslist = {};
  
  radius = Sqrt[2.] Sin[\[Pi]/4. - (\[Phi] Degree )/2. ];
  
  If[Length[styles] > 1, AppendTo[graphicslist, styles]];
  AppendTo[graphicslist, 
   Tooltip[Disk[{0, 0}, radius], "Friction Circle" ]];
  Graphics[Flatten[graphicslist]]
  
  ]

EqualAreaMarklandPlot[planes_List, slope_, \[Phi]_, opts___] :=
 
 (* Input expected in DIP and DIP DIRECTION! *)
 
 Module[{len, arcs1, arcs2, frictioncircle, frictionradius, x, y, 
   minx, maxx, r, circle1, circle2, criticalzone, interpolatedarc, 
   sltn},
  
  len = Length[planes];
  
  arcs1 = Table[
    EqualAreaPlaneArc[planes[[i]]], {i, len}
    ];
  
  arcs2 = EqualAreaPlaneArc[slope];
  
  frictioncircle = 
   EqualAreaFrictionCircle[\[Phi], 
    Style -> {EdgeForm[Black], Opacity[0.25], Yellow}];
  
  frictionradius = Sqrt[2] Sin[\[Pi]/4. - (\[Phi] Degree )/2. ];
  
  circle1 = x^2 + y^2 <= frictionradius^2;
  
  (*
  Find maxx by calculating apparent dip at the point where the \
friction circle intersects the slope arc
  *)
  
  (*
  Off[Solve::ifun];
  sltn=Solve[ApparentDip[{slope[[1]],\[Alpha]}]==\[Phi],\[Alpha]];
  maxx=Round[Abs[\[Alpha] /. sltn[[1]]]];
  On[Solve::ifun];
  *)
  
  maxx = N[
    Abs[ArcCos[Cot[slope[[1]] Degree] Tan[\[Phi] Degree]]]/Degree];
  
  circle2 =
   Table[
    EqualAreaXYCoords[{ApparentDip[{slope[[
         1]], \[Alpha] }], \[Alpha]}], {\[Alpha], -maxx, maxx, 0.5}
    ];
  
  
   circle1 =
   Table[EqualAreaXYCoords[{\[Phi], \[Theta]}], {\[Theta], -maxx, 
     maxx, 1}];
  
  criticalzone =
   {Red, Opacity[0.5], Polygon[Join[circle1, Reverse[circle2]]]};
  
  
  
  criticalzone =
   Graphics[
    Rotate[Tooltip[criticalzone, 
      "Unstable Zone"], -slope[[2]] Degree, {0, 0}]
    ];
  
  
  Show[EqualAreaMesh[], frictioncircle, criticalzone, 
   Graphics[{Blue, arcs1, Thick, Red, arcs2}], OuterCircle[]]
  
  
  ]


End[]

EndPackage[]
