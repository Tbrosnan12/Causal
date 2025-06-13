(* ::Package:: *)

(* ::Title:: *)
(*Causal Sets Package*)


(* ::Section::Closed:: *)
(*Definitions*)


BeginPackage["Causet`"];
(* Make CS Package symbols appear in green *)
SetOptions[$FrontEnd, 
  AutoStyleOptions -> {"SymbolContextStyles" ->
     {"Causet`" -> RGBColor[121/255, 151/255, 63/255]}
   }
]
SetOptions[Graphics, BaseStyle -> {FontFamily -> "Helvetica"}];
SetOptions[ListPlot, BaseStyle -> {FontFamily -> "Helvetica"}];


(* ::Item:: *)
(*Access:*)


CSize::usage   = "CSize[Sprinkling] returns the cardinality of Sprinkling";
CType::usage   = "CType[Sprinkling] returns the type of Sprinkling";
CDim::usage    = "CDim[Sprinkling] returns the spacetime dimension of Sprinkling";
CCoords::usage = "CCoords[sprinkling, i] returns the list of coordinate-i values of the points in Sprinkling";
CPoints::usage = "CPoints[Sprinkling] returns the list of sprinkled points";
CTime::usage   = "CTimes[Sprinkling] returns the list of time-coordinate values of the sprinkled points";
CSpace::usage  = "CSpace[Sprinkling] returns the list of spatial coordinate values of the sprinkled points";
CParams::usage = "CParams[Sprinkling] returns the list of the parameters that specify the geometry of the sprinkling";
CVolume::usage = "CVolume[Sprinkling] returns the spacetime volume of Sprinkling";


(* ::Item:: *)
(*Geometries:*)


CSMinkowski2Rectangle::usage      = "CSMinkowski2Rectangle[width_:1, height_:1, size] returns a Sprinkling into a rectangle in Minkowski2";
CSMinkowski3Rectangle::usage      = "CSMinkowski2Rectangle[width_:1, depth_: 1., height_:1, size] returns a Sprinkling into a rectangle in Minkowski3";
CSMinkowski2Diamond::usage        = "CSMinkowski2Diamond[height_:1, size] returns a Sprinkling into a causal diamond in Minkowski2";
CSMinkowski3Cylinder::usage       = "Returns a Sprinkling into a cylinder in Minkowski3";
CSMinkowski3Diamond::usage        = "Returns a Sprinkling into a causal diamond in Minkowski3";
CSConformalSpace2Rectangle::usage = "CSConformalSpace2Rectangle[function,...] returns a Sprinkling into a conformally flat spacetime in two dimensions whose volume measure (=function) depends only on the space coordinate.";
CSConformalTime2Rectangle::usage  = "CSConformalTime2Rectangle[function,...]  returns a Sprinkling into a conformally flat spacetime in two dimensions whose volume measure (=function) depends only on the time coordinate.";
CSDeSitterDGlobalSlab::usage      = "Returns a Sprinkling into a slab (d-1)-sphere x [0, height] of d-dimensional de Sitter space in closed global coordinates";
CSDeSitterDGlobalSlabFull::usage      = "Returns a Sprinkling into a slab (d-1)-sphere x [-height, height] of d-dimensional de Sitter space in closed global coordinates";
CSMinkowski4Diamond::usage        = "Returns a Sprinkling into a causal diamond in Minkowski4";
(*CSMinkowskiDDiamond::usage      = "CSMinkowskiDDiamond[D, Height:1., Size] returns a Sprinkling into a D-dimensional diamond."*)


(* ::Item:: *)
(*Order & Connectivity:*)


CausalMatrix::usage = "CausalMatrix[Sprinkling] returns the causal matrix of the sprinkling.";
LinkMatrix::usage   = "LinkMatrix[Sprinkling] returns the link matrix of the sprinkling. Dependencies: CausalMatrix[Sprinkling].";
Causal::usage = "Causal[Sprinkling, i, j] returns True/False according to whether the ith element is to the causal past of the jth element";
DistanceMatrix::usage = "DistanceMatrix[sprinkling] returns a matrix whose (ij)th entry gives the geodesic distance between atoms i and j.";
TestTau::usage = "TestTau[sprinkling] returns a matrix whose (ij)th entry gives the geodesic distance between atoms i and j.";
IntervalCardinalityMatrix::usage = "IntervalCardinalityMatrix[Sprinkling] returns a matrix whose (i, j\!\(\*SuperscriptBox[\()\), \(th\)]\) entry is the cardinality of the (exclusive) order interval [i, j]";
IntervalAbundances::usage = "IntervalAbundance[Sprinkling] returns a list of pairs {k, N[k]} where N[k] is the number of cardinality-k exclusive order intervals in the sprinkling. IntervalAbundance[sprinkling, k] returns the number of cardinality-k exclusive order intervals in the sprinkling.";
LongestChainLength::usage = "LongestChainLength[Sprinkling, i, j] returns the length of the longest chain between the ith and jth element in Sprinkling. Returns -\[Infinity] if i and j are not related.";


(* ::Item:: *)
(*Subsets*)


InRegion::usage  = "InRegion[Sprinkling, Function, Format] returns a list of elements in the sprinkling for whose sprinkling-coordinates Function is True. Format specifies the output format: ''Coordinates'' for coordinate values, ''Labels'' for labels. Default is ''Labels''.";
BelowSurface::usage = "BelowSurface[Sprinkling, Function] takes a Function T[x] of the spatial coordinates specifying a spacelike surface and returns a list of integers corresponding to the causal set elements below that surface.";
AboveSurface::usage = "BelowSurface[Sprinkling, Function] takes a Function T[x] of the spatial coordinates specifying a spacelike surface and returns a list of integers corresponding to the causal set elements above that surface.";
MinimalAboveSurface::usage = "MinimalAboveSurface[Sprinkling, Function] takes a Function t[x] of the spatial coordinates specifying a spacelike surface and returns a list of integers corresponding to the maximal causal set elements below that surface.";
MaximalBelowSurface::usage = "MaximalBelowSurface[Sprinkling, Function] takes a Function t[x] of the spatial coordinates specifying a spacelike surface and returns a list of integers corresponding to the maximal causal set elements below that surface.";
PastSet::usage = "PastSet[Sprinkling, element] counts the cardinality of the past set of element";
PastSets::usage = "PastSets[Sprinkling] counts the cardinality of the past sets of all elements in Sprinkling";
FutureSet::usage = "FutureSet[Sprinkling, element] counts the cardinality of the future set of the element";
FutureSets::usage = "FutureSets[Sprinkling] counts the cardinality of the future sets of all elements in Sprinkling";
MinimalElements::usage = "Minimal[Sprinkling] returns all the minimal elements in the sprinkling.";
MaximalElements::usage = "Maximal[Sprinkling] returns all the maximal elements in the sprinkling.";


(* ::Item:: *)
(*Visualisation:*)


CPlot::usage      = "CPlot[CSPrinkling] plots a sprinkled causal set in two or three dimensions.";
CEmbedding::usage = "CEmbedding[Sprinkling] gives the coordinates of some Euclidean embedding of the sprinkling";
LinkLines::usage  = "LinkLines[Sprinkling] gives a list of the lines corresponding to links for use in graphics";
(*Deprecated: CLabels::usage    = "CLabels[Sprinkling] gives a list of labels with coordinates for use in plots";*)
ShowLinks::usage  = "Option for CPlot. Set to True to display links";
ShowLabels::usage = "Option for CPlot. Set to True to display labels";
Highlight::usage  = "Option for CPlot. Highlight -> i highlights element i. Hightlight -> {\!\(\*SubscriptBox[\(i\), \(1\)]\), \!\(\*SubscriptBox[\(i\), \(2\)]\),...} highlights a list of elements.";
HighlightStyle::usage = "Option for CPlot to be used with Highlight. E.g. HighlightStyle -> {Blue, PointSize[10]}.";


(* ::Item:: *)
(*Conversion to Mathematica Graph:*)


CGraph::usage  = "CGraph[Sprinkling] makes a Mathematica Graph obtained from the sprinkling";


(* ::Item:: *)
(*Propagators*)


Retarded2D::usage = "Retarded2D[Sprinkling, m] returns the two-dimensional Johnston retarded propagator for a scalar field of mass m";
Wightman2D::usage = "Wightman2D[Sprinkling, m] returns the two-dimensional Sorkin-Johnston Wightman function for a scalar field of mass m";


(* ::Item:: *)
(*Nonlocal d'Alembertian*)


RetardedDAlembertian::usage = "RetardedDAlembertian[sprinkling_, lk_: 5] computes the retarded d'Alembertian for the sprinkling";
PlaneWave::usage = "PlaneWave[sprinkling_, amp_: 1, w_: 10, dir_: in, subset_List: All] discretises plane wave on a sprinkling";
NonlocalPropagator::usage = "NonlocalPropagator[sprinkling_] returns the retarded propagator of the retarded d'Alembetian";
NonlocalWightman::usage = "NonlocalWightman2D[sprinkling_] returns nonlocal retarded Green function and its SJ Wightman function";


(* ::Item:: *)
(*Action*)


Action::usage = "Returns action of spacetime region";


(* ::Item:: *)
(*Miscellaneous:*)


(*Private:
blue::usage         = "Custom color";
orange::usage       = "Custom color";
turquoise::usage    = "Custom color";
lightgreen::usage   = "Custom color";
Rows::usage         = "Rows[v] builds a square matrix from Length@v rows of v";
Cols::usage         = "Cols[v] builds a square matrix from Length@v columns of v";
ran::usage          = "ran[n] is shorthand for RandomReal[{n}]";
ranSpherical::usage = "ranSpherical[n, d] produces n uniformly distributed points on the d-sphere.";
ArcDistance::usage  = "ArcDistance[X,Y] gives the lenght of the minor arc segment between X and Y on the d-sphere given in embedding Euclidean coordinates.";*)


Begin["`Private`"];


(* ::Section:: *)
(*Geometries*)


(* ::Text:: *)
(*This is where the different sprinkling geometries are defined. A sprinkling is a list: {list of points, type, dimension, cardinality, volume, geometry parameters}.*)
(**)
(*For each geometry, you need to add the following:*)
(*1. A function CSGeometry that produces a list of points randomly sprinkled into that geometry*)
(*2. A method in CausalMatrix[] that produces the causal matrix for that geometry*)
(*3. A method in CEmbedding[] that produces a list of Cartesian coordinates that can be plotted in two or three dimensions.*)


(*Using Module is the fastest way I could find. Defining local static 
  variables w, h, s makes the List operations much faster somehow. 
  Should sprinkle 2^20 points in O(0.1s).*)
CSMinkowski2Rectangle[width_: 1. , height_: 1. , size_Integer] :=
Module[{w = width, h = height, s = size, coords},
coords = Transpose@{w * RandomReal[1, {s}], h * Sort@RandomReal[1, {s}]};
{coords,"CSMinkowski2Rectangle", 2, s, w * h, {w, h}}
]

CSMinkowski3Rectangle[width_: 1. , depth_: 1., height_: 1. , size_Integer] :=
Module[{w = width, h = height, d = depth, s = size, coords},
coords = Transpose@{w * RandomReal[1, {s}], d * RandomReal[1, {s}], h * Sort@RandomReal[1, {s}]};
{coords,"CSMinkowski3Rectangle", 3, s, w * d * h, {w, d, h}}
]


CSMinkowski2Diamond[height_: 1. , size_Integer] := 
Module[{h = height, s = size, coords, r},
 coords = h * Transpose@{RandomReal[1, {s}], Sort@RandomReal[1, {s}]};
 r = N@RotationMatrix[45Degree];
 coords = r . # & /@ coords;
 coords = SortBy[coords, Last];
{coords, "CSMinkowski2Diamond", 2, s, h * h, h}
]

CSMinkowski3Cylinder[radius_:1., height_:1., size_Integer] := 
Module[{r = radius, h = height, s = size, coords},
 coords = Transpose@{r * Sqrt@RandomReal[1, {s}], 2\[Pi] * RandomReal[1, {s}], h * Sort@RandomReal[1, {s}]};
{coords,"CSMinkowski3Cylinder", 3, s, \[Pi] * r^2 * h, {r, h}}
]

CSConformalTime2Rectangle[fn_, width_: 1., height_: 1., size_Integer] :=
Module[{w = width, h = height, s = size, norm, pdf, coords, t},
norm = NIntegrate[fn[t], {t, 0, h}];
pdf =  ProbabilityDistribution[fn[t] / norm, {t, 0., h}];
coords = Transpose@{RandomReal[{0, w},{s}], Sort@RandomVariate[pdf, s]};
{coords, "CSConformalTime2Rectangle", 2, s, w * h, {w, h}}
]

CSConformalSpace2Rectangle[fn_, width_: 1., height_: 1., size_Integer] :=
Module[{w = width, h = height, s = size, norm, pdf, coords, x},
norm   = NIntegrate[fn[x], {x, 0, w}];
pdf    = ProbabilityDistribution[fn[x] / norm, {x, 0., w}];
coords = Transpose@{RandomVariate[pdf, s], Sort@RandomReal[{0, w}, {s}]};
{coords, "CSConformalSpace2Rectangle", 2, s, h * norm, {w, h}}
]

(*
Method taken from S. V. Stefanescu "Generating uniform random points inside a cone"
This method could easily be generalised for the diamond in general d I think.

CSMinkowski3Diamond[height_: 1., size_Integer] :=
Module[{h = height, s = size, shalf, coords, Y, T, upcone, downcone},
Y:={h * T * Sqrt@RandomReal[1, {shalf}], 2\[Pi] * RandomReal[1, {shalf}]};
shalf = Floor[s/2];
T = Reverse@Sort@(RandomReal[1., {shalf}]^(1 / 3.));
upcone = Y~Join~{h * (1 - T)};
shalf = Ceiling[s/2];
T = Sort[RandomReal[1., {shalf}]^(1 / 3.)];
downcone = Y~Join~{h * (T - 1)};
coords = Transpose[Join[downcone, upcone, 2]];
{coords, "CSMinkowski3Diamond", 3, s, \[Pi] * h^3 / 6.0,  h}
]
*)

(* NOT AS EFFICIENT AS MIKE'S CODE, PROBABLY *)
CSMinkowski3Diamond[height_: 1., size_Integer] :=
Module[{h = height, s = size, p, coords},
coords={};
While[Length@coords<s,
	p= {RandomReal[{-h/2,h/2}],RandomReal[{-h/2,h/2}],RandomReal[{-h/2,h/2}]};
	If[(p[[3]]-h/2)^2>p[[1]]^2+p[[2]]^2&&(p[[3]]+h/2)^2>p[[1]]^2+p[[2]]^2,
	AppendTo[coords,p]
	]
];
coords=Sort[coords,#1[[3]]<#2[[3]]&];
{coords, "CSMinkowski3Diamond", 3, s, \[Pi] * h^3/12.0, h}
]

CSMinkowski4Diamond[height_: 1., size_Integer] :=
Module[{h = height, s = size, p, coords},
coords={};
While[Length@coords<s,
	p= {RandomReal[{-h/2,h/2}],RandomReal[{-h/2,h/2}],RandomReal[{-h/2,h/2}],RandomReal[{-h/2,h/2}]};
	If[(p[[4]]-h/2)^2>p[[1]]^2+p[[2]]^2+p[[3]]^2&&(p[[4]]+h/2)^2>p[[1]]^2+p[[2]]^2+p[[3]]^2,
	AppendTo[coords,p]
	]
];
coords=Sort[coords,#1[[4]]<#2[[4]]&];
{coords, "CSMinkowski4Diamond", 4, s, \[Pi] * h^4/24.0, h}
]

CSMinkowskiDDiamond[dimension_:2, height_: 1., size_Integer] :=
Module[{h = height, s = size, d = dimension, coords},
Y := {h * T * Sqrt@RandomReal[1, {shalf}], 2\[Pi] * RandomReal[1, {shalf}]};
shalf = Floor[s/2];
T = Reverse@Sort@(RandomReal[1., {shalf}]^(1 / 3.));
upcone = Y~Join~{h * (1 - T)};
shalf = Ceiling[s/2];
T= Sort[RandomReal[1., {shalf}]^(1 / 3.)];
downcone = Y~Join~{h * (T - 1)};
coords = Transpose[Join[downcone, upcone, 2]];
{coords, "CSMinkowski3Diamond", 3, s, \[Pi] * h^3 / 6.0,  h}
]

(*Spacetime dimension d*)
CSDeSitterDGlobalSlab[height_: 1., dimension_Integer, size_Integer] :=
Module[{h = height, s = size, d = dimension, pdf, norm, \[Eta], coords},
norm = Hypergeometric2F1[0.5, 0.5 * (1 + d), 1.5, Sin[h]^2] * Sin[h];
pdf = ProbabilityDistribution[Sec[\[Eta]]^d / norm, {\[Eta], 0, h}];
coords = Transpose[Transpose[Causet`Private`ranSpherical[s, d - 1]]~Join~{Sort@RandomVariate[pdf, {s}]}];
{coords, "CSDeSitterDGlobalSlab", d, s, norm * d \[Pi]^(d/2) Gamma[1 + d / 2]^-1,  h}
]

CSDeSitterDGlobalSlabFull[height_: 1., dimension_Integer, size_Integer] :=
Module[{h = height, s = size, d = dimension, pdf, norm, \[Eta], coords},
norm = 2 Hypergeometric2F1[0.5, 0.5 * (1 + d), 1.5, Sin[h]^2] * Sin[h];
pdf = ProbabilityDistribution[Sec[\[Eta]]^d / norm, {\[Eta], -h , h}];
coords = Transpose[Transpose[Causet`Private`ranSpherical[s, d - 1]]~Join~{Sort@RandomVariate[pdf, {s}]}];
{coords, "CSDeSitterDGlobalSlabFull", d, s, norm * d \[Pi]^(d/2) Gamma[1 + d / 2]^-1,  h}
]


(*"de Sitter global slab sprinkling: 
norm is the spacetime volume element, dV=1/cos^d(t) 2 \[Pi]^(d/2)/Gamma[d/2] integrated from t=-h or t=0 to t=h,
the constant factors of d, \[Pi]^d/2 and Gamma are being omitted and tagged on at the end.
Since the spacetime volume doesn't depend on the spatial spherical part, a sphere is sprinkled with constant 
probability and 
then also a line for the time is sprinkled with 1/cos^d(t) probability distribution and the two sets of numbers 
are then paired together."*)



(* ::Section::Closed:: *)
(*Visualisation*)


(* ::Text:: *)
(*The CPlot function plots a visualisation of the causal set in two or three dimensions. CPlot automatically chooses the dimension (2 or 3) of the plot for any given sprinkling, which for any new geometry must be specified in the CEmbedding method below.*)


CEmbedding[sprinkling_List] := Module[{type = CType@sprinkling, embedding},
embedding = 
Which[
    MemberQ[{"CSMinkowski2Rectangle", "CSMinkowski2Diamond", "CSConformalSpace2Rectangle", "CSConformalTime2Rectangle"}, type],
    "Cartesian2",
    MemberQ[{"CSMinkowski3Cylinder", "CSMinkowski3Diamond"}, type],
    "Polar3",
    MemberQ[{"CSMinkowski3Rectangle"}, type],
    "Cartesian3"
];

Switch[embedding,
	"Cartesian2", CPoints[sprinkling],
	"Polar3",    {#[[1]] * Cos[#[[2]]], #[[1]] * Sin[#[[2]]], #[[3]]}& /@ CPoints[sprinkling],
	"Cartesian3", CPoints[sprinkling]
]
]


LinkLines[sprinkling_List] := 
Module[{C = CEmbedding@sprinkling, L = LinkMatrix@sprinkling},
Pick[Flatten@Outer[Line[{#1,#2}]&, C, C, 1], Flatten@L, 1]
]

(*Deprecated: CLabels[sprinkling_List] := Text[#, CEmbedding[sprinkling][[#]]]& /@ Range[CSize@sprinkling]*)


(* Here are default options for CPlot. They can be overwritten in CPlot[]. The pattern is described in http://bit.ly/1hPag7u *)
Options[CPlot] = {Frame ->  True, ShowLinks -> False, ShowLabels -> False, Highlight -> None, HighlightStyle -> {}}~Join~Options@ListPlot~Join~Options@ListPointPlot3D;

(* Use SyntaxInformation to make sure that custom options in CPlot are correctly highlighted *)
SyntaxInformation[CPlot] = {"ArgumentsPattern" -> {__, OptionsPattern[]}};

(* The parameter \[Epsilon] determines the fraction by which the plot range exceeds the range of sprinkled points. *)
CPlot[sprinkling_List, opts:OptionsPattern[{ListPlot, ListPointPlot3D, CPlot}]] :=
Module[{plot, points = CEmbedding@sprinkling, size = CSize@sprinkling, dim = CDim@sprinkling, \[Epsilon] = 0.05, max, min, rangerule, stylerule, pointsize, display, allrules2D, allrules3D},
(*If[!OptionValue[Highlight]===None, points = Delete[points, {#} & /@ OptionValue@Highlight]]; (* Remove any highlighted points from the list of main points*)
*)
pointsize = Which[size >= 300 && dim == 3, 0.01, size >= 500 && dim == 2, 0.0075, True, 0.0125];
max = Max /@ Transpose@points;
min = Min /@ Transpose@points;
(*TODO: rangerule might be deprecated. need to check that. *)
rangerule = PlotRange -> Transpose@{min - 1 * pointsize , max + 1 * pointsize};
stylerule = PlotStyle -> Directive[If[OptionValue[ShowLinks], Black, blue], PointSize[pointsize]];
allrules2D = {opts, stylerule, PlotRange -> All, Axes -> None, Frame -> True, AspectRatio -> #[[2]]/#[[1]]&@(max-min)};
allrules3D = {opts, stylerule, rangerule, BoxRatios -> max - min};
Which[
dim == 2, 
    plot = ListPlot[points, Evaluate@FilterRules[allrules2D, Options@ListPlot]];
	display = {plot};
	If[OptionValue[ShowLabels],        AppendTo[display, Graphics[Text[#[[1]], #[[2]] - 0.015 * (max-min)]& /@ Transpose@{Range@CSize@sprinkling, CEmbedding@sprinkling}]]];
    If[OptionValue[ShowLinks],        PrependTo[display, Graphics[{blue, Sequence@@LinkLines@sprinkling}]]];
    If[!OptionValue[Highlight]===None, AppendTo[display, Graphics[{Sequence@@Flatten@{PointSize[2/3 * pointsize], White, OptionValue@HighlightStyle} , Point /@ points[[{Sequence@@OptionValue@Highlight}]]}]]];
	Show[display, Evaluate@FilterRules[allrules2D, Options@Graphics]],
dim == 3, 
    plot = ListPointPlot3D[points, stylerule, Evaluate@FilterRules[allrules3D, Options@ListPlot]];
	display = {plot};
	If[OptionValue[ShowLabels],        AppendTo[display, Graphics3D[Text[#[[1]], #[[2]] - 0.015 * (max-min)]& /@ Transpose@{Range@CSize@sprinkling, CEmbedding@sprinkling}]]];  
	If[OptionValue[ShowLinks],        PrependTo[display, Graphics3D[{blue, Sequence@@LinkLines@sprinkling}]]];
    If[!OptionValue[Highlight]===None, AppendTo[display, Graphics3D[{Sequence@@Flatten@{PointSize[2/3 * pointsize], White, OptionValue@HighlightStyle}, Point /@ points[[{Sequence@@OptionValue@Highlight}]]}]]];
	Show[display, Evaluate@FilterRules[allrules3D, Options@Graphics3D]]
]
];


(* ::Section:: *)
(*Order & Connectivity*)


Causal[s_List, i_Integer, j_Integer] := (CausalMatrix[s, {i, j}][[1, 2]] == 1)

CausalMatrix[sprinkling_List] := CausalMatrix[sprinkling, All];

CausalMatrix[sprinkling_List, subset_] := 
Which[
    (* Conformally flat geometries with Cartesian coordiantes in d = 2*)
    MemberQ[{"CSMinkowski2Rectangle", "CSMinkowski2Diamond", "CSConformalSpace2Rectangle", "CSConformalTime2Rectangle"}, CType@sprinkling],
    Module[{t = CTime[sprinkling][[subset]], x = CSpace[sprinkling][[subset]]}, 
      UnitStep[Abs@Outer[Plus, -x, x]~Subtract~Outer[Plus, -t, t]]~BitXor~1
    ],

    (* 4D Minkowski with Cartesian coordiantes
    MemberQ[{"CSMinkowski4Diamond"}, CType@sprinkling],
    Module[{t = CTime[sprinkling][[subset]], x = CSpace[sprinkling][[subset]](*x = CCoords[sprinkling,1][[subset]], y = CCoords[sprinkling, 2][[subset]], z = CCoords[sprinkling,3][[subset]]*)}, 
    UnitStep[Map[Norm, Outer[Plus, x, -x, 1], {2}]~Subtract~Outer[Plus, -t, t]]~BitXor~1],*)

    (* Flat space with polar coordinates in d = 3 *)    
    MemberQ[{"CSMinkowski3Cylinder"(*, "CSMinkowski3Diamond"*)}, CType@sprinkling],
    Module[{r = CCoords[sprinkling, 1][[subset]], \[Theta] = CCoords[sprinkling, 2][[subset]], t = CTime[sprinkling][[subset]]}, 
    UnitStep[(Causet`Private`Rows[r]^2 + Causet`Private`Cols[r]^2 - 2.0 * Outer[Times, r, r] * Cos[Outer[Plus, -\[Theta], \[Theta]]])~Subtract~Outer[Plus, -t, t]]~BitXor~1],
    
    (* De Sitter space in global coordinates with conformal time coordinate*)
    MemberQ[{"CSDeSitterDGlobalSlab","CSDeSitterDGlobalSlabFull"}, CType@sprinkling],
    Module[{t = CTime[sprinkling][[subset]], x = CSpace[sprinkling][[subset]]}, 
    UnitStep[(Re@ArcCos@Outer[Dot, x, x, 1])~Subtract~Outer[Plus, -t, t]]~BitXor~1],
    
    (* Flat space in general dimensions in Cartesian coordinates *)
    MemberQ[{"CSMinkowskiDRectangle", "CSMinkowski3Rectangle", "CSMinkowski3Diamond", "CSMinkowski4Diamond"}, CType@sprinkling],
    Module[{t = CTime[sprinkling][[subset]], x = CSpace[sprinkling][[subset]]}, 
    UnitStep[Map[Norm, Outer[Plus, x, -x, 1], {2}]~Subtract~Outer[Plus, -t, t]]~BitXor~1]
    ]
(*For de Sitter: The arccos part gets the spatial distance or arc on the sphere since the radius has unit length, then this is 
compared to the temporal separation to determined the causal relation between the pair of points (timelike if the temporal difference
is greater than the spatial difference and spacelike if vice versa. UnitStep makes a 1 for anything \[GreaterEqual] 0 (spacelike or null) and 0 otherwise. 
~BitXor~1 then changes any 1's to 0's and any 0's to 1's (in the unlikely case of null separated pairs, we would get a wrong entry).)*)
(* Method: there is a link between i and j iff (i precedes j) and (there are no elements in [i,j])*)
LinkMatrix[sprinkling_List, subset_: All] :=
Module[{c = CausalMatrix[sprinkling, subset]},
c~BitAnd~UnitStep[-c . c]
]

DistanceMatrixC[sprinkling_List, subset_: All, restriction_String: "None"] :=
Module[{geodist, type = CType@sprinkling},
geodist = 
Which[
    (* Conformally flat geometries with Cartesian coordiantes in d = 2*)
    MemberQ[{"CSMinkowski2Rectangle", "CSMinkowski3Rectangle", "CSMinkowski2Diamond", "CSConformalSpace2Rectangle", "CSConformalTime2Rectangle"}, type],
    Module[{t = CTime[sprinkling][[subset]], x = CSpace[sprinkling][[subset]]}, 
    Sqrt[Abs@Outer[Plus, -x, x]~Subtract~Abs@Outer[Plus, -t, t]]],
    
	(* Flat space with polar coordinates in d = 3 *)    
    MemberQ[{"CSMinkowski3Cylinder", "CSMinkowski3Diamond"}, type],
    Module[{r = CCoords[sprinkling, 1][[subset]], \[Theta] = CCoords[sprinkling, 2][[subset]], t = CTime[sprinkling][[subset]]}, 
    Sqrt[(Causet`Private`Rows[r]^2 + Causet`Private`Cols[r]^2 - 2.0 * Outer[Times, r, r] * Cos[Outer[Plus, -\[Theta], \[Theta]]])~Subtract~Abs@Outer[Plus, -t, t]]],
    (* De Sitter space in global coordinates with conformal time coordinate*)
    
    MemberQ[{"CSDeSitterDGlobalSlab","CSDeSitterDGlobalSlabFull"}, CType@sprinkling],
    Module[{t = CTime[sprinkling][[subset]], x = CSpace[sprinkling][[subset]]}, 
    UnitStep[(Re@ArcCos@Outer[Dot, x, x, 1])~Subtract~Outer[Plus, -t, t]]],
    
    MemberQ[{"CSMinkowski3Rectangle", "CSMinkowskiDRectangle"}, CType@sprinkling],
    Module[{t = CTime[sprinkling][[subset]], x = CSpace[sprinkling][[subset]]}, 
    UnitStep[Map[Norm, Outer[Plus, x, -x, 1], {2}]~Subtract~Outer[Plus, -t, t]]]
    ];
Which[restriction == "None", Abs@geodist, restriction == "Timelike", Im@geodist, restriction == "Spacelike", Re@geodist]
]

TestTau[sprinkling_List, subset_: All, restriction_String: "None"] :=
Module[{geodist, type = CType@sprinkling},
geodist = 
Which[    
    MemberQ[{"CSDeSitterDGlobalSlab","CSDeSitterDGlobalSlabFull"}, CType@sprinkling],
    Module[{t = CTime[sprinkling][[subset]], x = CSpace[sprinkling][[subset]]}, 
   (ArcCos@Outer[Dot, x, x, 1])~Subtract~Outer[Plus, -t, t]]
    ];
Which[restriction == "None", Abs@geodist, restriction == "Timelike", Im@geodist, restriction == "Spacelike", Re@geodist]
]

IntervalCardinalityMatrix[sprinkling_List] := Dot[#, #]&@(IdentityMatrix[CSize@sprinkling] + CausalMatrix@sprinkling);

IntervalAbundances[sprinkling_List] := Module[{c = CausalMatrix@sprinkling, nonrelations, out},
out = {#, Count[Flatten@IntervalCardinalityMatrix@sprinkling, #]} & /@ Range[0, Length@sprinkling - 2];
(*The formula above doesn't work for links (0-intervals). The following lines remedy this: *)
nonrelations = Count[Flatten@c, 0];
out[[1, 2]] = out[[1, 2]] - nonrelations;
out
];

IntervalAbundances[causalmatrix_List, n_Integer] := Count[Flatten@IntervalCardinalityMatrix@causalmatrix, n];



(* ::Section::Closed:: *)
(*Subsets*)


GetLabels[sprinkling_List, subset_]:= Range[CSize@sprinkling][[subset]];

InRegion[sprinkling_List, fn_, format_:"Labels"] := Switch[format, "Labels", Pick[Range@CSize@sprinkling, fn /@ CPoints@sprinkling], "Coordinates", Select[CPoints@sprinkling, fn]];

BelowSurface[sprinkling_List, fn_] := Select[Range@CSize@sprinkling, CTime[sprinkling][[#]] < fn[Sequence@@CSpace[sprinkling][[#]]] & ];

AboveSurface[sprinkling_List, fn_] := Select[Range@CSize@sprinkling, CTime[sprinkling][[#]] > fn[Sequence@@CSpace[sprinkling][[#]]] & ];

MinimalElements[input_List, option_:""] := Module[{c, r},
If[option == "FromCausalMatrix", 
c = input; r = Range@Length@input, 
c = CausalMatrix@input; r = Range@CSize@input];
Pick[r, Map[Not, MemberQ[#, 1] &/@ Transpose@c]]
]

MaximalElements[input_List, option_:""] := Module[{c , r},
If[option == "FromCausalMatrix", 
c = input; r = Range@Length@input, 
c = CausalMatrix@input; r = Range@CSize@input];
Pick[r, Map[Not, MemberQ[#, 1] &/@ c]]
]

MaximalBelowSurface[sprinkling_List, fn_] := 
Module[{b = BelowSurface[sprinkling, fn], c},
If[b == {}, Return[{}]; Break[]];
c = CausalMatrix[sprinkling, b];
Pick[b, 0 == # & /@ Total@Transpose@c]
]

MinimalAboveSurface[sprinkling_List, fn_] := 
Module[{b = AboveSurface[sprinkling, fn], c},
If[b == {}, Return[{}]; Break[]];
c = CausalMatrix[sprinkling, b];
Pick[b, 0 == # & /@ Total@c]
]

PastSet[sprinkling_List, element_Integer] := Pick[Range[element], CausalMatrix[sprinkling, Range[element]][[All, -1]], 1]

PastSets[input_List, option_:""] := Module[{s,c},
If[option == "FromCausalMatrix", 
c = input; s = Length@input, 
c = CausalMatrix@input; s = CSize@input];
Pick[Table[Range[s],{s}], Transpose@c,1]
]

FutureSet[sprinkling_List, element_Integer] := Pick[Range[element, CSize@sprinkling], CausalMatrix[sprinkling, Range[element, CSize@sprinkling]][[1, All]], 1]

FutureSets[input_List, option_:""] := Module[{s, c},
If[option == "FromCausalMatrix", 
c = input; s = Length@input, 
c = CausalMatrix@input; s = CSize@input];
Pick[Table[Range[s],{s}], c, 1]
]

(*PastSet[sprinkling_List, element_Integer] := Total[CausalMatrix[sprinkling, Range[element]][[All, element]]]

PastSets[sprinkling_List] := Total@CausalMatrix@sprinkling

FutureSet[sprinkling_List, element_Integer] := Total[CausalMatrix[sprinkling, Range[element, CSize@sprinkling]][[element, All]]]

FutureSets[sprinkling_List] := Total@Transpose@CausalMatrix[sprinkling]*)


(* ::Section::Closed:: *)
(*Propagators*)


Wightman2D[sprinkling_, mass_Real] :=
Module[{m = mass, \[Rho] = CSize@sprinkling / CVolume@sprinkling, C=SparseArray@CausalMatrix@sprinkling, s = CSize@sprinkling, KR, i\[CapitalDelta], EE, VV, Q},
KR = 0.5 * N@C . Inverse[SparseArray@N@IdentityMatrix@s + 0.5 * m^2 \[Rho]^-1 C];
i\[CapitalDelta] = N@I * (KR - Transpose@KR);
{EE,VV} = Chop@Eigensystem@i\[CapitalDelta];
Q = Transpose[VV] . DiagonalMatrix[Clip[N[EE],{0, Max[N[EE]]}]] . Conjugate[VV] // Chop
]

Retarded2D[sprinkling_, mass_Real] :=
Module[{m = mass, \[Rho] = CSize@sprinkling / CVolume@sprinkling, C=SparseArray@CausalMatrix@sprinkling, s = CSize@sprinkling, KR, i\[CapitalDelta], EE, VV, Q},
KR = 0.5 * N@C . Inverse[SparseArray@N@IdentityMatrix@s + 0.5 * m^2 \[Rho]^-1 C]
]


(* ::Section::Closed:: *)
(*Nonlocal d'Alembertian*)


RetardedDAlembertian[sprinkling_, lk_: 5] :=
Module[{size = CSize@sprinkling, c = N@CausalMatrix[sprinkling], icm, dim = CDim@sprinkling, f, a, l, \[Epsilon]},
(*l = (CVolume@sprinkling / size)^(1./dim);*)
l = 1.;
\[Epsilon] = (l/lk)^dim;
Which[
\[NonBreakingSpace]\[NonBreakingSpace]dim == 2, 
	f = 4 l^-2 \[Epsilon]^2 (1-\[Epsilon])^(#-2) (1-2 \[Epsilon] (#-2)/(1-\[Epsilon])+\[Epsilon]^2 (#-2) (#-3)/(2 (1-\[Epsilon])^2))&; 
	a = -2 \[Epsilon]/l^2 ;,
  dim == 3,
	a = -(\[Pi]/(3Sqrt[2]))^(2/3)/(l^2 N[Gamma[5/3]]);
	f = -a \[Epsilon]^(5/3) (1-\[Epsilon])^(#-2) (1-27 \[Epsilon] (#-2)/(8(1-\[Epsilon]))+9\[Epsilon]^2 (#-2) (#-3)/(8 (1-\[Epsilon])^2))&;,
\[NonBreakingSpace]\[NonBreakingSpace]dim == 4, 
	f = l^-2 4. / Sqrt[6.] \[Epsilon]^1.5 (1 - \[Epsilon])^(# - 2) (1 - 9\[Epsilon] (#-2)/(1-\[Epsilon]) +8\[Epsilon]^2(#-2)(#-3)/(1-\[Epsilon])^2-4\[Epsilon]^3(#-2)(#-3)(#-4)/(3(1-\[Epsilon])^3))&; 
	a = - 4. \[Epsilon]^(1/2) /(l^2 Sqrt[6.]);
];

icm = SparseArray@Dot[#, #]&@(IdentityMatrix[size] + c);
{a IdentityMatrix[size] + c * f[icm], \[Epsilon]}(* Map[f, icm, {2}] *)
]

PlaneWave[sprinkling_List, amp_: 1, w_: 10., dir_String: "in", subset_List: All] :=
Module[{size = CSize@sprinkling, dim = CDim@sprinkling, t = CTime[sprinkling][[subset]], x = CSpace[sprinkling][[subset]], X, Y, Z, r},
Which[
	dim == 2,
	Which[
	dir == "in", amp Exp[-I w(t-x)],
	dir == "out", amp Exp[-I w(t+x)]
	],
	dim > 2,
	X = CCoords[sprinkling,1][[subset]]; 
	(*Y = CCoords[sprinkling,2][[subset]]; 
	Z = CCoords[sprinkling,3][[subset]];*)
	r = Map[Norm,CSpace[sprinkling]];
	Which[
	dir == "in", amp Exp[-I w(t-X)](*amp Exp[-I w(t-r)]/r*),
	dir == "out", amp Exp[-I w(t+r)]/r
	]
]
]

NonlocalPropagator[sprinkling_, lk_: 5] :=
Module[{Bk=SparseArray@Transpose@RetardedDAlembertian[sprinkling, lk][[1]], GR},
GR = N@Inverse[Bk]
]

NonlocalWightman[sprinkling_, lk_: 5] :=
Module[{GR = NonlocalPropagator[sprinkling, lk], s = CSize@sprinkling, i\[CapitalDelta], EE, VV, Q},
i\[CapitalDelta] = N@I*(GR - Transpose@GR);
{EE,VV} = Chop@Eigensystem@i\[CapitalDelta];
{GR, Q = Transpose[VV] . DiagonalMatrix[Clip[N[EE],{0, Max[N[EE]]}]] . Conjugate[VV]// Chop}
]


(* ::Section::Closed:: *)
(*Action*)


Action[sprinkling_, lk_: 5] := 
Module[{Bk=SparseArray@Transpose@RetardedDAlembertian[sprinkling, lk][[1]], size = CSize[sprinkling], R, S},
R = Dot[Bk,ConstantArray[-2,size]];
{R, S = Total[R]}
]


(* ::Section:: *)
(*Access To Sprinklings*)


CType[sprinkling_List]   := sprinkling[[2]]

CDim[sprinkling_List]    := sprinkling[[3]]

CSize[sprinkling_List]   := sprinkling[[4]]

CVolume[sprinkling_List] := sprinkling[[5]]

CParams[sprinkling_List] := sprinkling[[6]]

CTime[sprinkling_List]   := sprinkling[[1]][[All,-1]]

CCoords[sprinkling_List, n_Integer] := sprinkling[[1]][[All, n]]

CCoords[sprinkling_List] := Transpose@sprinkling[[1]]

CPoints[sprinkling_List] := sprinkling[[1]]

CSpace[sprinkling_List]  := If[CType[sprinkling] == "CSDeSitterDGlobalSlab"||CType[sprinkling] == "CSDeSitterDGlobalSlabFull", sprinkling[[1]][[All, 1 ;; -2]],
	
	If[CDim[sprinkling] == 2, Flatten@sprinkling[[1]][[All, 1 ;; -2]], sprinkling[[1]][[All, 1 ;; -2]]]
]


(* ::Section::Closed:: *)
(*Link to Mathematica Graphs*)


SyntaxInformation[CGraph] = {"ArgumentsPattern" -> {__, OptionsPattern[]}};
Options[CGraph] = {SprinklingCoordinates -> False, TransitiveReduction -> False}~Join~Options@Graph;
CGraph[sprinkling_List, opts:OptionsPattern[]] := 
Module[{coords = CEmbedding@sprinkling, m = If[OptionValue[TransitiveReduction], LinkMatrix[sprinkling], CausalMatrix[sprinkling]]},
    If[OptionValue[SprinklingCoordinates], 
       AdjacencyGraph[m, Evaluate@FilterRules[{opts}, Options@Graph], VertexCoordinates -> CEmbedding@sprinkling],
       AdjacencyGraph[m, Evaluate@FilterRules[{opts}, Options@Graph]]
    ]
]


(* ::Section:: *)
(*Miscellaneous and Private Functions*)


(* Colours *)
blue       = RGBColor[66/255,106/255,181/255];
orange     = RGBColor[238/255,133/255,37/255];
turquoise  = RGBColor[26/255,178/255,204/255];
lightgreen = RGBColor[121/255,151/255,63/255];

ran[n_Integer] := RandomReal[1., {n}]

ran[] := RandomReal[]

(* Using David's trick X/Norm[X] where X is a (d + 1)-vector of standard normals *)
ranSpherical[n_Integer, d_Integer]:= # / Norm[#] &/@ RandomVariate[NormalDistribution[], {n, d + 1}];

ranCircular[n_Integer]:=2 * \[Pi] * RandomReal[1.0,{n}];

Rows[a_List] :=  Transpose@Table[a, {Length@a}];

Cols[a_List] :=  Table[a, {Length@a}];

ArcDistance[x_List, y_List] := ArcCos[x . y];

End[];

EndPackage[];

