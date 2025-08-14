SetDirectory["~/Documents/Causal/"];(*Directory should be set to where the causet2.wl package is saved*)
<<"causet2.wl";(*importing the package causet2.wl which has things like the sprinkling and causal/link matrices defined in it*)
m:=0;
h:=1.3;
n:=10;
d:=3;
l:=1;
"causal events" n
"~ causal relations" n^2
A=CSDeSitterDGlobalSlabFull[h,d,n];
\[Rho]:=CSize[A] / CVolume[A];
"density" \[Rho]
Subscript[C, 1]=SparseArray[CausalMatrix[A]];



 L=LinkMatrix[A];
maxlength=SparseArray[Subscript[C, 1]]["NonzeroValues"];
positions=SparseArray[Subscript[C, 1]//Flatten]["NonzeroPositions"]//Flatten;
zeroCompare=Table[0,{i,1,n},{j,1,n}];
counter=1;
chains=L;

Do[
counter++;

chains=chains.L;(*finding chains of length i+1 composed of links*)

If[chains==zeroCompare,Break[]];(*if chains = 0, we can stop the loop, since no longer chains will be found*)

length=Flatten[chains][[positions]];(*flattening the matrix to a 1d array and finding only the entries corresponding to the positions of nonzero causal matrix entries*)

comparePositions=SparseArray[length//Flatten]["NonzeroPositions"]//Flatten;(*finding the positions in the new array that are nonzero*)

length=SparseArray[length]["NonzeroValues"];(*keeping only the nonzero entries of the new array*)

length=counter*HeavisideTheta[length];(*currently the length array has the number of link chains between two points, but we want the length of the chain*)
(*since the chains are longer in each iteration of the loop, we can replace the corresponding entries in the max length array with the updated values *)
maxlength[[comparePositions]]=length
(*looping up to NN-1 (i.e. L^(NN-1)), since the longest possible chain in a causet of size NN has NN-1 links*)
,{i,1,n-2}];

flattenedMaxLength=Table[0,{i,1,n*n}];(*NN^2 array of zeros*)

Do[flattenedMaxLength[[ positions[[i]] ]]=maxlength[[i]],{i,1,Length[positions]}];(*replace entries with max chain lengths according to the position of nonzero causal matrix entries in flattened causal matrix*)

MaxLength=Table[flattenedMaxLength[[1+(i-1)*n;;n+(i-1)*n]],{i,1,n}];(*unflatten to a square NN x NN matrix*)
m3=1.8701503305151406;
a=((\[Rho] \[Pi])/12)^(1/3) m3/(2\[Pi]);
K = Table[If[Subscript[C, 1][[i,j]]==1,(a*Cos[m/(2 \[Pi] a)*MaxLength[[i,j]]])/MaxLength[[i,j]],0],{i,n},{j,n}];
(*K=K0.Inverse[IdentityMatrix[n]+m^2/\[Rho]K0];*)

MfdTime =Table[Re[ArcCosh[1/(Cos[A[[1,i,4]]]*Cos[A[[1,j,4]]])*(A[[1,i,1]]*A[[1,j,1]]+A[[1,i,2]]*A[[1,j,2]]+A[[1,i,3]]*A[[1,j,3]]-Sin[A[[1,i,4]]]*Sin[A[[1,j,4]]])]*Subscript[C, 1][[i,j]]],{i,n},{j,n}];

(*MfdSpace=Table[Im[ArcCosh[1/(Cos[A[[1,i,4]]]*Cos[A[[1,j,4]]])*(A[[1,i,1]]*A[[1,j,1]]+A[[1,i,2]]*A[[1,j,2]]+A[[1,i,3]]*A[[1,j,3]]-Sin[A[[1,i,4]]]*Sin[A[[1,j,4]]])]],{i,n},{j,n}];
MfdSpace= MfdSpace*(UpperTriangularize[ConstantArray[1,{n,n}],1]-Subscript[C, 1]);*)

(*time=Transpose[{Flatten[MfdSpace],Flatten[N[K]]}];*)

(*time=Transpose[{Flatten[MfdTime],Flatten[K]}];
time=DeleteCases[Chop[time],{0,_}];*)

(*timeVol=Transpose[{Flatten[MfdTime],Flatten[K0]}];
timevol=DeleteCases[Chop[timeVol],{0,_}];*)







i\[CapitalDelta]=N[I*(K-Transpose[K])];
{EE,VV}=Chop[Eigensystem[i\[CapitalDelta]]];
Q = Transpose[VV] . DiagonalMatrix[Clip[N[EE],{0, Max[N[Abs[EE]]]}]] . Conjugate[VV] // Chop;

pointsTime=Transpose[{Flatten[MfdTime],Flatten[Re[Q]]}];
pointsTime=DeleteCases[Chop[pointsTime],{0,_}];
Length[pointsTime]"No. of non-zero relations"
(*pointsSpace=Transpose[{Flatten[MfdSpace],Flatten[Re[Q]]}];
*)


Subscript[m, c]:=Sqrt[3]/2;
M:= Sqrt[m^2+Subscript[m, c]^2];
d:=3;
H1=(d-1)/2+Sqrt[((d-1)/2)^2-M^2 l^2];
H2=(d-1)/2-Sqrt[((d-1)/2)^2-M^2 l^2];
norm=Gamma[H1]*Gamma[H2]/((4*\[Pi])^(d/2)*l^2*Gamma[d/2]);
W[\[Tau]_]:=norm*Hypergeometric2F1[H1,H2,d/2,1/2*(1+Cosh[\[Tau]])];
Wpodal[\[Tau]_]:=norm*Hypergeometric2F1[H1,H2,d/2,1/2*(1-Cosh[\[Tau]])];
Walpha[\[Tau]_,\[Alpha]_,\[Beta]_]:= Cosh[2\[Alpha]]*W[\[Tau]]+Sinh[2\[Alpha]]*(Cos[\[Beta]]*Re[Wpodal[\[Tau]]]-Sin[\[Beta]]*Im[Wpodal[\[Tau]]]);
(*
Walpha[1,0,1]
Plot[Re[Walpha[\[Tau],1,0]],{\[Tau],0,8},PlotStyle->Orange]


Plot[Walpha[\[Tau],1,1],{\[Tau],0,8},PlotStyle->Orange]
1/2ArcTanh[Abs[Sin[\[Pi]*Sqrt[((d-1)/2)^2-M^2 l^2]]]]
(d/2+HeavisideTheta[-Sin[\[Pi]*Sqrt[((d-1)/2)^2-M^2 l^2]]])\[Pi]

N[Re[norm*Hypergeometric2F1[H1,H2,d/2,1/2*(1+Cosh[1])]]]
N[Limit[Re[norm*Hypergeometric2F1[H1,H2,d/2,1/2*(1+Cosh[\[Epsilon]])]],{\[Epsilon]->0}]]
*)






sample:=2;
plotLength := 4;
PlottingPoints=RandomSample[pointsTime,sample];
Export["/home/thomas/Documents/Causal/scripts/plot.png", Show[ListPlot[PlottingPoints,PlotStyle->PointSize[0.00001],ImageSize->600],Plot[Re[W[\[Tau]]],{\[Tau],0,plotLength},PlotStyle->Red,PlotRange->{{0,plotLength},{-.2,.2}}],AxesOrigin->{0,0},PlotRange->{{0,plotLength},{-.2,.2}},AxesLabel->{"Proper time \[Tau]","Re[W]"}]]
(*,Plot[Re[Walpha[\[Tau],alpha,beta]],{\[Tau],0,8},PlotStyle->Orange,PlotRange->{{0,8},{-.1,.1}}]*)
