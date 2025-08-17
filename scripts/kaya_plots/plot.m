d = 3;
l = 1;
sample := 100000;
plotLength := 8;

(*Import CSV filename*)
dataFileName = $CommandLine[[4]]

(*Extract n and m from filename*)
n = ToExpression@
   First@StringCases[dataFileName, 
     "n=" ~~ x : DigitCharacter .. ~~ "_" :> x];
m = ToExpression@
   First@StringCases[dataFileName, 
     "m=" ~~ x : DigitCharacter .. ~~ ".csv" :> x];

pointsTime = Import[dataFileName];


Subscript[m, c] := Sqrt[3]/2;
M := Sqrt[m^2 + Subscript[m, c]^2];
d := 3;
H1 = (d - 1)/2 + Sqrt[((d - 1)/2)^2 - M^2 l^2];
H2 = (d - 1)/2 - Sqrt[((d - 1)/2)^2 - M^2 l^2];
norm = Gamma[H1]*Gamma[H2]/((4*\[Pi])^(d/2)*l^2*Gamma[d/2]);
W[\[Tau]_] := 
  norm*Hypergeometric2F1[H1, H2, d/2, 1/2*(1 + Cosh[\[Tau]])];
Wpodal[\[Tau]_] := 
  norm*Hypergeometric2F1[H1, H2, d/2, 1/2*(1 - Cosh[\[Tau]])];
Walpha[\[Tau]_, \[Alpha]_, \[Beta]_] := 
  Cosh[2 \[Alpha]]*W[\[Tau]] + 
   Sinh[2 \[Alpha]]*(Cos[\[Beta]]*Re[Wpodal[\[Tau]]] - 
      Sin[\[Beta]]*Im[Wpodal[\[Tau]]]);


pdfFileName = 
  StringJoin["~/Documents/Causal/scripts/kaya_plots/pdfs/plot_n=", 
   ToString[n], "_m=", ToString[m], "_Re[two_point].pdf"];
pdfName =  StringJoin["plot_n=", 
   ToString[n], "_m=", ToString[m], "_Re[two_point].pdf"];


PlottingPoints = RandomSample[pointsTime, sample];
Export[pdfFileName,Show[ListPlot[PlottingPoints,PlotStyle->\
PointSize[0.00001],ImageSize->600],Plot[Re[W[\[Tau]]],{\[Tau],0,\
plotLength},PlotStyle->Red,PlotRange->{{0,plotLength},{-.2,.2}}],\
AxesOrigin->{0,0},PlotRange->{{0,plotLength},{-.2,.2}},AxesLabel->{\
"Proper time \[Tau]","Re[W]"}]];

Print[pdfName];
