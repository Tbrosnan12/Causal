d = 3;
l = 1;
sample := 1000000;
plotLength := 8;

(*Import CSV filename*)
dataFileName = $CommandLine[[4]]
Print[dataFileName];

(*Extract n and m from filename*)
n = ToExpression@First@StringCases[dataFileName,"n=" ~~ x : DigitCharacter .. ~~ "_" :> x];
m = ToExpression@First@StringCases[dataFileName,"m=" ~~ x : DigitCharacter .. ~~ "_" :> x];
h = ToExpression@First@StringCases[dataFileName,"h=" ~~ x : (DigitCharacter .. ~~ ("." ~~ DigitCharacter ..) ...) :> x];


pointsTime = Import[dataFileName];

pdfFileName = 
  StringJoin["~/Documents/Causal/scripts/pdfs/plot_n=", 
   ToString[n], "_m=", ToString[m], "_h=", ToString[h],  "_Re[two_point].pdf"];
pdfName =  StringJoin["plot_n=", 
   ToString[n], "_m=", ToString[m], "_h=", ToString[h],  "_Re[two_point].pdf"];
    
Print[pdfName];
