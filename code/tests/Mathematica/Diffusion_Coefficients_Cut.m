(* ::Package:: *)

(* ::Title::Initialization:: *)
(*(*(*(*(*(*(*(*(*Diffusion coefficients at fixed sma*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
col=ColorData[10]/@Range[20]
(*SetDirectory[NotebookDirectory[]]*)


(* ::Section::Initialization:: *)
(*(*(*(*(*(*(*(*(*Import data*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
tabj=Import["../../data/Dump_Diffusion_Coefficients_Cut.hf5",{"Datasets","tabj"}];
tabDRRjj=Import["../../data/Dump_Diffusion_Coefficients_Cut.hf5",{"Datasets","tabDRRjj"}]; 

tabDNRjj=Import["../../data/Dump_Diffusion_Coefficients_Cut.hf5",{"Datasets","tabDNRjj"}]; 
DRRjjTable={tabj,tabDRRjj}//Transpose;

DNRjjTable={tabj,tabDNRjj}//Transpose;



(* ::Section::Initialization:: *)
(*(*(*(*(*(*(*(*(*Plot data*)*)*)*)*)*)*)*)*)


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*(*SRR diffusion coefficients*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
pRRjj=ListPlot[DRRjjTable,PlotRange->{{0.05,1.0},{10^-5,1}},ScalingFunctions->{"Log10","Log10"},Joined->True,InterpolationOrder->1,PlotLegends->Automatic,Frame->True,AspectRatio->1,ImageSize->Medium,FrameLabel->{Style["j",Medium,Bold],Style["Djj (1/Myr)",Medium,Bold]}];


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*(*NR diffusion coefficients*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
pNRjj=ListPlot[DNRjjTable,PlotRange->{{0.05,1},{10^-5,1}},ScalingFunctions->{"Log10","Log10"},Joined->True,InterpolationOrder->1,PlotLegends->Automatic,Frame->True,AspectRatio->1,ImageSize->Medium,FrameLabel->{Style["j",Medium,Bold],Style["Djj (1/Myr)",Medium,Bold]}];


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*(*RR and NR diffusion coefficients together*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
pjj=ListPlot[{DRRjjTable,DNRjjTable},PlotRange->{{0.05,1},{10^-5,1}},ScalingFunctions->{"Log10","Log10"},Joined->True,InterpolationOrder->1,PlotLegends->{"\!\(\*SubscriptBox[SuperscriptBox[\(D\), \(RR\)], \(jj\)]\)","\!\(\*SubscriptBox[SuperscriptBox[\(D\), \(NR\)], \(jj\)]\)"},Frame->True,AspectRatio->1,ImageSize->Medium,FrameLabel->{Style["j",Medium,Bold],Style["Djj (1/Myr)",Medium,Bold]}];


(* ::Section::Initialization:: *)
(*(*(*(*(*(*(*(*(*Save data*)*)*)*)*)*)*)*)*)


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*(*Both RR and NR diffusion coefficients together*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
Export["../../graphs/Mathematica/DjjCut.png",pjj]



