(* ::Package:: *)

(* ::Title::Initialization:: *)
(*Diffusion coefficients in orbital space*)


(* ::Section::Initialization:: *)
(*Import data*)


(* ::Input::Initialization:: *)
tabaj=Import["../../data/Dump_Diffusion_Coefficients.hf5",{"Datasets","tabaj"}];
tabDRRjj=Import["../../data/Dump_Diffusion_Coefficients.hf5",{"Datasets","tabDRRjj"}]; 
tabDNRjj=Import["../../data/Dump_Diffusion_Coefficients.hf5",{"Datasets","tabDNRjj"}]; 
tabDjj=Import["../../data/Dump_Diffusion_Coefficients.hf5",{"Datasets","tabDjj"}]; 
tabaj={tabaj[[All,2]],tabaj[[All,1]]}//Transpose;
DRRjjTable=Partition[{tabaj,tabDRRjj}//Transpose//Flatten,3];
DNRjjTable=Partition[{tabaj,tabDNRjj}//Transpose//Flatten,3];
DjjTable=Partition[{tabaj,tabDjj}//Transpose//Flatten,3];


(* ::Section::Initialization:: *)
(*Plot data*)


(* ::Subsection::Initialization:: *)
(*SRR coefficients*)


(* ::Input::Initialization:: *)
pRRjj=ListContourPlot[DRRjjTable,Contours->30,PlotRange->{All,All,All},ScalingFunctions->{"Log10","Log10","Log10"},ColorFunction->"Rainbow",PlotLegends->Automatic,Frame->True,AspectRatio->1,ImageSize->Medium,FrameLabel->{Style["j",Medium,Bold],Style["a (mpc)",Medium,Bold]}];


(* ::Subsection::Initialization:: *)
(*NR diffusion coefficients*)


(* ::Input::Initialization:: *)
pNRjj=ListContourPlot[DNRjjTable,Contours->30,PlotRange->{All,All,All},ScalingFunctions->{"Log10","Log10","Log10"},ColorFunction->"Rainbow",PlotLegends->Automatic,Frame->True,AspectRatio->1,ImageSize->Medium,FrameLabel->{Style["j",Medium,Bold],Style["a (mpc)",Medium,Bold]}];


(* ::Subsection::Initialization:: *)
(*Total diffusion coefficients*)


(* ::Input::Initialization:: *)
pjj=ListContourPlot[DjjTable,Contours->30,PlotRange->{All,All,All},ScalingFunctions->{"Log10","Log10","Log10"},ColorFunction->"Rainbow",PlotLegends->Automatic,Frame->True,AspectRatio->1,ImageSize->Medium,FrameLabel->{Style["j",Medium,Bold],Style["a (mpc)",Medium,Bold]}];


(* ::Section::Initialization:: *)
(*Save data*)


(* ::Subsection::Initialization:: *)
(*SRR diffusion coefficients*)


(* ::Input::Initialization:: *)
Export["../../graphs/Mathematica/DRRjj.png",pRRjj];


(* ::Subsection::Initialization:: *)
(*NR diffusion coefficients *)


(* ::Input::Initialization:: *)
Export["../../graphs/Mathematica/DNRjj.png",pNRjj];


(* ::Subsection::Initialization:: *)
(*Total diffusion coefficients*)


(* ::Input::Initialization:: *)
Export["../../graphs/Mathematica/Djj.png",pjj];
