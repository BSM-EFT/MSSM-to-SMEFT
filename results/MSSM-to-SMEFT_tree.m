(* ::Package:: *)

(* ::Title:: *)
(*MSSM-to-SMEFT matching at one-loop order*)


(* ::Subtitle:: *)
(*Works  on  Matchete  git commit "f8ed8d53"  (and later):*)
(*[ https://gitlab.com/matchete/matchete/-/commit/f8ed8d533e17003baf3006163c6e9c8f2acd8a11 ]*)


(* ::Chapter::Closed:: *)
(*Preparations*)


(* ::Text:: *)
(*Set current directory:*)


SetDirectory@NotebookDirectory[];


(* ::Text:: *)
(*Path  to  the MSSM  model  file:*)


$ModelFileMSSM= FileNameJoin[{NotebookDirectory[], "MSSM.m"}];


(* ::Text:: *)
(*Directory  where  Matchete  looks  for  model  files:*)


$ModelFileMatchete= FileNameJoin[{$MatchetePath, "Models", "MSSM.m"}];


(* ::Text:: *)
(*Copy the MSSM model file to the Matchete directory so that LoadModel finds it.*)
(*Possibly add the model file to the files ignored by git (add the line "Models/MSSM.m"  to the file .git/info/exclude in the Matchete repository).*)


CopyFile[
	$ModelFileMSSM,
	$ModelFileMatchete
	,
	OverwriteTarget->True
];


(* ::Chapter:: *)
(*The computation*)


(* ::Section:: *)
(*Loading the MSSM*)


(* ::Text:: *)
(*Loading the model:*)


(*\[ScriptCapitalL]MSSM= LoadModel[
	"MSSM",
	(* renaming of some default parameters names *)
	ModelParameters->{
		"gs"  -> g3 (*strong coupling*)
		,"gL" -> g2 (*weak coupling*)
		,"gY" -> g1 (*hypercharge coupling*)
		,"Yu" -> yu
		,"Yd" -> yd
		,"Ye" -> ye
	}
];*)


\[ScriptCapitalL]MSSM= LoadModel["MSSM"];


(* ::Subsubsection::Closed:: *)
(*IntroduceEffectiveCouplings *)


(*
\[ScriptCapitalL]MassFree = \[ScriptCapitalL]MSSM/.{Coupling[mq|mu|md|ml|me|m\[CapitalPhi]|\[Mu]H|mB|mW|mG,___]->0};
\[ScriptCapitalL]Mass     = \[ScriptCapitalL]MSSM-\[ScriptCapitalL]MassFree;
\[ScriptCapitalL]MSSM     = \[ScriptCapitalL]Mass+IntroduceEffectiveCouplings[\[ScriptCapitalL]MassFree, EffectiveCouplingSymbol->"\[Kappa]"];
*)


(* ::Subsection:: *)
(*Restricting to sUp, sCharm, sTop*)



\[ScriptCapitalL]MSSM = \[ScriptCapitalL]MSSM/.{Field[Alternatives@@{\[CapitalSigma],lt,et,qt,ut,dt,Gt,Wt,Bt},___]->0};



(*
\[ScriptCapitalL]MSSM = \[ScriptCapitalL]MSSM/.{Field[Alternatives@@{\[CapitalPhi],lt,et,dt,Gt,Wt},___]->0};
*)


(* ::Subsection:: *)
(*Print the MSSM Lagrangian:*)


Print["The MSSM Lagrangian is:"];
Print@Format[HcSimplify@\[ScriptCapitalL]MSSM, NiceForm]



(*
Print["--- where we have defined ---"];
PrintEffectiveCouplings[\[ScriptCapitalL]MSSM];
*)


(* ::Section::Closed:: *)
(*Tree-level matching*)


(* ::Text:: *)
(*We first perform only the tree-level matching. Calculating the tree-level result separately simplifies the proper treatment of evanescent operators, which only appear at tree level. (One-loop generated evanescent operators only affect two-loop matching calculations.) *)


\[ScriptCapitalL]EFT0= EOMSimplify[
	Match[\[ScriptCapitalL]MSSM, EFTOrder->6, LoopOrder->0]
	,
	ReductionIdentities->dDimensional
];


(* ::Section:: *)
(*Mapping to Warsaw basis*)


(* ::Text:: *)
(*In this section the MSSM matching results are mapped onto the Warsaw basis,*)


(* this is just a performance improvement *)
SetOptions[EOMSimplify,DummyCoefficients->True,ReductionIdentities->dDimensional];


(* ::Text:: *)
(*Load the Warsaw basis Lagrangian and Fierz it to the "Matchete basis". This only Fierzes the \!\(TraditionalForm\`*)
(*SubsuperscriptBox[*)
(*StyleBox["Q", "TI"], *)
(*StyleBox[*)
(*RowBox[{"l", "e", "q", "u"}], "TI"], *)
(*RowBox[{"(", "3", ")"}]]\) operator to a leptoquark-like basis.*)


\[ScriptCapitalL]Warsaw= EchoTiming[EOMSimplify[
	LoadModel["SMEFT", ModelParameters->{"gs"->g3,"gL"->g2,"gY"->g1}],
	ReductionIdentities->FourDimensional
]
,"GreensSimplify"];


(* ::Text:: *)
(*Remove all Baryon number violating terms, due to some issues with Fierz identities*)


\[ScriptCapitalL]Warsaw= \[ScriptCapitalL]Warsaw/.{Coupling[cllHH|cduu|cqqq|cqqu|cduq,___]->0};


(* ::Text:: *)
(*Determine the matching conditions in the Warsaw basis.*)


(* remove all poles *)
\[ScriptCapitalL]EFT$OnShell$simpl = CollectOperators[Matchete`PackageScope`BetterExpand[\[ScriptCapitalL]EFT$OnShell]/.\[Epsilon]^-1->0];


MatchingCondition= EchoTiming[
	MapEffectiveCouplings[
		\[ScriptCapitalL]EFT0,
		ReplaceEffectiveCouplings[\[ScriptCapitalL]Warsaw]
		,AppendEffectiveCouplingsDefs -> True
		,EOMSimplify                  -> False
		,ReductionIdentities          -> dDimensional
		,ShiftRenCouplings            -> True(*False*)
	]
,"Determining matching conditions"];


(* ::Text:: *)
(*Export the results*)


Export[FileNameJoin[{"results","MSSM-matching-conditions.m"}],MatchingCondition];
