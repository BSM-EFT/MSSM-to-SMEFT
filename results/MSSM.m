(* ::Package:: *)

(* ::Title:: *)
(*MSSM model file*)


(* ::Subtitle:: *)
(*Works  on  Matchete  git commit "f8ed8d53"  (and later):*)
(*[ https://gitlab.com/matchete/matchete/-/commit/f8ed8d533e17003baf3006163c6e9c8f2acd8a11 ]*)


(* ::Section:: *)
(*Lagrangian*)


MatcheteLagrangianParameters["MSSM"] = 
{
	"SU3c", "g3", "G", 
	"SU2L", "g2", "W", 
	"U1Y", "g1", "B", 
	"Flavor", 
	"q", "u", "d", "l", "e", "H",
	"yd", "yu", "ye", (*"\[Lambda]", "\[Mu]2",*)
	"Nf"-> 3,
	"\[Phi]", "mH2", "\[CapitalPhi]", "m\[CapitalPhi]",
	(*"Hu", "MHu", "Hd", "MHd", "Hdc",*)
	"\[CapitalSigma]",
	(*"Hut", "Hdt",*)
	"qt", "mq", "ut", "mu", "dt", "md", "lt", "ml", "et", "me",
	"Gt", "mG", "Wt", "mW", "Bt", "mB",
	"\[Mu]t", "au", "ad", "ae", "b",
	"s\[Beta]", "s2\[Beta]", "s4\[Beta]", "c\[Beta]", "c2\[Beta]", "c4\[Beta]"
};


MatcheteLagrangianAlphabets["MSSM"] = 
{
	"SU3cFundAlphabet" -> {"a","b","c","d","e","f"}, "SU3cAdjAlphabet" -> {"A","B","C","D","E","F"},
	"SU2LFundAlphabet" -> {"i","j","k","l","m","n"}, "SU2LAdjAlphabet" -> {"I","J","K","L","M","N"},
	"FlavorAlphabet"   -> {"p","r","s","t","u","v"}
};


MatcheteLagrangian["MSSM",ModParam_,IndAlphabet_]:= MatcheteLagrangian["MSSM",ModParam,IndAlphabet]=
Module[{\[ScriptCapitalL]free, \[ScriptCapitalL]gauge, \[ScriptCapitalL]SuperPotential, \[ScriptCapitalL]SoftBreaking, \[ScriptCapitalL]MSSM, relabel=Join[ModParam,IndAlphabet,{"Nf"->3}]},
	(* loading the SM definitions, the predefined SM Lagrangian is not used due to the 2HDM structure of the MSSM *)
	(*LoadModel["SM",ModelParameters->ModParam[[1;;22]],IndexAlphabet->IndAlphabet];*)
	ReleaseHold[MSSM`DefineSM[]/.relabel];

	(* define superpartner fields *)
	ReleaseHold[MSSM`DefineMSSMFields[]/.relabel];

	(* define MSSM couplings *)
	ReleaseHold[MSSM`DefineMSSMCouplings[]/.relabel];

	(* derive free MSSM Lagrangian *)
	\[ScriptCapitalL]free= ReleaseHold[MSSM`FreeLagrangianMSSM[]/.relabel];

	(* derive gauge Lagrangian *)
	\[ScriptCapitalL]gauge= ReleaseHold[MSSM`GaugeLagrangianMSSM[]/.relabel];

	(* derive gauge Lagrangian *)
	\[ScriptCapitalL]SuperPotential= ReleaseHold[MSSM`SuperPotentialMSSM[]/.relabel];

	(* soft SUSY breaking terms *)
	\[ScriptCapitalL]SoftBreaking= ReleaseHold[MSSM`SoftSUSYbreakingMSSM[]/.relabel];

	(* full MSSM Lagrangian *)
	\[ScriptCapitalL]MSSM= \[ScriptCapitalL]free + \[ScriptCapitalL]gauge + \[ScriptCapitalL]SuperPotential + \[ScriptCapitalL]SoftBreaking;

	(* rotate to Higgs basis, or the doublet mass basis before EWSB *)
	(*\[ScriptCapitalL]MSSM= ReleaseHold[MSSM`ToHiggsBasis[\[ScriptCapitalL]MSSM]/.relabel];*)
	\[ScriptCapitalL]MSSM= ReleaseHold[MSSM`ToHiggsDoubletMassBasis[\[ScriptCapitalL]MSSM]/.relabel];

	(* rotate to Higgsino mass basis *)
	\[ScriptCapitalL]MSSM= ReleaseHold[MSSM`ToHiggsinoMassBasis[\[ScriptCapitalL]MSSM]/.relabel];

	(*
	(* simplify trigonometric functions *)
	\[ScriptCapitalL]MSSM= TrigReduce/@ \[ScriptCapitalL]MSSM;
	*)
	
	(* rotate to Higgsino mass basis *)
	\[ScriptCapitalL]MSSM= ReleaseHold[MSSM`ReplaceRotationAngles[\[ScriptCapitalL]MSSM]/.relabel];

	(* simplify the results *)
	\[ScriptCapitalL]MSSM= GreensSimplify[\[ScriptCapitalL]MSSM,ReductionIdentities->dDimensional];

	(* simplify some terms using trigonometric idnetities *)
	\[ScriptCapitalL]MSSM= ReleaseHold[MSSM`TrigSimplifiy[\[ScriptCapitalL]MSSM]];

	(* remove intermediate auxiliary fields and couplings and properly set the field properties*)
	ReleaseHold[MSSM`RemoveTMPdef[]/.relabel]; (* without this line there are additional power-type supertraces computed, which are unnecessary... *)

	(* delete identities for operators containing superpartner fields *)
	Matchete`PackageScope`ResetOperatorAssociations[];

	Return[\[ScriptCapitalL]MSSM]
]


(* ::Chapter:: *)
(*Start MSSM context*)


Begin["MSSM`"];


(* ::Chapter:: *)
(*SM definitions*)


(* ::Section::Closed:: *)
(*Minimal definition of SM required for MSSM*)


DefineSM[]:=Hold@Module[{},
	DefineGaugeGroup["SU3c", SU@3, "g3", "G", FundAlphabet -> "SU3cFundAlphabet", AdjAlphabet -> "SU3cAdjAlphabet", NiceForm->{"\!\(\*SubscriptBox[\(g\), \(s\)]\)",Default}];
	DefineGaugeGroup["SU2L", SU@2, "g2", "W", FundAlphabet -> "SU2LFundAlphabet", AdjAlphabet -> "SU2LAdjAlphabet", NiceForm->{"\!\(\*SubscriptBox[\(g\), \(L\)]\)",Default}];
	DefineGaugeGroup["U1Y", U1, "g1", "B", NiceForm->{"\!\(\*SubscriptBox[\(g\), \(Y\)]\)",Default}];
	
	DefineFlavorIndex["Flavor","Nf",IndexAlphabet-> "FlavorAlphabet"];
	
	DefineField["q",Fermion,Indices->{"SU3c"[fund],"SU2L"[fund],"Flavor"},Charges->{"U1Y"[1/6]},Chiral-> LeftHanded,Mass->0];
	DefineField["u",Fermion,Indices->{"SU3c"[fund],"Flavor"},Charges->{"U1Y"[2/3]},Chiral-> RightHanded,Mass->0];
	DefineField["d",Fermion,Indices->{"SU3c"[fund],"Flavor"},Charges->{"U1Y"[-1/3]},Chiral-> RightHanded,Mass->0];
	DefineField["l",Fermion,Indices->{"SU2L"[fund],"Flavor"},Charges->{"U1Y"[-1/2]},Chiral-> LeftHanded,Mass->0];
	DefineField["e",Fermion,Indices->{"Flavor"},Charges->{"U1Y"[-1]},Chiral-> RightHanded,Mass->0];
	
	(* Light Higgs boson is defiend later in totation to mass basis *)

	DefineCoupling["yu",Indices->{"Flavor","Flavor"},NiceForm->"\!\(\*SubscriptBox[\(y\), \(u\)]\)"];
	DefineCoupling["yd",Indices->{"Flavor","Flavor"},NiceForm->"\!\(\*SubscriptBox[\(y\), \(d\)]\)"];
	DefineCoupling["ye",Indices->{"Flavor","Flavor"},NiceForm->"\!\(\*SubscriptBox[\(y\), \(e\)]\)"];
]


(* ::Chapter:: *)
(*MSSM definitions*)


(* ::Section:: *)
(*Definitions of fields & couplings*)


(* ::Subsection::Closed:: *)
(*Defining the fields*)


DefineMSSMFields[]:=Hold[
	(* 2HDM *)
	DefineField[
		ToExpression["Hu"],
		Scalar,
		Indices->{ToExpression["SU2L"][fund]},
		Charges->{ToExpression["U1Y"][1/2]},
		Mass->{Heavy,ToExpression["MHu"]},
		NiceForm->"\!\(\*SubscriptBox[\(H\), \(u\)]\)"
	];
	DefineField[
		ToExpression["Hd"],
		Scalar,
		Indices->{ToExpression["SU2L"][fund]},
		Charges->{ToExpression["U1Y"][-1/2]},
		Mass->{Heavy,ToExpression["MHd"]},
		NiceForm->"\!\(\*SubscriptBox[\(H\), \(d\)]\)"
	];
	
	(* sFermions *)
	DefineField[
		ToExpression["qt"],
		Scalar,
		Indices->{ToExpression["SU3c"][fund],ToExpression["SU2L"][fund],ToExpression["Flavor"]},
		Charges->{ToExpression["U1Y"][1/6]},
		Mass->{Heavy,ToExpression["mq"],{ToExpression["Flavor"]}},
		NiceForm->{"\!\(\*OverscriptBox[\(q\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), \(q\)]\)"}
	];
	DefineField[
		ToExpression["ut"],
		Scalar,
		Indices->{ToExpression["SU3c"][fund],ToExpression["Flavor"]},
		Charges->{ToExpression["U1Y"][2/3]},
		Mass->{Heavy,ToExpression["mu"],{ToExpression["Flavor"]}},
		NiceForm->{"\!\(\*OverscriptBox[\(u\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), \(u\)]\)"}
	];
	DefineField[
		ToExpression["dt"],
		Scalar,
		Indices->{ToExpression["SU3c"][fund],ToExpression["Flavor"]},
		Charges->{ToExpression["U1Y"][-1/3]},
		Mass->{Heavy,ToExpression["md"],{ToExpression["Flavor"]}},
		NiceForm->{"\!\(\*OverscriptBox[\(d\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), \(d\)]\)"}
	];
	DefineField[
		ToExpression["lt"],
		Scalar,
		Indices->{ToExpression["SU2L"][fund],ToExpression["Flavor"]},
		Charges->{ToExpression["U1Y"][-1/2]},
		Mass->{Heavy,ToExpression["ml"],{ToExpression["Flavor"]}},
		NiceForm->{"\!\(\*OverscriptBox[\(l\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), \(l\)]\)"}
	];
	DefineField[
		ToExpression["et"],
		Scalar,
		Indices->{ToExpression["Flavor"]},
		Charges->{ToExpression["U1Y"][-1]},
		Mass->{Heavy,ToExpression["me"],{ToExpression["Flavor"]}},
		NiceForm->{"\!\(\*OverscriptBox[\(e\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), \(e\)]\)"}
	];
	
	(* gauginos *)
	DefineField[
		ToExpression["Gt"],
		Fermion,
		Indices->{ToExpression["SU3c"][adj]},
		SelfConjugate->True,
		Mass->{Heavy,ToExpression["mG"]},
		NiceForm->{"\!\(\*OverscriptBox[\(G\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), \(G\)]\)"}
	];
	DefineField[
		ToExpression["Wt"],
		Fermion,
		Indices->{ToExpression["SU2L"][adj]},
		SelfConjugate->True,
		Mass->{Heavy,ToExpression["mW"]},
		NiceForm->{"\!\(\*OverscriptBox[\(W\), \(~\)]\)","\!\(\*SubscriptBox[\(M\), \(W\)]\)"}
	];
	DefineField[
		ToExpression["Bt"],
		Fermion,
		SelfConjugate->True,
		Mass->{Heavy,ToExpression["mB"]},
		NiceForm->{"\!\(\*OverscriptBox[\(B\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), \(B\)]\)"}
	];
	
	(* Higgsinos *)
	DefineField[
		ToExpression["Hut"],
		Fermion,
		Indices->{ToExpression["SU2L"][fund]},
		Charges->{ToExpression["U1Y"][1/2]},
		Chiral->LeftHanded,
		Mass->0, (* mass defined separately later *)
		NiceForm->"\!\(\*SubscriptBox[OverscriptBox[\(H\), \(~\)], \(u\)]\)"
	];
	DefineField[
		ToExpression["Hdt"],
		Fermion,
		Indices->{ToExpression["SU2L"][fund]},
		Charges->{ToExpression["U1Y"][-1/2]},
		Chiral->LeftHanded,
		Mass->0, (* mass defined separately later *)
		NiceForm->"\!\(\*SubscriptBox[OverscriptBox[\(H\), \(~\)], \(d\)]\)"
	];
]


(* ::Subsection::Closed:: *)
(*Defining the couplings*)


DefineMSSMCouplings[]:=Hold[
	DefineCoupling[
		ToExpression["\[Mu]t"], 
		EFTOrder->0, 
		SelfConjugate->True,
		NiceForm->"\!\(\*OverscriptBox[\(\[Mu]\), \(~\)]\)"
	];
	DefineCoupling[
		ToExpression["au"], 
		Indices->{ToExpression["Flavor"],ToExpression["Flavor"]}, 
		EFTOrder->0, 
		NiceForm->"\!\(\*SubscriptBox[\(a\), \(u\)]\)"
	];
	DefineCoupling[
		ToExpression["ad"], 
		Indices->{ToExpression["Flavor"],ToExpression["Flavor"]}, 
		EFTOrder->0, 
		NiceForm->"\!\(\*SubscriptBox[\(a\), \(d\)]\)"
	];
	DefineCoupling[
		ToExpression["ae"], 
		Indices->{ToExpression["Flavor"],ToExpression["Flavor"]}, 
		EFTOrder->0, 
		NiceForm->"\!\(\*SubscriptBox[\(a\), \(e\)]\)"
	];
	DefineCoupling[
		ToExpression["b"],
		EFTOrder->0,
		SelfConjugate->True
	];
]


(* ::Section:: *)
(*Writing the MSSM Lagrangian*)


(* ::Subsection::Closed:: *)
(*Free Lagrangian*)


FreeLagrangianMSSM[]:=Hold[
	FreeLag[
		ToExpression["G"]   , ToExpression["W"]  , ToExpression["B"]  ,
		ToExpression["Gt"]  , ToExpression["Wt"] , ToExpression["Bt"] ,
		ToExpression["q"]   , ToExpression["u"]  , ToExpression["d"]  , ToExpression["l"]  , ToExpression["e"]  ,
		ToExpression["qt"]  , ToExpression["ut"] , ToExpression["dt"] , ToExpression["lt"] , ToExpression["et"] ,
		ToExpression["Hu"]  , ToExpression["Hd"] ,
		ToExpression["Hut"] , ToExpression["Hdt"]
	]
]


(* ::Subsection::Closed:: *)
(*Gauge interactions (not encoded in covariant derivatives)*)


GaugeLagrangianMSSM[]:=Hold@Module[
	{
		(*auxiliary terms*)
		\[CapitalOmega]1,\[CapitalOmega]2,\[CapitalOmega]3,\[CapitalLambda]1,\[CapitalLambda]2,\[CapitalLambda]3,
		(*generators*)
		y,TSU2L,TSU3c,
		(*index names*)
		A,a,c,
		J,i,j,
		p
	}
	,
	
	(* define shorthand for generators *)
	y[\[Phi]_]:=FirstCase[GetFields[\[Phi]][Charges],ToExpression["U1Y"][x_]:>x,0,All];
	TSU2L[Jadj_,ifund_,jfund_]:= CG[gen[ToExpression["SU2L"][fund]],{Jadj,ifund,jfund}];
	TSU3c[Aadj_,\[Alpha]fund_,\[Beta]fund_]:=CG[gen[ToExpression["SU3c"][fund]],{Aadj,\[Alpha]fund,\[Beta]fund}];
	
	
	(* fermionic gauge interactions *)
	\[CapitalOmega]3= Sqrt[2] ToExpression["g3"][]RelabelIndices[
		Bar@ToExpression["q"][a,i,p]**PR**ToExpression["Gt"][A]TSU3c[A,a,c]ToExpression["qt"][c,i,p]-
		Bar@ToExpression["u"][a,p]**PL**ToExpression["Gt"][A]TSU3c[A,a,c]ToExpression["ut"][c,p]-
		Bar@ToExpression["d"][a,p]**PL**ToExpression["Gt"][A]TSU3c[A,a,c]ToExpression["dt"][c,p]
	];
	
	\[CapitalOmega]2= Sqrt[2] ToExpression["g2"][]RelabelIndices[
		Bar@ToExpression["q"][a,i,p]**PR**ToExpression["Wt"][J]TSU2L[J,i,j]ToExpression["qt"][a,j,p]+
		Bar@ToExpression["l"][i,p]**PR**ToExpression["Wt"][J]TSU2L[J,i,j]ToExpression["lt"][j,p]+
		Bar@ToExpression["Hut"][i]**PR**ToExpression["Wt"][J]TSU2L[J,i,j]ToExpression["Hu"][j]+
		Bar@ToExpression["Hdt"][i]**PR**ToExpression["Wt"][J]TSU2L[J,i,j]ToExpression["Hd"][j]
	];
	
	\[CapitalOmega]1= Sqrt[2] ToExpression["g1"][]RelabelIndices[
		y[ToExpression["q"]]Bar@ToExpression["q"][a,i,p]**PR**ToExpression["Bt"][]ToExpression["qt"][a,i,p]-
		y[ToExpression["u"]]Bar@ToExpression["u"][a,p]**PL**ToExpression["Bt"][]ToExpression["ut"][a,p]-
		y[ToExpression["d"]]Bar@ToExpression["d"][a,p]**PL**ToExpression["Bt"][]ToExpression["dt"][a,p]+
		y[ToExpression["l"]]Bar@ToExpression["l"][i,p]**PR**ToExpression["Bt"][]ToExpression["lt"][i,p]-
		y[ToExpression["e"]]Bar@ToExpression["e"][p]**PL**ToExpression["Bt"][]ToExpression["et"][p]+
		y[ToExpression["Hu"]]Bar@ToExpression["Hut"][i]**PR**ToExpression["Bt"][]ToExpression["Hu"][i]+
		y[ToExpression["Hd"]]Bar@ToExpression["Hdt"][i]**PR**ToExpression["Bt"][]ToExpression["Hd"][i]
	];
	
	
	(* scalar gauge interactions *)
	\[CapitalLambda]3[A_]:=ToExpression["g3"][]RelabelIndices[
		Bar@ToExpression["qt"][a,i,p]TSU3c[A,a,c]ToExpression["qt"][c,i,p]-
		Bar@ToExpression["ut"][a,p]TSU3c[A,a,c]ToExpression["ut"][c,p]-
		Bar@ToExpression["dt"][a,p]TSU3c[A,a,c]ToExpression["dt"][c,p]
	,
		Unique->True
	];
	
	\[CapitalLambda]2[I_]:=ToExpression["g2"][]RelabelIndices[
		Bar@ToExpression["qt"][a,i,p]TSU2L[I,i,j]ToExpression["qt"][a,j,p]+
		Bar@ToExpression["lt"][i,p]TSU2L[I,i,j]ToExpression["lt"][j,p]+
		Bar@ToExpression["Hu"][i]TSU2L[I,i,j]ToExpression["Hu"][j]+
		Bar@ToExpression["Hd"][i]TSU2L[I,i,j]ToExpression["Hd"][j]
	,
		Unique->True
	];
	
	\[CapitalLambda]1[]:=ToExpression["g1"][]RelabelIndices[
		y[ToExpression["q"]]Bar@ToExpression["qt"][a,i,p]ToExpression["qt"][a,i,p]-
		y[ToExpression["u"]]Bar@ToExpression["ut"][a,p]ToExpression["ut"][a,p]-
		y[ToExpression["d"]]Bar@ToExpression["dt"][a,p]ToExpression["dt"][a,p]+
		y[ToExpression["l"]]Bar@ToExpression["lt"][i,p]ToExpression["lt"][i,p]-
		y[ToExpression["e"]]Bar@ToExpression["et"][p]ToExpression["et"][p]+
		y[ToExpression["Hu"]]Bar@ToExpression["Hu"][i]ToExpression["Hu"][i]+
		y[ToExpression["Hd"]]Bar@ToExpression["Hd"][i]ToExpression["Hd"][i]
	,
		Unique->True
	];
	
	
	(* full gauge interactions *)
	-RelabelIndices[
		PlusHc[\[CapitalOmega]3+\[CapitalOmega]2+\[CapitalOmega]1]+
		1/2 \[CapitalLambda]3[A]\[CapitalLambda]3[A]+1/2 \[CapitalLambda]2[J]\[CapitalLambda]2[J]+1/2 \[CapitalLambda]1[]\[CapitalLambda]1[]
	]
]


(* ::Subsection::Closed:: *)
(*Superpotential (& F-terms)*)


SuperPotentialMSSM[]:=Hold@Module[
	{
		\[CurlyEpsilon],
		\[ScriptCapitalL]Wscalar, \[ScriptCapitalL]Wfermion,
		(*index names*)
		a,c,
		i,j,
		p,r,s,t
	}
	,
	(* define SUSubscript[(2), L] \[Epsilon] tensor *)
	\[CurlyEpsilon][i_,j_]:=CG[eps@ToExpression["SU2L"],{i,j}];
	
	
	(* scalar contribution to superpotential *)
	\[ScriptCapitalL]Wscalar = -RelabelIndices@Expand[
     	(* scalar bilinears *)
     	Bar@ToExpression["\[Mu]t"][] ToExpression["\[Mu]t"][] (Bar@ToExpression["Hu"][i] ToExpression["Hu"][i] + Bar@ToExpression["Hd"][i] ToExpression["Hd"][i]) -
      	(* scalar trilinears *)
      	PlusHc[Bar@ToExpression["\[Mu]t"][] ( 
         		Bar@ToExpression["Hd"][i] Bar@ToExpression["ut"][a, p] Bar@ToExpression["yu"][r, p] ToExpression["qt"][a, i, r] +
          		Bar@ToExpression["Hu"][i] Bar@ToExpression["dt"][a, p] Bar@ToExpression["yd"][r, p] ToExpression["qt"][a, i, r] +
          		Bar@ToExpression["Hu"][i] Bar@ToExpression["et"][p] Bar@ToExpression["ye"][r, p] ToExpression["lt"][i, r]
         	)] +
      	(* 2 sfermions and 2 Higgs *)
      	(Bar@ToExpression["Hu"][i] ToExpression["Hu"][i]) (Bar@ToExpression["ut"][a, p] Bar@ToExpression["yu"][r, p] ToExpression["yu"][r, s] ToExpression["ut"][a, s]) +
      	(Bar@ToExpression["Hd"][i] ToExpression["Hd"][i]) (Bar@ToExpression["dt"][a, p] Bar@ToExpression["yd"][r, p] ToExpression["yd"][r, s] ToExpression["dt"][a, s]) +
      	(Bar@ToExpression["Hd"][i] ToExpression["Hd"][i]) (Bar@ToExpression["et"][p] Bar@ToExpression["ye"][r, p] ToExpression["ye"][r, s] ToExpression["et"][s]) +
      	(Bar@ToExpression["Hd"][i] ToExpression["Hd"][i]) (Bar@ToExpression["lt"][j, p] ToExpression["ye"][p, r] Bar@ToExpression["ye"][s, r] ToExpression["lt"][j, s]) +
      	(Bar@ToExpression["Hd"][i] ToExpression["Hd"][i]) (Bar@ToExpression["qt"][a, j, p] ToExpression["yd"][p, r] Bar@ToExpression["yd"][s, r] ToExpression["qt"][a, j, s]) +
      	(Bar@ToExpression["Hu"][i] ToExpression["Hu"][i]) (Bar@ToExpression["qt"][a, j, p] ToExpression["yu"][p, r] Bar@ToExpression["yu"][s, r] ToExpression["qt"][a, j, s]) -
      	PlusHc[(Bar@ToExpression["Hd"][i] ToExpression["Hu"][i]) (Bar@ToExpression["ut"][a, p] Bar@ToExpression["yu"][r, p] ToExpression["yd"][r, s] ToExpression["dt"][a, s])] -
      	Bar@ToExpression["qt"][a, i, p] ToExpression["Hu"][i] ToExpression["yu"][p, r] Bar@ToExpression["yu"][s, r] Bar@ToExpression["Hu"][j] ToExpression["qt"][a, j, s] -
      	Bar@ToExpression["qt"][a, i, p] ToExpression["Hd"][i] ToExpression["yd"][p, r] Bar@ToExpression["yd"][s, r] Bar@ToExpression["Hd"][j] ToExpression["qt"][a, j, s] -
      	Bar@ToExpression["lt"][i, p] ToExpression["Hd"][i] ToExpression["ye"][p, r] Bar@ToExpression["ye"][s, r] Bar@ToExpression["Hd"][j] ToExpression["lt"][j, s] +
      	(* 4 sfermions *)
      	(Bar@ToExpression["ut"][a, p] Bar@ToExpression["yu"][r, p] ToExpression["qt"][a, i, r]) Bar[Bar@ToExpression["ut"][c, s] Bar@ToExpression["yu"][t, s] ToExpression["qt"][c, i, t]] +
      	(Bar@ToExpression["dt"][a, p] Bar@ToExpression["yd"][r, p] ToExpression["qt"][a, i, r]) Bar[Bar@ToExpression["dt"][c, s] Bar@ToExpression["yd"][t, s] ToExpression["qt"][c, i, t]] +
      	(Bar@ToExpression["et"][p] Bar@ToExpression["ye"][r, p] ToExpression["lt"][i, r]) Bar[Bar@ToExpression["et"][s] Bar@ToExpression["ye"][t, s] ToExpression["lt"][i, t]] +
      	PlusHc[(Bar@ToExpression["dt"][a, p] Bar@ToExpression["yd"][r, p] ToExpression["qt"][a, i, r]) Bar[Bar@ToExpression["et"][s] Bar@ToExpression["ye"][t, s] ToExpression["lt"][i, t]]]
     ];
     
     
     (* fermonic contribution to superpotential *)
     \[ScriptCapitalL]Wfermion=-PlusHc@RelabelIndices[
		(* Higgsino mass *)
		ToExpression["\[Mu]t"][]Bar@\[CurlyEpsilon][i,j]Bar@CConj@ToExpression["Hut"][i]**PL**ToExpression["Hdt"][j]+
		(* SM Yukawas *)
		Bar@ToExpression["yu"][r,p]Bar@ToExpression["u"][a,p]**PL**ToExpression["q"][a,i,r]Bar@\[CurlyEpsilon][i,j]ToExpression["Hu"][j]-
		Bar@ToExpression["yd"][r,p]Bar@ToExpression["d"][a,p]**PL**ToExpression["q"][a,i,r]Bar@\[CurlyEpsilon][i,j]ToExpression["Hd"][j]-
		Bar@ToExpression["ye"][r,p]Bar@ToExpression["e"][p]**PL**ToExpression["l"][i,r]Bar@\[CurlyEpsilon][i,j]ToExpression["Hd"][j]+
		(* Remaining trilinears *)
		Bar@ToExpression["yu"][r,p]Bar@ToExpression["u"][a,p]**PL**ToExpression["Hut"][j]ToExpression["qt"][a,i,r]Bar@\[CurlyEpsilon][i,j]-
		Bar@ToExpression["yd"][r,p]Bar@ToExpression["d"][a,p]**PL**ToExpression["Hdt"][j]ToExpression["qt"][a,i,r]Bar@\[CurlyEpsilon][i,j]-
		Bar@ToExpression["ye"][r,p]Bar@ToExpression["e"][p]**PL**ToExpression["Hdt"][j]ToExpression["lt"][i,r]Bar@\[CurlyEpsilon][i,j]+
		(**)
		Bar@ToExpression["yu"][r,p]Bar@ToExpression["ut"][a,p]Bar@CConj@ToExpression["q"][a,i,r]**PL**ToExpression["Hut"][j]Bar@\[CurlyEpsilon][i,j]-
		Bar@ToExpression["yd"][r,p]Bar@ToExpression["dt"][a,p]Bar@CConj@ToExpression["q"][a,i,r]**PL**ToExpression["Hdt"][j]Bar@\[CurlyEpsilon][i,j]-
		Bar@ToExpression["ye"][r,p]Bar@ToExpression["et"][p]Bar@CConj@ToExpression["l"][i,r]**PL**ToExpression["Hdt"][j]Bar@\[CurlyEpsilon][i,j]
	];
	
	
	(* full super potential *)
	\[ScriptCapitalL]Wscalar+\[ScriptCapitalL]Wfermion
]


(* ::Subsection::Closed:: *)
(*Soft SUSY breaking terms*)


SoftSUSYbreakingMSSM[]:=Hold@Module[
	{\[CurlyEpsilon], \[ScriptCapitalL]SoftBreak,
	i,j,p,r,a}
	,
	(* define SUSubscript[(2), L] \[Epsilon] tensor *)
	\[CurlyEpsilon][i_,j_]:=CG[eps@ToExpression["SU2L"],{i,j}];
	
	
	(* soft SUSY breaking *)
	\[ScriptCapitalL]SoftBreak= -RelabelIndices@PlusHc[
		ToExpression["b"][]Bar@\[CurlyEpsilon][i,j]ToExpression["Hu"][i]ToExpression["Hd"][j]+
		Bar@ToExpression["qt"][a,i,p]ToExpression["au"][p,r]ToExpression["ut"][a,r]\[CurlyEpsilon][i,j]Bar@ToExpression["Hu"][j]-
		Bar@ToExpression["qt"][a,i,p]ToExpression["ad"][p,r]ToExpression["dt"][a,r]\[CurlyEpsilon][i,j]Bar@ToExpression["Hd"][j]-
		Bar@ToExpression["lt"][i,p]ToExpression["ae"][p,r]ToExpression["et"][r]\[CurlyEpsilon][i,j]Bar@ToExpression["Hd"][j]
	];
	
	\[ScriptCapitalL]SoftBreak
]


(* ::Section:: *)
(*Rotation to the mass basis*)


(* ::Subsection::Closed:: *)
(*Rotation to Higgs basis [not in use]*)


ToHiggsBasis[L_]:=Hold@Module[
	{
		Lag= L,
		To2HDM, ToHiggsBasis,
		mLight, mHeavy, mLightHeavy, mHeavyLight,
		\[CurlyEpsilon],
		j
	}
	,
	(* define anti-symmetric SUSubscript[(2), L] tensor *)
	\[CurlyEpsilon][i_,j_]:=CG[eps@ToExpression["SU2L"],{i,j}];
	
	(* * * * * * * * * * *)
	
	(* define conjugate of Hd *)
	DefineField[ToExpression["Hdc"],Scalar,
		Indices->{ToExpression["SU2L"][fund]},
		Charges->{ToExpression["U1Y"][1/2]},
		Mass->Heavy
	];
	
	(* define light Higgs doublet to be idetified with SMEFT Higgs doublet *)
	DefineField["H",Scalar,
		Indices->{"SU2L"[fund]},
		Charges->{"U1Y"[1/2]},
		Mass->0
	];
	
	(* define heavy doublet to be integarated out *)
	DefineField[ToExpression["\[CapitalPhi]"],Scalar,
		Indices->{ToExpression["SU2L"][fund]},
		Charges->{ToExpression["U1Y"][1/2]},
		Mass->{Heavy,ToExpression["$M\[CapitalPhi]"]}
	];
	
	(* * * * * * * * * * *)
	
	(* rule to rotate from Hd to Hdc to obtain the usual 2HDM Lagrangian *)
	To2HDM={Field[ToExpression["Hd"],Scalar,{Index[i_,ToExpression["SU2L"][fund]]},devs_]:>RelabelIndices[-\[CurlyEpsilon][i,j]CD[devs,Bar@ToExpression["Hdc"][j]],Unique->True]};
	
	(* rule to rotate to Higgs basis *)
	ToHiggsBasis={
		Field[ToExpression["Hu"],Scalar,{Index[i_,ToExpression["SU2L"][fund]]},devs_]:>(Sin[ToExpression[ToExpression["\[Beta]"]]] CD[devs,ToExpression["H"][i]]+Cos[ToExpression["\[Beta]"]] CD[devs,ToExpression["\[CapitalPhi]"][i]]),
		Field[ToExpression["Hdc"],Scalar,{Index[i_,ToExpression["SU2L"][fund]]},devs_]:>(-Cos[ToExpression["\[Beta]"]] CD[devs,ToExpression["H"][i]]+Sin[ToExpression["\[Beta]"]] CD[devs,ToExpression["\[CapitalPhi]"][i]])
	};
	
	(* apply rotations to Higgs basis*)
	Lag= CollectOperators[
		TrigReduce@ RelabelIndices@ Contract@ ContractCGs[Lag/.To2HDM/.ToHiggsBasis]
	];
	
	Lag= CollectOperators[
		TrigReduce/@Lag,
		Matchete`PackageScope`NormalForm->False
	];
	
	(* * * * * * * * * * *)
	
	(* define masses in Higgs basis *)
	DefineCoupling[ToExpression["mH2"], EFTOrder->2, SelfConjugate->True, NiceForm->"\!\(\*SubsuperscriptBox[\(m\), \(H\), \(2\)]\)"];
	DefineCoupling[ToExpression["m\[CapitalPhi]"], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(m\), \(\[CapitalPhi]\)]\)"];
	
	(* extract current masses *)
	(* light mass *)
	mLight=CollectOperators[
		SelectOperatorClass[Lag,{Bar@ToExpression["H"],ToExpression["H"]},0],
		Matchete`PackageScope`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	(*heavy mass*)
	mHeavy=CollectOperators[
		SelectOperatorClass[Lag,{Bar@ToExpression["\[CapitalPhi]"],ToExpression["\[CapitalPhi]"]},0],
		Matchete`PackageScope`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	(*mixed mass*)
	mLightHeavy=CollectOperators[
		HcSimplify@SelectOperatorClass[Lag,{Bar@ToExpression["\[CapitalPhi]"],ToExpression["H"]},0]/.HcTerms->Identity,
		Matchete`PackageScope`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	mHeavyLight=CollectOperators[
		HcSimplify@SelectOperatorClass[Lag,{Bar@ToExpression["\[CapitalPhi]"],ToExpression["H"]},0]/.HcTerms->Bar,
		Matchete`PackageScope`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	
	(* redefine all the Higgs masses *)
	Lag= RelabelIndices@Matchete`PackageScope`NormalForm[
		Lag/.{
			mLight->-ToExpression["mH2"][],
			mHeavy->-ToExpression["m\[CapitalPhi]"][]^2,
			mLightHeavy|mHeavyLight->+ToExpression["mH2"][] * Tan[2*ToExpression["\[Beta]"]]
		}
	];
	
	Lag
]


(* ::Subsection:: *)
(*Rotation to mass basis for the Higgs doublets*)


ToHiggsDoubletMassBasis[L_]:=Hold@Module[
	{
		Lag= L,
		To2HDM, ToHiggsDoubletMassBasis,
		mLight, mHeavy, mLightHeavy, mHeavyLight,
		\[CurlyEpsilon],
		j
	}
	,
	(* define anti-symmetric SUSubscript[(2), L] tensor *)
	\[CurlyEpsilon][i_,j_]:=CG[eps@ToExpression["SU2L"],{i,j}];
	
	(* * * * * * * * * * *)
	
	(* define conjugate of Hd *)
	DefineField[ToExpression["Hdc"],Scalar,
		Indices->{ToExpression["SU2L"][fund]},
		Charges->{ToExpression["U1Y"][1/2]},
		Mass->Heavy
	];
	
	(* define light Higgs doublet to be idetified with SMEFT Higgs doublet *)
	DefineField[ToExpression["H"],Scalar,
		Indices->{ToExpression["SU2L"][fund]},
		Charges->{ToExpression["U1Y"][1/2]},
		Mass->0
	];
	
	(* define heavy doublet *)
	DefineField[ToExpression["\[CapitalPhi]"],Scalar,
		Indices->{ToExpression["SU2L"][fund]},
		Charges->{ToExpression["U1Y"][1/2]},
		Mass->{Heavy,ToExpression["$M\[CapitalPhi]"]}
	];
	
	(* * * * * * * * * * *)
	
	(* rule to rotate from Hd to Hdc to obtain the usual 2HDM Lagrangian *)
	To2HDM={Field[ToExpression["Hd"],Scalar,{Index[i_,ToExpression["SU2L"][fund]]},devs_]:>RelabelIndices[-\[CurlyEpsilon][i,j]CD[devs,Bar@ToExpression["Hdc"][j]],Unique->True]};
	
	(* rule to rotate to doublet mass basis *)
	ToHiggsDoubletMassBasis={
		Field[ToExpression["Hu"],Scalar,{Index[i_,ToExpression["SU2L"][fund]]},devs_]:>(Sin[ToExpression[ToExpression["\[Beta]"]]] CD[devs,ToExpression["H"][i]]+Cos[ToExpression["\[Beta]"]] CD[devs,ToExpression["\[CapitalPhi]"][i]]),
		Field[ToExpression["Hdc"],Scalar,{Index[i_,ToExpression["SU2L"][fund]]},devs_]:>(-Cos[ToExpression["\[Beta]"]] CD[devs,ToExpression["H"][i]]+Sin[ToExpression["\[Beta]"]] CD[devs,ToExpression["\[CapitalPhi]"][i]])
	};
	
	(* apply rotations to doublet mass basis*)
	Lag= CollectOperators[
		TrigReduce@ RelabelIndices@ Contract@ ContractCGs[Lag/.To2HDM/.ToHiggsDoubletMassBasis]
	];

	Lag= CollectOperators[
		TrigReduce/@Lag,
		Matchete`PackageScope`NormalForm->False
	];
	
	(* * * * * * * * * * *)
	
	(* define masses in Higgs basis *)
	DefineCoupling[ToExpression["mH2"], EFTOrder->2, SelfConjugate->True, NiceForm->"\!\(\*SubsuperscriptBox[\(m\), \(H\), \(2\)]\)"];
	DefineCoupling[ToExpression["m\[CapitalPhi]"], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(m\), \(\[CapitalPhi]\)]\)"];

	(* extract current masses *)
	(* light mass *)
	mLight=CollectOperators[
		SelectOperatorClass[Lag,{Bar@ToExpression["H"],ToExpression["H"]},0],
		Matchete`PackageScope`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	(*heavy mass*)
	mHeavy=CollectOperators[
		SelectOperatorClass[Lag,{Bar@ToExpression["\[CapitalPhi]"],ToExpression["\[CapitalPhi]"]},0],
		Matchete`PackageScope`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	(*mixed mass*)
	mLightHeavy=CollectOperators[
		HcSimplify@SelectOperatorClass[Lag,{Bar@ToExpression["\[CapitalPhi]"],ToExpression["H"]},0]/.HcTerms->Identity,
		Matchete`PackageScope`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	mHeavyLight=CollectOperators[
		HcSimplify@SelectOperatorClass[Lag,{Bar@ToExpression["\[CapitalPhi]"],ToExpression["H"]},0]/.HcTerms->Bar,
		Matchete`PackageScope`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	
	(* redefine all the Higgs masses *)
	Lag= RelabelIndices@Matchete`PackageScope`NormalForm[
		Lag/.{
			mLight->-ToExpression["mH2"][],
			mHeavy->-ToExpression["m\[CapitalPhi]"][]^2,
			mLightHeavy|mHeavyLight->0
		}
	];

	Lag
]


(* ::Subsection::Closed:: *)
(*Rotation to Higgsino mass basis*)


ToHiggsinoMassBasis[L_]:=Hold@Module[
	{
		Lag=L,
		RotateHiggsinos, CombineToVectorlikeFermion,
		\[CurlyEpsilon],
		j
	}
	,
	(* define anti-symmetric SUSubscript[(2), L] tensor *)
	\[CurlyEpsilon][i_,j_]:=CG[eps@ToExpression["SU2L"],{i,j}];
	
	(* * * * * * * * * * *)
	
	(* introduce a vectorlike fermion *)
	
	(* define vectorlike fermion *)
	DefineField[
		ToExpression["\[CapitalSigma]"],
		Fermion,
		Indices->{ToExpression["SU2L"][fund]},
		Charges->{ToExpression["U1Y"][1/2]},
		Mass->0 (* the mass is defined manually *)
	];
	
	(* replace Hut & Hdt by \[CapitalSigma] *)
	RotateHiggsinos={
		Field[ToExpression["Hut"],Fermion,{Index[i_,ToExpression["SU2L"][fund]]},dev_List]:>PL**CD[dev,ToExpression["\[CapitalSigma]"][i]],
		Field[ToExpression["Hdt"],Fermion,{Index[i_,ToExpression["SU2L"][fund]]},dev_List]:>RelabelIndices[-\[CurlyEpsilon][i,j]PL**CD[dev,CConj@ToExpression["\[CapitalSigma]"][j]],Unique->True]
	};
	
	(* combine chiral to vectorlike kinetic and mass terms *)
	CombineToVectorlikeFermion={
		(* kinetic terms *)
		Evaluate[
			I Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{i_},{}]]**DiracProduct[GammaM[\[Mu]_]]**PL**Field[ToExpression["\[CapitalSigma]"],Fermion,{i_},{\[Mu]_}]+
			I Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{j_},{}]]**DiracProduct[GammaM[\[Nu]_]]**PR**Field[ToExpression["\[CapitalSigma]"],Fermion,{j_},{\[Nu]_}]
		]:>I Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{i},{}]]**DiracProduct[GammaM[\[Nu]]]**Field[ToExpression["\[CapitalSigma]"],Fermion,{i},{\[Nu]}]
		,
		Evaluate[
			I Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{i_},{}]]**DiracProduct[GammaM[\[Mu]_]]**PL**Field[ToExpression["\[CapitalSigma]"],Fermion,{i_},{\[Mu]_}]-
			I Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{j_},{\[Nu]_}]]**DiracProduct[GammaM[\[Nu]_]]**PR**Field[ToExpression["\[CapitalSigma]"],Fermion,{j_},{}]
		]:>I Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{i},{}]]**DiracProduct[GammaM[\[Nu]]]**Field[ToExpression["\[CapitalSigma]"],Fermion,{i},{\[Nu]}]
		,
		Evaluate[
			-I Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{i_},{\[Mu]_}]]**DiracProduct[GammaM[\[Mu]_]]**PL**Field[ToExpression["\[CapitalSigma]"],Fermion,{i_},{}]+
			I Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{j_},{}]]**DiracProduct[GammaM[\[Nu]_]]**PR**Field[ToExpression["\[CapitalSigma]"],Fermion,{j_},{\[Nu]_}]
		]:>I Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{i},{}]]**DiracProduct[GammaM[\[Nu]]]**Field[ToExpression["\[CapitalSigma]"],Fermion,{i},{\[Nu]}]
		,
		(* mass terms *)
		Evaluate[
			Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{i_},{}]]**DiracProduct[Proj[-1]]**Field[ToExpression["\[CapitalSigma]"],Fermion,{i_},{}]
		]:>(Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{i},{}]]**Field[ToExpression["\[CapitalSigma]"],Fermion,{i},{}]-Bar[Field[ToExpression["\[CapitalSigma]"],Fermion,{i},{}]]**DiracProduct[Proj[+1]]**Field[ToExpression["\[CapitalSigma]"],Fermion,{i},{}])
	};
	
	(* * * * * * * * * * *)
	
	(* apply the rules introducing \[CapitalSigma] *)
	Lag = RelabelIndices@Contract@ContractCGs[Lag /. RotateHiggsinos];
	
	(* vectoralize mass terms*)
	Lag = Expand[Lag//.CombineToVectorlikeFermion];
	
	Lag
]


(* ::Subsection::Closed:: *)
(*Replace rotation angles*)


(* ::Text:: *)
(*It is more efficient to define several different couplings for sin(\[Beta]), sin(2\[Beta]), sin(4\[Beta]), cos(\[Beta]), cos(2\[Beta]), cos(4\[Beta]). *)
(*If only a couplings \[Beta] would be defined, the coefficients would become complicated functions of \[Beta].*)


ReplaceRotationAngles[L_]:=Hold@Module[
	{Lag = L, ReplaceTrig\[Beta]}
	,
	(* define couplings for all trigonomatric functions appearing in the expression *)
	DefineCoupling[ToExpression["s\[Beta]"],  EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(s\), \(\[Beta]\)]\)"];
	DefineCoupling[ToExpression["c\[Beta]"],  EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(c\), \(\[Beta]\)]\)"];
	DefineCoupling[ToExpression["s2\[Beta]"], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(s\), \(2\[Beta]\)]\)"];
	DefineCoupling[ToExpression["c2\[Beta]"], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(c\), \(2\[Beta]\)]\)"];
	DefineCoupling[ToExpression["s4\[Beta]"], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(s\), \(4\[Beta]\)]\)"];
	DefineCoupling[ToExpression["c4\[Beta]"], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(c\), \(4\[Beta]\)]\)"];
	
	(* define replacement rules *)
	ReplaceTrig\[Beta]={
		Sin[ToExpression["\[Beta]"]]  -> ToExpression["s\[Beta]"][],
		Cos[ToExpression["\[Beta]"]]  -> ToExpression["c\[Beta]"][],
		Sin[2ToExpression["\[Beta]"]] -> ToExpression["s2\[Beta]"][],
		Cos[2ToExpression["\[Beta]"]] -> ToExpression["c2\[Beta]"][],
		Sin[4ToExpression["\[Beta]"]] -> ToExpression["s4\[Beta]"][],
		Cos[4ToExpression["\[Beta]"]] -> ToExpression["c4\[Beta]"][],
		Tan[2ToExpression["\[Beta]"]] :> (
			DefineCoupling[ToExpression["t2\[Beta]"], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(t\), \(2\[Beta]\)]\)"];
			ToExpression["t2\[Beta]"][]
		)
	};
	
	(* apply replacement rules *)
	Lag = Contract@ContractCGs[Lag /.ReplaceTrig\[Beta]];
	
	Lag
]


TrigSimplifiy[expr_]:=Hold[expr/.{ToExpression["s\[Beta]"][]*ToExpression["c\[Beta]"][]->ToExpression["s2\[Beta]"][]/2}]


(* ::Section:: *)
(*Cleaning up*)


(* ::Subsection::Closed:: *)
(*Remove temporarily defined fields and couplings*)


RemoveTMPdef[]:=Hold[
	(* remove all fields defined intermediately *)
	RemoveField[ToExpression["Hu"]];
	RemoveField[ToExpression["Hd"]];
	RemoveField[ToExpression["Hdc"]];
	RemoveField[ToExpression["Hut"]];
	RemoveField[ToExpression["Hdt"]];
	RemoveCoupling[ToExpression["MHu"]];
	RemoveCoupling[ToExpression["MHd"]];
	
	(* correct the masses associated to fields *)
	Matchete`PackageScope`$FieldAssociation[ToExpression["\[CapitalPhi]"]][Mass]  = ToExpression["m\[CapitalPhi]"];
	RemoveCoupling[ToExpression["$M\[CapitalPhi]"]];
	Matchete`PackageScope`$FieldAssociation[ToExpression["\[CapitalSigma]"]][Mass]  = ToExpression["\[Mu]t"];
	Matchete`PackageScope`$FieldAssociation[ToExpression["\[CapitalSigma]"]][Heavy] = True;
]


(* ::Section:: *)
(*End*)


End[];
