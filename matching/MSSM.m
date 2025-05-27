(* ::Package:: *)

(* ::Title:: *)
(*MSSM model file*)


(* ::Subtitle:: *)
(*Works  on  Matchete  git commit "46ec78a7":*)
(*[ https://gitlab.com/matchete/matchete/-/commit/46ec78a79dfd73d1ef238dfcea8df37a3b8744be ]*)


(* ::Chapter:: *)
(*Definitions*)


(* ::Section:: *)
(*Gauge groups*)


(* ::Subsection::Closed:: *)
(*QCD*)


DefineGaugeGroup[SU3c, SU@3, g3, G, 
	FundAlphabet-> CharacterRange["a", "f"], 
	AdjAlphabet -> CharacterRange["A", "F"], 
	NiceForm->{"\!\(\*SubscriptBox[\(g\), \(3\)]\)", Default}
];


(* ::Subsection::Closed:: *)
(*Weak*)


DefineGaugeGroup[SU2L, SU@2, g2, W, 
	FundAlphabet -> CharacterRange["i", "n"], 
	AdjAlphabet -> CharacterRange["I", "N"], 
	NiceForm->{"\!\(\*SubscriptBox[\(g\), \(2\)]\)", Default}
];


(* ::Subsection::Closed:: *)
(*Hypercharge*)


DefineGaugeGroup[U1Y, U1, g1, B, 
	NiceForm->{"\!\(\*SubscriptBox[\(g\), \(1\)]\)", Default}
];


(* ::Section:: *)
(*Flavor index*)


(* ::Subsection::Closed:: *)
(*Flavor index*)


ParameterDefault[Nf-> 3]


DefineFlavorIndex[Flavor, Nf,
	IndexAlphabet-> {"p", "r", "s", "t", "u", "v"}
];


(* ::Section:: *)
(*Fields*)


(* ::Subsection::Closed:: *)
(*SM fermions*)


DefineField[q, Fermion, 
	Indices -> {SU3c[fund], SU2L[fund], Flavor}, 
	Charges -> {U1Y[1/6]}, 
	Chiral  -> LeftHanded, 
	Mass    -> 0
];

DefineField[u, Fermion,
	Indices -> {SU3c[fund], Flavor},
	Charges -> {U1Y[2/3]}, 
	Chiral  -> RightHanded,
	Mass    -> 0
];
	
DefineField[d, Fermion,
	Indices -> {SU3c[fund], Flavor},
	Charges -> {U1Y[-1/3]},
	Chiral  -> RightHanded,
	Mass    -> 0
];
	
DefineField[l, Fermion,
	Indices -> {SU2L[fund], Flavor},
	Charges -> {U1Y[-1/2]},
	Chiral  -> LeftHanded, 
	Mass    -> 0
];
	
DefineField[e, Fermion,
	Indices -> {Flavor},
	Charges -> {U1Y[-1]},
	Chiral-> RightHanded,
	Mass-> 0
];


(* ::Subsection::Closed:: *)
(*Squarks and sleptons*)


DefineField[qt, Scalar,
	Indices  -> {SU3c[fund], SU2L[fund], Flavor},
	Charges  -> {U1Y[1/6]},
	Mass     -> {Heavy, mqt, {Flavor}},
	NiceForm -> {"\!\(\*OverscriptBox[\(q\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), OverscriptBox[\(q\), \(~\)]]\)"}
];

DefineField[ut, Scalar,
	Indices  -> {SU3c[fund], Flavor},
	Charges  -> {U1Y[2/3]},
	Mass     -> {Heavy, mut, {Flavor}},
	NiceForm -> {"\!\(\*OverscriptBox[\(u\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), OverscriptBox[\(u\), \(~\)]]\)"}
];

DefineField[dt, Scalar,
	Indices  -> {SU3c[fund], Flavor},
	Charges  -> {U1Y[-1/3]},
	Mass     -> {Heavy, mdt, {Flavor}},
	NiceForm -> {"\!\(\*OverscriptBox[\(d\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), OverscriptBox[\(d\), \(~\)]]\)"}
];
DefineField[lt, Scalar,
	Indices  -> {SU2L[fund], Flavor},
	Charges  -> {U1Y[-1/2]},
	Mass     -> {Heavy, mlt, {Flavor}},
	NiceForm -> {"\!\(\*OverscriptBox[\(l\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), OverscriptBox[\(l\), \(~\)]]\)"}
];

DefineField[et, Scalar,
	Indices  -> {Flavor},
	Charges  -> {U1Y[-1]},
	Mass     -> {Heavy, met, {Flavor}},
	NiceForm -> {"\!\(\*OverscriptBox[\(e\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), OverscriptBox[\(e\), \(~\)]]\)"}
];


(* ::Subsection::Closed:: *)
(*Gauginos*)


DefineField[Gt, Fermion,
	Indices       -> {SU3c[adj]},
	SelfConjugate -> True,
	Mass          -> {Heavy, m3},
	NiceForm      -> {"\!\(\*OverscriptBox[\(G\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), \(3\)]\)"}
];
DefineField[Wt, Fermion,
	Indices       -> {SU2L[adj]},
	SelfConjugate -> True,
	Mass          -> {Heavy, m2},
	NiceForm      -> {"\!\(\*OverscriptBox[\(W\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), \(2\)]\)"}
];
DefineField[Bt, Fermion,
	SelfConjugate -> True,
	Mass          -> {Heavy, m1},
	NiceForm      -> {"\!\(\*OverscriptBox[\(B\), \(~\)]\)","\!\(\*SubscriptBox[\(m\), \(1\)]\)"}
];


(* ::Subsection::Closed:: *)
(*Higgs sector (2HDM type-II)*)


(* ::Text:: *)
(*Higgs fields of initial Lagrangian Hu and Hd are removed in the end.*)


DefineField[Hu, Scalar,
	Indices  -> {SU2L[fund]},
	Charges  -> {U1Y[1/2]},
	Mass     -> {Heavy, MHu},
	NiceForm -> "\!\(\*SubscriptBox[\(H\), \(u\)]\)"
];

DefineField[Hd, Scalar,
	Indices  -> {SU2L[fund]},
	Charges  -> {U1Y[-1/2]},
	Mass     -> {Heavy, MHd},
	NiceForm -> "\!\(\*SubscriptBox[\(H\), \(d\)]\)"
];


(* ::Text:: *)
(*Conjugate of Hd defined for convenience and deleted eventually.*)


DefineField[Hdc, Scalar,
	Indices -> {SU2L[fund]},
	Charges -> {U1Y[1/2]},
	Mass    -> Heavy
];


(* ::Text:: *)
(*Definition of the light SMEFT-like Higgs doublet H and the heavy BSM-like Higgs doublet \[CapitalPhi]*)


DefineField[H, Scalar,
	Indices -> {SU2L[fund]},
	Charges -> {U1Y[1/2]},
	Mass    -> 0
];

DefineField[\[CapitalPhi], Scalar,
	Indices -> {SU2L[fund]},
	Charges -> {U1Y[1/2]},
	Mass    -> {Heavy, $M\[CapitalPhi]}
];


(* ::Subsection::Closed:: *)
(*Higgsinos*)


(* ::Text:: *)
(*Chiral Higgsino fields present in initial Lagrangian but removed later on*)


DefineField[Hut, Fermion,
	Indices  -> {SU2L[fund]},
	Charges  -> {U1Y[1/2]},
	Chiral   -> LeftHanded,
	Mass     -> 0, (* mass defined separately later *)
	NiceForm -> "\!\(\*SubscriptBox[OverscriptBox[\(H\), \(~\)], \(u\)]\)"
];
DefineField[Hdt, Fermion,
	Indices  -> {SU2L[fund]},
	Charges  -> {U1Y[-1/2]},
	Chiral   -> LeftHanded,
	Mass     -> 0, (* mass defined separately later *)
	NiceForm -> "\!\(\*SubscriptBox[OverscriptBox[\(H\), \(~\)], \(d\)]\)"
];


(* ::Text:: *)
(*Vectorlike Higgsino field present in the final MSSM Lagrangian*)


(* define vectorlike fermion *)
DefineField[\[CapitalSigma], Fermion,
	Indices -> {SU2L[fund]},
	Charges -> {U1Y[1/2]},
	Mass    -> 0 (* the mass is later defined manually *)
];


(* ::Section:: *)
(*Couplings*)


(* ::Subsection::Closed:: *)
(*Yukawas*)


DefineCoupling[yu,
	Indices  -> {Flavor, Flavor},
	NiceForm -> "\!\(\*SubscriptBox[\(y\), \(u\)]\)"
];

DefineCoupling[yd, 
	Indices  -> {Flavor, Flavor},
	NiceForm -> "\!\(\*SubscriptBox[\(y\), \(d\)]\)"
];

DefineCoupling[ye,
	Indices  -> {Flavor, Flavor},
	NiceForm -> "\!\(\*SubscriptBox[\(y\), \(e\)]\)"
];


(* ::Subsection::Closed:: *)
(*SUSY-Yukawas*)


DefineCoupling[au, 
	Indices  -> {Flavor, Flavor}, 
	EFTOrder -> 0, 
	NiceForm -> "\!\(\*SubscriptBox[\(a\), \(u\)]\)"
];

DefineCoupling[ad, 
	Indices  -> {Flavor, Flavor}, 
	EFTOrder -> 0, 
	NiceForm -> "\!\(\*SubscriptBox[\(a\), \(d\)]\)"
];

DefineCoupling[ae, 
	Indices  -> {Flavor, Flavor}, 
	EFTOrder -> 0, 
	NiceForm -> "\!\(\*SubscriptBox[\(a\), \(e\)]\)"
];


(* ::Subsection::Closed:: *)
(*Superpotential masses*)


DefineCoupling[\[Mu]t, 
	EFTOrder      -> 0, 
	SelfConjugate -> True,
	NiceForm      -> "\!\(\*OverscriptBox[\(\[Mu]\), \(~\)]\)"
];

DefineCoupling[b,
	EFTOrder      -> 0,
	SelfConjugate -> True
];


(* ::Subsection::Closed:: *)
(*Higgs masses*)


DefineCoupling[mH2, 
	EFTOrder      -> 2,
	SelfConjugate -> True,
	NiceForm      -> "\!\(\*SubsuperscriptBox[\(m\), \(H\), \(2\)]\)"
];
(* only for light fields, the mass the squared masses is directly defined *)
DefineCoupling[m\[CapitalPhi],
	EFTOrder      -> 0,
	SelfConjugate -> True,
	NiceForm      -> "\!\(\*SubscriptBox[\(m\), \(\[CapitalPhi]\)]\)"
];


(* ::Subsection::Closed:: *)
(*Rotation angle to soft-SUSY mass basis and to Higgs basis*)


(* ::Text:: *)
(*Define short-hand notation for trigonometric functions in soft-SUSY mass basis with rotation angle \[Gamma]: [default]*)


DefineCoupling[s\[Gamma],  EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(s\), \(\[Gamma]\)]\)"];
DefineCoupling[c\[Gamma],  EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(c\), \(\[Gamma]\)]\)"];
DefineCoupling[s2\[Gamma], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(s\), \(2\[Gamma]\)]\)"];
DefineCoupling[c2\[Gamma], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(c\), \(2\[Gamma]\)]\)"];
DefineCoupling[s4\[Gamma], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(s\), \(4\[Gamma]\)]\)"];
DefineCoupling[c4\[Gamma], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(c\), \(4\[Gamma]\)]\)"];


(* ::Text:: *)
(*Define short-hand notation for trigonometric functions in Higgs basis with rotation angle \[Beta]: [not in use by default]*)


DefineCoupling[s\[Beta],  EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(s\), \(\[Beta]\)]\)"];
DefineCoupling[c\[Beta],  EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(c\), \(\[Beta]\)]\)"];
DefineCoupling[s2\[Beta], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(s\), \(2\[Beta]\)]\)"];
DefineCoupling[c2\[Beta], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(c\), \(2\[Beta]\)]\)"];
DefineCoupling[s4\[Beta], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(s\), \(4\[Beta]\)]\)"];
DefineCoupling[c4\[Beta], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(c\), \(4\[Beta]\)]\)"];


DefineCoupling[t2\[Beta], EFTOrder->0, SelfConjugate->True, NiceForm->"\!\(\*SubscriptBox[\(t\), \(2\[Beta]\)]\)"];


(* ::Chapter:: *)
(*Functions for constructing the Lagrangian*)


(* ::Section:: *)
(*Lagrangian in MSSM basis*)


(* ::Subsection::Closed:: *)
(*Free Lagrangian*)


MSSM`FreeLagrangianMSSM[] := FreeLag[
	G   , W  , B  ,
	Gt  , Wt , Bt ,
	q   , u  , d  , l  , e  ,
	qt  , ut , dt , lt , et ,
	Hu  , Hd ,
	Hut , Hdt
]


(* ::Subsection::Closed:: *)
(*Gauge interactions*)


(* ::Text:: *)
(*Only gauge interactions that are not encoded in covariant derivatives*)


MSSM`GaugeLagrangianMSSM[] := Module[
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
	y[\[Phi]_]                      := FirstCase[GetFields[\[Phi]][Charges], U1Y[x_]:>x, 0, All];
	TSU2L[Jadj_,ifund_,jfund_] := CG[gen[SU2L@fund], {Jadj,ifund,jfund}];
	TSU3c[Aadj_,\[Alpha]fund_,\[Beta]fund_] := CG[gen[SU3c@fund], {Aadj,\[Alpha]fund,\[Beta]fund}];
	
	
	(* fermionic gauge interactions *)
	\[CapitalOmega]3= Sqrt[2] g3[]RelabelIndices[
		Bar@q[a,i,p]**PR**Gt[A]TSU3c[A,a,c]qt[c,i,p]-
		Bar@u[a,p]**PL**Gt[A]TSU3c[A,a,c]ut[c,p]-
		Bar@d[a,p]**PL**Gt[A]TSU3c[A,a,c]dt[c,p]
	];
	
	\[CapitalOmega]2= Sqrt[2] g2[]RelabelIndices[
		Bar@q[a,i,p]**PR**Wt[J]TSU2L[J,i,j]qt[a,j,p]+
		Bar@l[i,p]**PR**Wt[J]TSU2L[J,i,j]lt[j,p]+
		Bar@Hut[i]**PR**Wt[J]TSU2L[J,i,j]Hu[j]+
		Bar@Hdt[i]**PR**Wt[J]TSU2L[J,i,j]Hd[j]
	];
	
	\[CapitalOmega]1= Sqrt[2] g1[]RelabelIndices[
		y[q]Bar@q[a,i,p]**PR**Bt[]qt[a,i,p]-
		y[u]Bar@u[a,p]**PL**Bt[]ut[a,p]-
		y[d]Bar@d[a,p]**PL**Bt[]dt[a,p]+
		y[l]Bar@l[i,p]**PR**Bt[]lt[i,p]-
		y[e]Bar@e[p]**PL**Bt[]et[p]+
		y[Hu]Bar@Hut[i]**PR**Bt[]Hu[i]+
		y[Hd]Bar@Hdt[i]**PR**Bt[]Hd[i]
	];
	
	
	(* scalar gauge interactions *)
	\[CapitalLambda]3[A_]:=g3[]RelabelIndices[
		Bar@qt[a,i,p]TSU3c[A,a,c]qt[c,i,p]-
		Bar@ut[a,p]TSU3c[A,a,c]ut[c,p]-
		Bar@dt[a,p]TSU3c[A,a,c]dt[c,p]
	,
		Unique->True
	];
	
	\[CapitalLambda]2[I_]:=g2[]RelabelIndices[
		Bar@qt[a,i,p]TSU2L[I,i,j]qt[a,j,p]+
		Bar@lt[i,p]TSU2L[I,i,j]lt[j,p]+
		Bar@Hu[i]TSU2L[I,i,j]Hu[j]+
		Bar@Hd[i]TSU2L[I,i,j]Hd[j]
	,
		Unique->True
	];
	
	\[CapitalLambda]1[]:=g1[]RelabelIndices[
		y[q]Bar@qt[a,i,p]qt[a,i,p]-
		y[u]Bar@ut[a,p]ut[a,p]-
		y[d]Bar@dt[a,p]dt[a,p]+
		y[l]Bar@lt[i,p]lt[i,p]-
		y[e]Bar@et[p]et[p]+
		y[Hu]Bar@Hu[i]Hu[i]+
		y[Hd]Bar@Hd[i]Hd[i]
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


MSSM`SuperPotentialMSSM[] := Module[
	{
		\[CurlyEpsilon],
		\[ScriptCapitalL]Wscalar, \[ScriptCapitalL]Wfermion,
		(*index names*)
		a,c,
		i,j,
		p,r,s,t
	}
	,
	(* define SU2L \[Epsilon] tensor *)
	\[CurlyEpsilon][i_,j_] := CG[eps@SU2L, {i,j}];
	
	
	(* scalar contribution to superpotential *)
	\[ScriptCapitalL]Wscalar = -RelabelIndices@Expand[
     	(* scalar bilinears *)
     	Bar@\[Mu]t[]\[Mu]t[](Bar@Hu[i]Hu[i] + Bar@Hd[i]Hd[i]) -
      	(* scalar trilinears *)
      	PlusHc[Bar@\[Mu]t[]( 
         		Bar@Hd[i]Bar@ut[a,p]Bar@yu[r,p]qt[a,i,r] +
          		Bar@Hu[i]Bar@dt[a,p]Bar@yd[r,p]qt[a,i,r] +
          		Bar@Hu[i]Bar@et[p]Bar@ye[r,p]lt[i,r]
         	)] +
      	(* 2 sfermions and 2 Higgs *)
      	(Bar@Hu[i]Hu[i]) (Bar@ut[a,p]Bar@yu[r,p]yu[r,s]ut[a,s]) +
      	(Bar@Hd[i]Hd[i]) (Bar@dt[a,p]Bar@yd[r,p]yd[r,s]dt[a,s]) +
      	(Bar@Hd[i]Hd[i]) (Bar@et[p]Bar@ye[r,p]ye[r,s]et[s]) +
      	(Bar@Hd[i]Hd[i]) (Bar@lt[j,p]ye[p,r]Bar@ye[s,r]lt[j,s]) +
      	(Bar@Hd[i]Hd[i]) (Bar@qt[a,j,p]yd[p,r]Bar@yd[s,r]qt[a,j,s]) +
      	(Bar@Hu[i]Hu[i]) (Bar@qt[a,j,p]yu[p,r]Bar@yu[s,r]qt[a,j,s]) -
      	PlusHc[(Bar@Hd[i]Hu[i]) (Bar@ut[a,p]Bar@yu[r,p]yd[r,s]dt[a,s])] -
      	Bar@qt[a,i,p]Hu[i]yu[p,r]Bar@yu[s,r]Bar@Hu[j]qt[a,j,s] -
      	Bar@qt[a,i,p]Hd[i]yd[p,r]Bar@yd[s,r]Bar@Hd[j]qt[a,j,s] -
      	Bar@lt[i,p]Hd[i]ye[p,r]Bar@ye[s,r]Bar@Hd[j]lt[j,s] +
      	(* 4 sfermions *)
      	(Bar@ut[a,p]Bar@yu[r,p]qt[a,i,r]) Bar[Bar@ut[c,s]Bar@yu[t,s]qt[c,i,t]] +
      	(Bar@dt[a,p]Bar@yd[r,p]qt[a,i,r]) Bar[Bar@dt[c,s]Bar@yd[t,s]qt[c,i,t]] +
      	(Bar@et[p]Bar@ye[r,p]lt[i,r]) Bar[Bar@et[s]Bar@ye[t,s]lt[i,t]] +
      	PlusHc[(Bar@dt[a,p]Bar@yd[r,p]qt[a,i,r]) Bar[Bar@et[s]Bar@ye[t,s]lt[i,t]]]
     ];
     
     
     (* fermonic contribution to superpotential *)
     \[ScriptCapitalL]Wfermion=-PlusHc@RelabelIndices[
		(* Higgsino mass *)
		\[Mu]t[]Bar@\[CurlyEpsilon][i,j]Bar@CConj@Hut[i]**PL**Hdt[j]+
		(* SM Yukawas *)
		Bar@yu[r,p]Bar@u[a,p]**PL**q[a,i,r]Bar@\[CurlyEpsilon][i,j]Hu[j]-
		Bar@yd[r,p]Bar@d[a,p]**PL**q[a,i,r]Bar@\[CurlyEpsilon][i,j]Hd[j]-
		Bar@ye[r,p]Bar@e[p]**PL**l[i,r]Bar@\[CurlyEpsilon][i,j]Hd[j]+
		(* Remaining trilinears *)
		Bar@yu[r,p]Bar@u[a,p]**PL**Hut[j]qt[a,i,r]Bar@\[CurlyEpsilon][i,j]-
		Bar@yd[r,p]Bar@d[a,p]**PL**Hdt[j]qt[a,i,r]Bar@\[CurlyEpsilon][i,j]-
		Bar@ye[r,p]Bar@e[p]**PL**Hdt[j]lt[i,r]Bar@\[CurlyEpsilon][i,j]+
		(**)
		Bar@yu[r,p]Bar@ut[a,p]Bar@CConj@q[a,i,r]**PL**Hut[j]Bar@\[CurlyEpsilon][i,j]-
		Bar@yd[r,p]Bar@dt[a,p]Bar@CConj@q[a,i,r]**PL**Hdt[j]Bar@\[CurlyEpsilon][i,j]-
		Bar@ye[r,p]Bar@et[p]Bar@CConj@l[i,r]**PL**Hdt[j]Bar@\[CurlyEpsilon][i,j]
	];
	
	
	(* full super potential *)
	\[ScriptCapitalL]Wscalar+\[ScriptCapitalL]Wfermion
]


(* ::Subsection::Closed:: *)
(*Soft SUSY breaking terms*)


MSSM`SoftSUSYbreakingMSSM[] := Module[
	{
		\[CurlyEpsilon],
		\[ScriptCapitalL]SoftBreak,
		i,j,p,r,a
	}
	,
	(* define SU2L \[Epsilon] tensor *)
	\[CurlyEpsilon][i_,j_]:=CG[eps@SU2L, {i,j}];
	
	
	(* soft SUSY breaking *)
	\[ScriptCapitalL]SoftBreak= -RelabelIndices@PlusHc[
		b[]Bar@\[CurlyEpsilon][i,j]Hu[i]Hd[j] +
		Bar@qt[a,i,p]au[p,r]ut[a,r]\[CurlyEpsilon][i,j]Bar@Hu[j] -
		Bar@qt[a,i,p]ad[p,r]dt[a,r]\[CurlyEpsilon][i,j]Bar@Hd[j] -
		Bar@lt[i,p]ae[p,r]et[r]\[CurlyEpsilon][i,j]Bar@Hd[j]
	];
	
	\[ScriptCapitalL]SoftBreak
]


(* ::Section:: *)
(*Rotation of Higgs sector to the soft-SUSY mass basis*)


(* ::Text:: *)
(*With soft-SUSY mass basis we refer to the mass eigenstate basis after soft SUSY breaking, but before EWSB.*)


(* ::Subsection::Closed:: *)
(*Rotation to mass basis for the Higgs doublets [default]*)


MSSM`ToSoftSUSYMassBasis[L_] := Module[
	{
		Lag = L,
		To2HDM, ToHiggsDoubletMassBasis,
		mLight, mHeavy, mLightHeavy, mHeavyLight,
		\[CurlyEpsilon],
		j
	}
	,
	(* define anti-symmetric SUSubscript[(2), L] tensor *)
	\[CurlyEpsilon][i_,j_] := CG[eps@SU2L, {i,j}];
	
	(* rule to rotate from Hd to Hdc to obtain the usual 2HDM Lagrangian *)
	To2HDM = {Field[Hd, Scalar, {Index[i_,SU2L[fund]]}, devs_] :> RelabelIndices[-\[CurlyEpsilon][i,j]CD[devs,Bar@Hdc[j]], Unique->True]};
	
	(* rule to rotate to doublet mass basis *)
	ToHiggsDoubletMassBasis = {
		Field[Hu,  Scalar, {Index[i_,SU2L[fund]]}, devs_] :> (Sin[\[Gamma]] CD[devs,H[i]] + Cos[\[Gamma]] CD[devs,\[CapitalPhi][i]]),
		Field[Hdc, Scalar, {Index[i_,SU2L[fund]]}, devs_] :> (-Cos[\[Gamma]] CD[devs,H[i]] + Sin[\[Gamma]] CD[devs,\[CapitalPhi][i]])
	};
	
	(* apply rotations to doublet mass basis and simplify trigonometric functions *)
	Lag = CollectOperators[
		TrigReduce@ RelabelIndices@ Contract@ ContractCGs[Lag/.To2HDM/.ToHiggsDoubletMassBasis]
	];

	Lag = CollectOperators[
		TrigReduce/@Lag,
		Matchete`Simplifications`PackagePrivate`NormalForm->False
	];
	
	(* * * * * * * * * * *)

	(* extract current masses *)
	(* light mass *)
	mLight = CollectOperators[
		SelectOperatorClass[Lag,{Bar@H,H},0],
		Matchete`Simplifications`PackagePrivate`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	(*heavy mass*)
	mHeavy = CollectOperators[
		SelectOperatorClass[Lag,{Bar@\[CapitalPhi],\[CapitalPhi]},0],
		Matchete`Simplifications`PackagePrivate`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	(*mixed mass*)
	mLightHeavy = CollectOperators[
		HcSimplify@SelectOperatorClass[Lag,{Bar@\[CapitalPhi],H},0]/.HcTerms->Identity,
		Matchete`Simplifications`PackagePrivate`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	mHeavyLight = CollectOperators[
		HcSimplify@SelectOperatorClass[Lag,{Bar@\[CapitalPhi],H},0]/.HcTerms->Bar,
		Matchete`Simplifications`PackagePrivate`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	
	(* redefine all the Higgs masses *)
	Lag = RelabelIndices@Matchete`PackageScope`OperatorToNormalForm[
		Lag/.{
			mLight                  -> -mH2[],
			mHeavy                  -> -m\[CapitalPhi][]^2,
			mLightHeavy|mHeavyLight -> 0
		}
	];

	Lag
]


(* ::Subsection::Closed:: *)
(*Rotation to Higgs basis [ALTERNATIVE : not default]*)


MSSM`ToUnbrokenPhaseHiggsBasis[L_] := Module[
	{
		Lag = L,
		To2HDM, ToHiggsBasis,
		mLight, mHeavy, mLightHeavy, mHeavyLight,
		\[CurlyEpsilon],
		j
	}
	,
	(* define anti-symmetric SUSubscript[(2), L] tensor *)
	\[CurlyEpsilon][i_,j_] := CG[eps@SU2L, {i,j}];
	
	(* rule to rotate from Hd to Hdc to obtain the usual 2HDM Lagrangian *)
	To2HDM = {Field[Hd, Scalar, {Index[i_,SU2L[fund]]}, devs_] :> RelabelIndices[-\[CurlyEpsilon][i,j]CD[devs,Bar@Hdc[j]], Unique->True]};
	
	(* rule to rotate to Higgs basis *)
	ToHiggsBasis = {
		Field[Hu,  Scalar, {Index[i_,SU2L[fund]]}, devs_] :> (Sin[\[Beta]] CD[devs,H[i]] + Cos[\[Beta]] CD[devs,\[CapitalPhi][i]]),
		Field[Hdc, Scalar, {Index[i_,SU2L[fund]]}, devs_] :> (-Cos[\[Beta]] CD[devs,H[i]] + Sin[\[Beta]] CD[devs,\[CapitalPhi][i]])
	};
	
	(* apply rotations to Higgs basis*)
	Lag = CollectOperators[
		TrigReduce@ RelabelIndices@ Contract@ ContractCGs[Lag/.To2HDM/.ToHiggsBasis]
	];
	
	Lag = CollectOperators[
		TrigReduce/@Lag,
		Matchete`Simplifications`PackagePrivate`NormalForm->False
	];
	
	(* extract current masses *)
	(* light mass *)
	mLight = CollectOperators[
		SelectOperatorClass[Lag,{Bar@H,H},0],
		Matchete`Simplifications`PackagePrivate`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	(*heavy mass*)
	mHeavy = CollectOperators[
		SelectOperatorClass[Lag,{Bar@\[CapitalPhi],\[CapitalPhi]},0],
		Matchete`Simplifications`PackagePrivate`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	(*mixed mass*)
	mLightHeavy = CollectOperators[
		HcSimplify@SelectOperatorClass[Lag,{Bar@\[CapitalPhi],H},0]/.HcTerms->Identity,
		Matchete`Simplifications`PackagePrivate`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	mHeavyLight = CollectOperators[
		HcSimplify@SelectOperatorClass[Lag,{Bar@\[CapitalPhi],H},0]/.HcTerms->Bar,
		Matchete`Simplifications`PackagePrivate`NormalForm->False
	]/._Matchete`PackageScope`Operator->1;
	
	(* redefine all the Higgs masses *)
	Lag = RelabelIndices@Matchete`PackageScope`OperatorToNormalForm[
		Lag/.{
			mLight                  -> -mH2[],
			mHeavy                  -> -m\[CapitalPhi][]^2,
			mLightHeavy|mHeavyLight -> +mH2[]*Tan[2*\[Beta]]
		}
	];
	
	Lag
]


(* ::Subsection::Closed:: *)
(*Rotation to Higgsino mass basis*)


MSSM`ToHiggsinoMassBasis[L_] := Module[
	{
		Lag=L,
		RotateHiggsinos, CombineToVectorlikeFermion,
		\[CurlyEpsilon],
		j
	}
	,
	(* define anti-symmetric SUSubscript[(2), L] tensor *)
	\[CurlyEpsilon][i_,j_] := CG[eps@SU2L, {i,j}];
	
	(* replace Hut & Hdt by \[CapitalSigma] *)
	RotateHiggsinos={
		Field[Hut, Fermion, {Index[i_,SU2L[fund]]}, dev_List] :> PL**CD[dev,\[CapitalSigma][i]],
		Field[Hdt, Fermion, {Index[i_,SU2L[fund]]}, dev_List] :> RelabelIndices[-\[CurlyEpsilon][i,j]PL**CD[dev,CConj@\[CapitalSigma][j]], Unique->True]
	};
	
	(* combine chiral to vectorlike kinetic and mass terms *)
	CombineToVectorlikeFermion={
		(* kinetic terms *)
		Evaluate[
			I Bar[Field[\[CapitalSigma], Fermion, {i_}, {}]]**DiracProduct[GammaM[\[Mu]_]]**PL**Field[\[CapitalSigma], Fermion, {i_}, {\[Mu]_}] +
			I Bar[Field[\[CapitalSigma], Fermion, {j_}, {}]]**DiracProduct[GammaM[\[Nu]_]]**PR**Field[\[CapitalSigma], Fermion, {j_}, {\[Nu]_}]
		] :> I Bar[Field[\[CapitalSigma], Fermion, {i}, {}]]**DiracProduct[GammaM[\[Nu]]]**Field[\[CapitalSigma], Fermion, {i}, {\[Nu]}]
		,
		Evaluate[
			I Bar[Field[\[CapitalSigma], Fermion, {i_}, {}]]**DiracProduct[GammaM[\[Mu]_]]**PL**Field[\[CapitalSigma], Fermion, {i_}, {\[Mu]_}] -
			I Bar[Field[\[CapitalSigma], Fermion, {j_}, {\[Nu]_}]]**DiracProduct[GammaM[\[Nu]_]]**PR**Field[\[CapitalSigma], Fermion, {j_}, {}]
		] :> I Bar[Field[\[CapitalSigma], Fermion, {i}, {}]]**DiracProduct[GammaM[\[Nu]]]**Field[\[CapitalSigma], Fermion, {i}, {\[Nu]}]
		,
		Evaluate[
			-I Bar[Field[\[CapitalSigma], Fermion, {i_}, {\[Mu]_}]]**DiracProduct[GammaM[\[Mu]_]]**PL**Field[\[CapitalSigma], Fermion, {i_}, {}] +
			I Bar[Field[\[CapitalSigma], Fermion, {j_}, {}]]**DiracProduct[GammaM[\[Nu]_]]**PR**Field[\[CapitalSigma], Fermion, {j_}, {\[Nu]_}]
		] :> I Bar[Field[\[CapitalSigma], Fermion, {i}, {}]]**DiracProduct[GammaM[\[Nu]]]**Field[\[CapitalSigma], Fermion, {i}, {\[Nu]}]
		,
		(* mass terms *)
		Evaluate[
			Bar[Field[\[CapitalSigma], Fermion, {i_}, {}]]**DiracProduct[Proj[-1]]**Field[\[CapitalSigma], Fermion, {i_}, {}]
		] :> (Bar[Field[\[CapitalSigma], Fermion, {i}, {}]]**Field[\[CapitalSigma], Fermion, {i}, {}]-Bar[Field[\[CapitalSigma], Fermion, {i}, {}]]**DiracProduct[Proj[+1]]**Field[\[CapitalSigma], Fermion, {i}, {}])
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
(*It is more convenient to define several different couplings for sin(\[Beta]), sin(2\[Beta]), sin(4\[Beta]), cos(\[Beta]), cos(2\[Beta]), cos(4\[Beta]). *)
(*If only a couplings \[Beta] would be defined, the coefficients would become complicated functions of \[Beta].*)


MSSM`ReplaceRotationAngles[L_] := Module[
	{Lag = L, ReplaceTrig}
	,

	(* define replacement rules *)
	ReplaceTrig = {
		(* soft-SUSY mass basis *)
		Sin[\[Gamma]]  -> s\[Gamma][],
		Cos[\[Gamma]]  -> c\[Gamma][],
		Sin[2\[Gamma]] -> s2\[Gamma][],
		Cos[2\[Gamma]] -> c2\[Gamma][],
		Sin[4\[Gamma]] -> s4\[Gamma][],
		Cos[4\[Gamma]] -> c4\[Gamma][],
		(* Higgs basis *)
		Sin[\[Beta]]  -> s\[Beta][],
		Cos[\[Beta]]  -> c\[Beta][],
		Sin[2\[Beta]] -> s2\[Beta][],
		Cos[2\[Beta]] -> c2\[Beta][],
		Sin[4\[Beta]] -> s4\[Beta][],
		Cos[4\[Beta]] -> c4\[Beta][],
		Tan[2\[Beta]] -> t2\[Beta][]
	};
	
	(* apply replacement rules *)
	Lag = Contract@ContractCGs[Lag /.ReplaceTrig];
	
	Lag
]


(* ::Text:: *)
(*Manual simplification of trigonometric functions.*)


MSSM`TrigSimplifiy[expr_] := expr/.{
	s\[Beta][]*c\[Beta][] -> s2\[Beta][]/2, 
	s\[Gamma][]*c\[Gamma][] -> s2\[Gamma][]/2
}


(* ::Section:: *)
(*Cleaning up*)


(* ::Subsection::Closed:: *)
(*Remove temporarily defined fields and couplings*)


MSSM`RemoveTMPdef[L_] := (
	(* remove all fields defined intermediately *)
	RemoveField[Hu];
	RemoveField[Hd];
	RemoveField[Hdc];
	RemoveField[Hut];
	RemoveField[Hdt];
	RemoveCoupling[MHu];
	RemoveCoupling[MHd];
	
	(* correct the masses associated to fields *)
	Matchete`PackageScope`$FieldAssociation[\[CapitalPhi]][Mass]  = m\[CapitalPhi];
	RemoveCoupling[$M\[CapitalPhi]];
	Matchete`PackageScope`$FieldAssociation[\[CapitalSigma]][Mass]  = \[Mu]t;
	Matchete`PackageScope`$FieldAssociation[\[CapitalSigma]][Heavy] = True;
	
	(* remove all rotation angles that do not appear in the given Lagrangian *)
	Do[
		If[FreeQ[L,coup,All],
			RemoveCoupling@First[coup]
		]
		,
		{coup,{s\[Gamma][],c\[Gamma][],s2\[Gamma][],c2\[Gamma][],s4\[Gamma][],c4\[Gamma][],s\[Beta][],c\[Beta][],s2\[Beta][],c2\[Beta][],s4\[Beta][],c4\[Beta][],t2\[Beta][]}}
	];
)


(* ::Chapter:: *)
(*MSSM Lagrangian*)


(* ::Section:: *)
(*Lagrangian construction*)


(* ::Text:: *)
(*Derive full MSSM Lagrangian:*)


Module[{\[ScriptCapitalL]free, \[ScriptCapitalL]gauge, \[ScriptCapitalL]SuperPotential, \[ScriptCapitalL]SoftBreaking, \[ScriptCapitalL]MSSM},
	(* derive free MSSM Lagrangian *)
	\[ScriptCapitalL]free = MSSM`FreeLagrangianMSSM[];

	(* gauge interactions *)
	\[ScriptCapitalL]gauge = MSSM`GaugeLagrangianMSSM[];

	(* superpotential *)
	\[ScriptCapitalL]SuperPotential = MSSM`SuperPotentialMSSM[];

	(* soft SUSY-breaking terms *)
	\[ScriptCapitalL]SoftBreaking = MSSM`SoftSUSYbreakingMSSM[];

	(* full MSSM Lagrangian *)
	\[ScriptCapitalL]MSSM = \[ScriptCapitalL]free + \[ScriptCapitalL]gauge + \[ScriptCapitalL]SuperPotential + \[ScriptCapitalL]SoftBreaking;

	(*
	Rotate MSSM Lagrangian to a basis that is convenient for the matching.
	Use either of the two options implemented: 
		(a) soft-SUSY mass basis;
		(b) unbroken-phase Higgs basis.
	*)
	(* (a): soft-SUSY mass basis *)
	\[ScriptCapitalL]MSSM = MSSM`ToSoftSUSYMassBasis[\[ScriptCapitalL]MSSM];
	(* (b): unbroken-phase Higgs basis *)
	(*\[ScriptCapitalL]MSSM = MSSM`ToUnbrokenPhaseHiggsBasis[\[ScriptCapitalL]MSSM];*)

	(* rotate to Higgsino mass basis *)
	\[ScriptCapitalL]MSSM = MSSM`ToHiggsinoMassBasis[\[ScriptCapitalL]MSSM];
	
	(* replace rotation angles *)
	\[ScriptCapitalL]MSSM = MSSM`ReplaceRotationAngles[\[ScriptCapitalL]MSSM];

	(* simplify the Lagrangian by grouping everything by operator class *)
	\[ScriptCapitalL]MSSM = GreensSimplify[\[ScriptCapitalL]MSSM, ReductionIdentities->dDimensional];

	(* simplify some terms using trigonometric idnetities *)
	\[ScriptCapitalL]MSSM = MSSM`TrigSimplifiy[\[ScriptCapitalL]MSSM];

	(* remove intermediate auxiliary fields and couplings and properly set the field properties *)
	MSSM`RemoveTMPdef[\[ScriptCapitalL]MSSM]; (* without this line there are additional power-type supertraces computed, which do not contribute to the final result *)

	(* delete internal simplification identities for operators containing superpartners *)
	Matchete`PackageScope`ResetOperatorAssociations[];

	\[ScriptCapitalL]MSSM
]
