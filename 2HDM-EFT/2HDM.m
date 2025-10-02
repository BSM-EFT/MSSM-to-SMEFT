(* ::Package:: *)

(* ::Title:: *)
(*2HDM model file*)


(* ::Section:: *)
(*Gauge groups and flavors*)


(* ::Text:: *)
(*Gauge groups*)


DefineGaugeGroup[SU3c, SU@3, gs, G, 
	FundAlphabet-> CharacterRange["a", "f"], 
	AdjAlphabet -> CharacterRange["A", "F"], 
	NiceForm->{"\!\(\*SubscriptBox[\(g\), \(s\)]\)", Default}
];
DefineGaugeGroup[SU2L, SU@2, gL, W, 
	FundAlphabet -> CharacterRange["i", "n"], 
	AdjAlphabet -> CharacterRange["I", "N"], 
	NiceForm->{"\!\(\*SubscriptBox[\(g\), \(L\)]\)", Default}
];
DefineGaugeGroup[U1Y, U1, gY, B, 
	NiceForm->{"\!\(\*SubscriptBox[\(g\), \(Y\)]\)", Default}
];


(* ::Text:: *)
(*Flavor index*)


ParameterDefault[Nf-> 3]


DefineFlavorIndex[Flavor, Nf,
	IndexAlphabet-> {"p", "r", "s", "t", "u", "v"}
];


(* ::Section:: *)
(*Fields*)


(* ::Text:: *)
(*Matter fields*)


DefineField[q, Fermion, 
	Indices-> {SU3c[fund], SU2L[fund], Flavor}, 
	Charges->{U1Y[1/6]}, 
	Chiral-> LeftHanded, 
	Mass-> 0
];
DefineField[u, Fermion,
	Indices-> {SU3c[fund], Flavor},
	Charges->{U1Y[2/3]}, 
	Chiral-> RightHanded, 
	Mass-> 0
];
DefineField[d, Fermion,
	Indices-> {SU3c[fund], Flavor},
	Charges->{U1Y[-1/3]}, 
	Chiral-> RightHanded, 
	Mass-> 0
];
DefineField[l, Fermion,
	Indices-> {SU2L[fund], Flavor},
	Charges->{U1Y[-1/2]}, 
	Chiral-> LeftHanded, 
	Mass-> 0
];
DefineField[e, Fermion,
	Indices-> {Flavor},
	Charges->{U1Y[-1]}, 
	Chiral-> RightHanded, 
	Mass-> 0
];


(* ::Text:: *)
(*Higgs fields*)


DefineField[\[Phi]1, Scalar, 
	Indices-> {SU2L[fund]}, 
	Charges-> {U1Y[1/2]}, 
	Mass-> 0,
	NiceForm-> "\!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\)"
];
DefineField[\[Phi]2, Scalar, 
	Indices-> {SU2L[fund]}, 
	Charges-> {U1Y[1/2]}, 
	Mass-> 0,
	NiceForm-> "\!\(\*SubscriptBox[\(\[Phi]\), \(2\)]\)"
];


SubsuperscriptBox[a,b,c]//DisplayForm


(* ::Section:: *)
(*Couplings*)


(* ::Text:: *)
(*Yukawa couplings*)


DefineCoupling[Y1u, 
	Indices-> {Flavor, Flavor},
	NiceForm-> "\!\(\*SubsuperscriptBox[\(Y\), \(u\), \(1\)]\)"
];
DefineCoupling[Y1d, 
	Indices-> {Flavor, Flavor},
	NiceForm-> "\!\(\*SubsuperscriptBox[\(Y\), \(d\), \(1\)]\)"
];
DefineCoupling[Y1e, 
	Indices-> {Flavor, Flavor},
	NiceForm-> "\!\(\*SubsuperscriptBox[\(Y\), \(e\), \(1\)]\)"
];
DefineCoupling[Y2u, 
	Indices-> {Flavor, Flavor},
	NiceForm-> "\!\(\*SubsuperscriptBox[\(Y\), \(u\), \(2\)]\)"
];
DefineCoupling[Y2d, 
	Indices-> {Flavor, Flavor},
	NiceForm-> "\!\(\*SubsuperscriptBox[\(Y\), \(d\), \(2\)]\)"
];
DefineCoupling[Y2e, 
	Indices-> {Flavor, Flavor},
	NiceForm-> "\!\(\*SubsuperscriptBox[\(Y\), \(e\), \(2\)]\)"
];


(* ::Text:: *)
(*Higgs mass terms*)


DefineCoupling[m12, 
	SelfConjugate-> True, 
	EFTOrder-> 2,
	NiceForm-> "\!\(\*TemplateBox[{\"m\", \"1\", \"2\"},\n\"Subsuperscript\"]\)"
];
DefineCoupling[m22, 
	SelfConjugate-> True, 
	EFTOrder-> 2,
	NiceForm-> "\!\(\*TemplateBox[{\"m\", \"2\", \"2\"},\n\"Subsuperscript\"]\)"
];
DefineCoupling[m122, 
	SelfConjugate-> False, 
	EFTOrder-> 2,
	NiceForm-> "\!\(\*TemplateBox[{\"m\", \"12\", \"2\"},\n\"Subsuperscript\"]\)"
];


(* ::Text:: *)
(*Higgs quartic couplings*)


DefineCoupling[\[Lambda]1, 
	SelfConjugate-> True, 
	EFTOrder-> 0,
	NiceForm-> "\!\(\*SubscriptBox[\(\[Lambda]\), \(1\)]\)"
];
DefineCoupling[\[Lambda]2, 
	SelfConjugate-> True, 
	EFTOrder-> 0,
	NiceForm-> "\!\(\*SubscriptBox[\(\[Lambda]\), \(2\)]\)"
];
DefineCoupling[\[Lambda]3, 
	SelfConjugate-> True, 
	EFTOrder-> 0,
	NiceForm-> "\!\(\*SubscriptBox[\(\[Lambda]\), \(3\)]\)"
];
DefineCoupling[\[Lambda]4, 
	SelfConjugate-> True, 
	EFTOrder-> 0,
	NiceForm-> "\!\(\*SubscriptBox[\(\[Lambda]\), \(4\)]\)"
];
DefineCoupling[\[Lambda]5, 
	SelfConjugate-> False, 
	EFTOrder-> 0,
	NiceForm-> "\!\(\*SubscriptBox[\(\[Lambda]\), \(5\)]\)"
];
DefineCoupling[\[Lambda]6, 
	SelfConjugate-> False, 
	EFTOrder-> 0,
	NiceForm-> "\!\(\*SubscriptBox[\(\[Lambda]\), \(6\)]\)"
];
DefineCoupling[\[Lambda]7, 
	SelfConjugate-> False, 
	EFTOrder-> 0,
	NiceForm-> "\!\(\*SubscriptBox[\(\[Lambda]\), \(6\)]\)"
];


(* ::Section:: *)
(*Lagrangian*)


Module[{p,r,s,i,j,\[Alpha],\[Beta],YukawaL,HiggsPotential},
	
	YukawaL =  Y1u[p,r] Bar@q[\[Alpha],i,p]**u[\[Alpha],r] CG[eps[SU2L],{i,j}]Bar@\[Phi]1[j] + 
		Y1d[p,r] Bar@q[\[Alpha],i,p]**d[\[Alpha],r]\[Phi]1[i] + 
		Y1e[p,r] Bar@l[i,p]**e[r]\[Phi]1[i] +
		Y2u[p,r] Bar@q[\[Alpha],i,p]**u[\[Alpha],r] CG[eps[SU2L],{i,j}]Bar@\[Phi]2[j] + 
		Y2d[p,r] Bar@q[\[Alpha],i,p]**d[\[Alpha],r]\[Phi]2[i] + 
		Y2e[p,r] Bar@l[i,p]**e[r]\[Phi]2[i];
		
	HiggsPotential = m12[]Bar@\[Phi]1[i]\[Phi]1[i] + 
		m22[]Bar@\[Phi]2[i]\[Phi]2[i] + 
		PlusHc[m122[]Bar@\[Phi]1[i]\[Phi]2[i]] +
		\[Lambda]1[]/2 (Bar@\[Phi]1[i]\[Phi]1[i])(Bar@\[Phi]1[j]\[Phi]1[j]) +
		\[Lambda]2[]/2 (Bar@\[Phi]2[i]\[Phi]2[i])(Bar@\[Phi]2[j]\[Phi]2[j]) +
		\[Lambda]3[] (Bar@\[Phi]1[i]\[Phi]1[i])(Bar@\[Phi]2[j]\[Phi]2[j]) +
		\[Lambda]4[] (Bar@\[Phi]1[i]\[Phi]2[i])(Bar@\[Phi]2[j]\[Phi]1[j]) +
		PlusHc[
			\[Lambda]5[]/2 (Bar@\[Phi]1[i]\[Phi]1[i])(Bar@\[Phi]1[j]\[Phi]1[j]) +
			\[Lambda]6[] (Bar@\[Phi]1[i]\[Phi]1[i])(Bar@\[Phi]1[j]\[Phi]2[j]) +
			\[Lambda]7[] (Bar@\[Phi]2[i]\[Phi]2[i])(Bar@\[Phi]1[j]\[Phi]2[j])
		];
	
	FreeLag[q, u, d, l, e, \[Phi]1, \[Phi]2, G, W, B] - PlusHc[YukawaL] - HiggsPotential //RelabelIndices
]
