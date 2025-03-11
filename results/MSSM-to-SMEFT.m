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


(* ::Subsubsection:: *)
(*IntroduceEffectiveCouplings *)


(*
\[ScriptCapitalL]MassFree = \[ScriptCapitalL]MSSM/.{Coupling[mq|mu|md|ml|me|m\[CapitalPhi]|\[Mu]H|mB|mW|mG,___]->0};
\[ScriptCapitalL]Mass     = \[ScriptCapitalL]MSSM-\[ScriptCapitalL]MassFree;
\[ScriptCapitalL]MSSM     = \[ScriptCapitalL]Mass+IntroduceEffectiveCouplings[\[ScriptCapitalL]MassFree, EffectiveCouplingSymbol->"\[Kappa]"];
*)


(* ::Subsection::Closed:: *)
(*Restricting to sUp, sCharm, sTop*)


(*
\[ScriptCapitalL]MSSM = \[ScriptCapitalL]MSSM/.{Field[Alternatives@@{(*\[CapitalSigma],*)\[CapitalPhi],lt,et,qt,ut,dt,Gt,Wt(*,Bt*)},___]->0};
*)


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


(* ::Section::Closed:: *)
(*Manual evanescent operators [NOT IN USE]*)


(* ::Text:: *)
(*In this section we treat the evanescent operators that are generated when applying Fierz identities (which are intrinsically four dimensional) to the tree-level EFT Lagrangian in a dimensional regularization calculation.*)


(* ::Text:: *)
(*Our results for the one-loop EFT Lagrangian are derived in a renormalization scheme that is an evanescent-free version of the \!\(\*OverscriptBox[\(MS\), \(_\)]\) scheme, as described in [2211.09144].*)


(* ::Subsection::Closed:: *)
(*Definition of  relevant Warsaw basis operators*)


(* ::Text:: *)
(*Here we define short hands for the operators of the Warsaw that receive a finite renormalization in the evanescent-free scheme.*)


(* ::Subsubsection::Closed:: *)
(*Some short hands for generators, ...*)


\[Tau]SU2L[Jadj_,ifund_,jfund_]:= 2 CG[gen[SU2L[fund]],{Jadj,ifund,jfund}];

\[Epsilon]SU2L[ifund_,jfund_]:= CG[eps[SU2L],{ifund,jfund}];

tSU3c[Aadj_,\[Alpha]fund_,\[Beta]fund_]:=CG[gen[SU3c[fund]],{Aadj,\[Alpha]fund,\[Beta]fund}];


(* ::Subsubsection::Closed:: *)
(*Warsaw basis operators*)


Qle[p_,r_,s_,t_]:=Module[{i,\[Mu]},
	RelabelIndices[
		(Bar@l[i,p]**\[Gamma][\[Mu]]**l[i,r])(Bar@e[s]**\[Gamma][\[Mu]]**e[t])
		,
		Unique->True
	]
]

Qqu1[p_,r_,s_,t_]:=Module[{a,b,i,\[Mu]},
	RelabelIndices[
		(Bar@q[a,i,p]**\[Gamma][\[Mu]]**q[a,i,r])(Bar@u[b,s]**\[Gamma][\[Mu]]**u[b,t])
		,
		Unique->True
	]
]

Qqd1[p_,r_,s_,t_]:=Module[{a,b,i,\[Mu]},
	RelabelIndices[
		(Bar@q[a,i,p]**\[Gamma][\[Mu]]**q[a,i,r])(Bar@d[b,s]**\[Gamma][\[Mu]]**d[b,t])
		,
		Unique->True
	]
]

Qqu8[p_,r_,s_,t_]:=Module[{a,b,c,d,A,i,\[Mu]},
	RelabelIndices[
		tSU3c[A,a,b]tSU3c[A,c,d](Bar@q[a,i,p]**\[Gamma][\[Mu]]**q[b,i,r])(Bar@u[c,s]**\[Gamma][\[Mu]]**u[d,t])
		,
		Unique->True
	]
]

Qqd8[p_,r_,s_,t_]:=Module[{a,b,c,e,A,i,\[Mu]},
	RelabelIndices[
		tSU3c[A,a,b]tSU3c[A,c,e](Bar@q[a,i,p]**\[Gamma][\[Mu]]**q[b,i,r])(Bar@d[c,s]**\[Gamma][\[Mu]]**d[e,t])
		,
		Unique->True
	]
]

Qee[p_,r_,s_,t_]:=Module[{\[Mu]},
	RelabelIndices[
		(Bar@e[p]**\[Gamma][\[Mu]]**e[r])(Bar@e[s]**\[Gamma][\[Mu]]**e[t])
		,
		Unique->True
	]
]

Quu[p_,r_,s_,t_]:=Module[{a,b,\[Mu]},
	RelabelIndices[
		(Bar@u[a,p]**\[Gamma][\[Mu]]**u[a,r])(Bar@u[b,s]**\[Gamma][\[Mu]]**u[b,t])
		,
		Unique->True
	]
]

Qdd[p_,r_,s_,t_]:=Module[{a,b,\[Mu]},
	RelabelIndices[
		(Bar@d[a,p]**\[Gamma][\[Mu]]**d[a,r])(Bar@d[b,s]**\[Gamma][\[Mu]]**d[b,t])
		,
		Unique->True
	]
]

Qud1[p_,r_,s_,t_]:=Module[{a,b,\[Mu]},
	RelabelIndices[
		(Bar@u[a,p]**\[Gamma][\[Mu]]**u[a,r])(Bar@d[b,s]**\[Gamma][\[Mu]]**d[b,t])
		,
		Unique->True
	]
]

Qud8[p_,r_,s_,t_]:=Module[{a,b,c,f,A,\[Mu]},
	RelabelIndices[
		tSU3c[A,a,b]tSU3c[A,c,f](Bar@u[a,p]**\[Gamma][\[Mu]]**u[b,r])(Bar@d[c,s]**\[Gamma][\[Mu]]**d[f,t])
		,
		Unique->True
	]
]

Qll[p_,r_,s_,t_]:=Module[{i,j,\[Mu]},
	RelabelIndices[
		(Bar@l[i,p]**\[Gamma][\[Mu]]**l[i,r])(Bar@l[j,s]**\[Gamma][\[Mu]]**l[j,t])
		,
		Unique->True
	]
]

Qqq1[p_,r_,s_,t_]:=Module[{a,b,i,j,\[Mu]},
	RelabelIndices[
		(Bar@q[a,i,p]**\[Gamma][\[Mu]]**q[a,i,r])(Bar@q[b,j,s]**\[Gamma][\[Mu]]**q[b,j,t])
		,
		Unique->True
	]
]

Qqq3[p_,r_,s_,t_]:=Module[{a,b,i,j,k,l,J,\[Mu]},
	RelabelIndices[
		\[Tau]SU2L[J,i,j]\[Tau]SU2L[J,k,l](Bar@q[a,i,p]**\[Gamma][\[Mu]]**q[a,j,r])(Bar@q[b,k,s]**\[Gamma][\[Mu]]**q[b,l,t])
		,
		Unique->True
	]
]

Qquqd1[p_,r_,s_,t_]:=Module[{i,j,a,b},
	RelabelIndices[
		\[Epsilon]SU2L[i,j](Bar@q[a,i,p]**u[a,r])(Bar@q[b,j,s]**d[b,t])
		,
		Unique->True
	]
]

Qquqd8[p_,r_,s_,t_]:=Module[{i,j,a,b,c,f,A},
	RelabelIndices[
		tSU3c[A,a,b]tSU3c[A,c,f]\[Epsilon]SU2L[i,j](Bar@q[a,i,p]**u[b,r])(Bar@q[c,j,s]**d[f,t])
		,
		Unique->True
	]
]

Qlequ1[p_,r_,s_,t_]:=Module[{i,j,a},
	RelabelIndices[
		\[Epsilon]SU2L[i,j](Bar@l[i,p]**e[r])(Bar@q[a,j,s]**u[a,t])
		,
		Unique->True
	]
]

Qledq[p_,r_,s_,t_]:=Module[{i,a},
	RelabelIndices[
		(Bar@l[i,p]**e[r])(Bar@d[a,s]**q[a,i,t])
		,
		Unique->True
	]
]

QeB[p_,r_]:=Module[{i,\[Mu],\[Nu]},
	RelabelIndices[
		(Bar@l[i,p]**\[Sigma][\[Mu],\[Nu]]**e[r]) H[i]FS[B,\[Mu],\[Nu]]
		,
		Unique->True
	]
]

QuB[p_,r_]:=Module[{a,i,j,\[Mu],\[Nu]},
	RelabelIndices[
		(Bar@q[a,i,p]**\[Sigma][\[Mu],\[Nu]]**u[a,r]) \[Epsilon]SU2L[i,j]Bar@H[j]FS[B,\[Mu],\[Nu]]
		,
		Unique->True
	]
]

QdB[p_,r_]:=Module[{a,i,\[Mu],\[Nu]},
	RelabelIndices[
		(Bar@q[a,i,p]**\[Sigma][\[Mu],\[Nu]]**d[a,r]) H[i]FS[B,\[Mu],\[Nu]]
		,
		Unique->True
	]
]

QeW[p_,r_]:=Module[{i,j,J,\[Mu],\[Nu]},
	RelabelIndices[
		(Bar@l[i,p]**\[Sigma][\[Mu],\[Nu]]**e[r])\[Tau]SU2L[J,i,j] H[j]FS[W,\[Mu],\[Nu],J]
		,
		Unique->True
	]
]

QuW[p_,r_]:=Module[{a,i,j,k,J,\[Mu],\[Nu]},
	RelabelIndices[
		(Bar@q[a,i,p]**\[Sigma][\[Mu],\[Nu]]**u[a,r])\[Tau]SU2L[J,i,j] \[Epsilon]SU2L[j,k]Bar@H[k]FS[W,\[Mu],\[Nu],J]
		,
		Unique->True
	]
]

QdW[p_,r_]:=Module[{a,i,j,J,\[Mu],\[Nu]},
	RelabelIndices[
		(Bar@q[a,i,p]**\[Sigma][\[Mu],\[Nu]]**d[a,r])\[Tau]SU2L[J,i,j] H[j]FS[W,\[Mu],\[Nu],J]
		,
		Unique->True
	]
]

QeH[p_,r_]:=Module[{i,j},
	RelabelIndices[
		(Bar@H[i]H[i])(Bar@l[j,p]**e[r])H[j]
		,
		Unique->True
	]
]

QuH[p_,r_]:=Module[{a,i,j,k},
	RelabelIndices[
		(Bar@H[i]H[i])(Bar@q[a,j,p]**u[a,r])\[Epsilon]SU2L[j,k]Bar@H[k]
		,
		Unique->True
	]
]

QdH[p_,r_]:=Module[{a,i,j},
	RelabelIndices[
		(Bar@H[i]H[i])(Bar@q[a,j,p]**d[a,r])H[j]
		,
		Unique->True
	]
]

Qye[p_,r_]:=Module[{i},
	RelabelIndices[
		(Bar@l[i,p]**e[r])H[i]
		,
		Unique->True
	]
]

Qyu[p_,r_]:=Module[{a,i,j},
	RelabelIndices[
		(Bar@q[a,i,p]**u[a,r])\[Epsilon]SU2L[i,j]Bar@H[j]
		,
		Unique->True
	]
]

Qyd[p_,r_]:=Module[{a,i},
	RelabelIndices[
		(Bar@q[a,i,p]**d[a,r])H[i]
		,
		Unique->True
	]
]


(* ::Subsection::Closed:: *)
(*Evanescent shifts (finite renormalizations)*)


(* ::Text:: *)
(*These are the shifts corresponding to the finite renormalizations in the evanescent-free scheme.*)
(*These results are extracted from the ancillary material of [2211.09144].*)


EvanescentShift= {
((Bar@l[i_,p_]**e[r_])(Bar@e[s_]**l[i_,t_])):>RelabelIndices[-(1/2)Qle[p,t,s,r]+hbar(1/4 Bar@ye[u,v]ye[t,s]Qle[p,u,v,r]+1/4 Bar@ye[p,r]ye[u,v]Qle[u,t,s,v]+3/8 g1[]Bar@ye[p,r]Bar@QeB[t,s]+3/8 g1[]ye[t,s]QeB[p,r]+1/2 Bar@ye[p,r]Bar@yu[u,v]Bar@Qlequ1[t,s,u,v]+1/2 ye[t,s]yu[u,v]Qlequ1[p,r,u,v]+Bar@ye[p,u]Bar@ye[v,r]ye[v,u]Bar@QeH[t,s]+QeH[p,r](Bar@ye[u,v]ye[t,v]ye[u,s]-1/2 \[Lambda][]ye[t,s])-1/8 g2[]Bar@ye[p,r]Bar@QeW[t,s]-1/8 g2[]ye[t,s]QeW[p,r]-1/4 Bar@ye[u,r]ye[t,v]Qle[p,u,s,v]-1/4 Bar@ye[p,u]ye[v,s]Qle[v,t,u,r]-1/4 Bar@ye[u,r]ye[v,s]Qll[v,u,p,t]-1/2 \[Lambda][]Bar@ye[p,r]Bar@QeH[t,s]-1/2 \[Mu]2[] Bar@ye[p,r]Bar@Qye[t,s]-1/2 Bar@ye[p,r]yd[v,u]Bar@Qledq[t,s,u,v]-1/2 Bar@ye[p,u]ye[t,v]Qee[u,r,s,v]-1/2 Bar@yd[u,v]ye[t,s]Qledq[p,r,v,u]-1/2 \[Mu]2[] ye[t,s]Qye[p,r]),Unique->True]
,
((Bar@q[a_,i_,p_]**u[a_,r_])(Bar@u[b_,s_]**q[b_,i_,t_])):>RelabelIndices[-(1/6)Qqu1[p,t,s,r]-Qqu8[p,t,s,r]+hbar(1/12 yd[t,u]yu[v,s]Qquqd1[v,r,p,u]+1/4 Bar@yu[u,v]yu[t,s]Qqu1[p,u,v,r]+1/4 Bar@yu[p,r]yu[u,v]Qqu1[u,t,s,v]+1/2 yd[t,u]yu[v,s]Qquqd8[v,r,p,u]+Bar@yd[p,u]Bar@yu[v,r](1/12 Bar@Qquqd1[v,s,t,u]+1/2 Bar@Qquqd8[v,s,t,u]-1/2 Bar@Qquqd1[t,s,v,u])+QuH[p,r](3Bar@yu[u,v]yu[t,v]yu[u,s]-3/2 \[Lambda][]yu[t,s])+3/2 Bar@ye[u,v]Bar@yu[p,r]Bar@Qlequ1[u,v,t,s]+3/2 ye[u,v]yu[t,s]Qlequ1[u,v,p,r]+3/2 Bar@yu[u,v]yu[t,s]Qqu8[p,u,v,r]+3/2 Bar@yu[p,r]yu[u,v]Qqu8[u,t,s,v]+3Bar@yu[p,u]Bar@yu[v,r]yu[v,u]Bar@QuH[t,s]-1/8 Bar@yu[u,r]yu[v,s]Qqq1[v,t,p,u]-1/8 Bar@yu[u,r]yu[v,s]Qqq3[v,t,p,u]-1/6 Bar@yd[p,u]yd[t,v]Qud1[s,r,u,v]-1/4 Bar@yu[u,r]yu[t,v]Qqu1[p,u,s,v]-1/4 Bar@yu[p,u]yu[v,s]Qqu1[v,t,u,r]-3/8 g2[]Bar@yu[p,r]Bar@QuW[t,s]-3/8 g2[]yu[t,s]QuW[p,r]-1/2 yd[t,u]yu[v,s]Qquqd1[p,r,v,u]-1/2 Bar@yu[p,u]yu[t,v]Quu[u,r,s,v]-5/8 g1[]Bar@yu[p,r]Bar@QuB[t,s]-5/8 g1[]yu[t,s]QuB[p,r]-Bar@yd[p,u]yd[t,v]Qud8[s,r,u,v]-3/2 Bar@yd[u,v]Bar@yu[p,r]Bar@Qquqd1[t,s,u,v]-3/2 \[Lambda][]Bar@yu[p,r]Bar@QuH[t,s]-3/2 \[Mu]2[] Bar@yu[p,r]Bar@Qyu[t,s]-3/2 yd[u,v]yu[t,s]Qquqd1[p,r,u,v]-3/2 \[Mu]2[] yu[t,s]Qyu[p,r]),Unique->True]
,
((Bar@q[a_,i_,p_]**d[a_,r_])(Bar@d[b_,s_]**q[b_,i_,t_])):>RelabelIndices[-(1/6)Qqd1[p,t,s,r]-Qqd8[p,t,s,r]+hbar(1/12 yd[u,s]yu[t,v]Qquqd1[p,v,u,r]+1/8 g1[]Bar@yd[p,r]Bar@QdB[t,s]+1/8 g1[]yd[t,s]QdB[p,r]+1/4 Bar@yd[u,v]yd[t,s]Qqd1[p,u,v,r]+1/4 Bar@yd[p,r]yd[u,v]Qqd1[u,t,s,v]+1/2 yd[u,s]yu[t,v]Qquqd8[p,v,u,r]+Bar@yd[u,r]Bar@yu[p,v](1/12 Bar@Qquqd1[t,v,u,s]+1/2 Bar@Qquqd8[t,v,u,s]-1/2 Bar@Qquqd1[u,v,t,s])+QdH[p,r](3Bar@yd[u,v]yd[t,v]yd[u,s]-3/2 \[Lambda][]yd[t,s])+Qquqd1[u,v,p,r](-(1/2)yd[u,s]yu[t,v]-3/2 yd[t,s]yu[u,v])+3/2 Bar@yd[u,v]yd[t,s]Qqd8[p,u,v,r]+3/2 Bar@yd[p,r]yd[u,v]Qqd8[u,t,s,v]+3Bar@yd[p,u]Bar@yd[v,r]yd[v,u]Bar@QdH[t,s]-1/8 Bar@yd[u,r]yd[v,s]Qqq1[v,t,p,u]-1/8 Bar@yd[u,r]yd[v,s]Qqq3[v,t,p,u]-1/6 Bar@yu[p,u]yu[t,v]Qud1[u,v,s,r]-1/4 Bar@yd[u,r]yd[t,v]Qqd1[p,u,s,v]-1/4 Bar@yd[p,u]yd[v,s]Qqd1[v,t,u,r]-3/8 g2[]Bar@yd[p,r]Bar@QdW[t,s]-3/8 g2[]yd[t,s]QdW[p,r]-1/2 Bar@yd[p,u]yd[t,v]Qdd[u,r,s,v]-Bar@yu[p,u]yu[t,v]Qud8[u,v,s,r]-3/2 Bar@yd[p,r]Bar@yu[u,v]Bar@Qquqd1[u,v,t,s]-3/2 \[Lambda][]Bar@yd[p,r]Bar@QdH[t,s]-3/2 \[Mu]2[] Bar@yd[p,r]Bar@Qyd[t,s]-3/2 Bar@ye[u,v]yd[t,s]Bar@Qledq[u,v,r,p]-3/2 Bar@yd[p,r]ye[u,v]Qledq[u,v,s,t]-3/2 \[Mu]2[] yd[t,s]Qyd[p,r]),Unique->True]
};


(* ::Subsection::Closed:: *)
(*Applying evanescent shifts to tree-level Lagrangian to change to the evanescent-free scheme*)


(* the shifts are in terms of Warsaw absis couplings and not the couplings in the matched Lagrangian *)
(*\[ScriptCapitalL]EFT$TreeLevel= GreensSimplify[
	RelabelIndices[Expand[\[ScriptCapitalL]EFT0]/.EvanescentShift] 
	,
	ReductionIdentities->dDimensional
];*)


(* ::Section::Closed:: *)
(*Automatic evanescent operators*)


\[ScriptCapitalL]EFT$TreeLevel= GreensSimplify[Expand[\[ScriptCapitalL]EFT0], ReductionIdentities -> EvanescenceFree]


(* ::Section::Closed:: *)
(*One-loop matching the full MSSM*)


(* ::Text:: *)
(*Perform the one-loop matching, excluding the tree-level results that were computed above.*)


\[ScriptCapitalL]EFT1= EchoTiming[CollectOperators@ Match[\[ScriptCapitalL]MSSM, EFTOrder->6, LoopOrder->{1} (* this excludes the tree-level results *)],"One-loop matching"];


(* ::Section:: *)
(*Simplifying the EFT Lagrangian*)


(* ::Subsection::Closed:: *)
(*Off-Shell*)


(* ::Text:: *)
(*For the one-loop results we can simply apply the four-dimensional Fierz identities, since the evanescent operators generated by this are of two-loop order and can be neglected here.*)


\[ScriptCapitalL]EFT$OneLoop= EchoTiming[GreensSimplify[\[ScriptCapitalL]EFT1,ReductionIdentities->FourDimensional],"Off-shell simplifications"];


(* ::Text:: *)
(*Combine tree-level and one-loop results, giving the full EFT Lagrangian.*)


\[ScriptCapitalL]EFT$OffShell= \[ScriptCapitalL]EFT$TreeLevel + \[ScriptCapitalL]EFT$OneLoop;


(* ::Text:: *)
(*Renormalize in evanescent-free MSbar scheme:*)


(* why is \[ScriptD] still appearing here? *)
\[ScriptCapitalL]EFT$OffShell= Matchete`PackageScope`LagrangianExpand[\[ScriptCapitalL]EFT$OffShell/.(\[ScriptD]->(4-2\[Epsilon]))]/.(\[Epsilon]^-1->0)/.(\[Epsilon]->0);


(* ::Text:: *)
(*Export the resulting EFT Lagrangian*)


Export[FileNameJoin[{"results","MSSM-EFT-Lagrangian_off-shell.m"}], \[ScriptCapitalL]EFT$OffShell];


(* ::Subsection::Closed:: *)
(*On-Shell*)


(* ::Text:: *)
(*Apply the necessary field redefinitions to map to an on-shell basis.*)


\[ScriptCapitalL]EFT$OnShell= EchoTiming[EOMSimplify[
	\[ScriptCapitalL]EFT$OffShell
	,
	DummyCoefficients->True,
	ReductionIdentities->dDimensional
],"Field redefinitions"];


(* ::Text:: *)
(*The field redefinitions can introduce operators that require Fierz identities to be mapped onto the Warsaw basis.*)


\[ScriptCapitalL]EFT$OnShell= EchoTiming[GreensSimplify[
	\[ScriptCapitalL]EFT$OnShell
	,
	ReductionIdentities->EvanescenceFree
],"Fierzing the final results"];


\[ScriptCapitalL]EFT$OnShell$replaced=GreensSimplify[
	EchoTiming[ReplaceEffectiveCouplings@\[ScriptCapitalL]EFT$OnShell, "ReplaceEffectiveCouplings"]
	,
	ReductionIdentities->dDimensional
];


(* ::Text:: *)
(*Export the resulting EFT Lagrangian*)


Export[FileNameJoin[{"results","MSSM-EFT-Lagrangian_on-shell.m"}], \[ScriptCapitalL]EFT$OnShell$replaced];
Export[FileNameJoin[{"results","MSSM-EFT-Lagrangian_on-shell-eff.m"}], \[ScriptCapitalL]EFT$OnShell];


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
		ReplaceEffectiveCouplings[\[ScriptCapitalL]EFT$OnShell$simpl],
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
