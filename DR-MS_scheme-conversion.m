(* ::Package:: *)

(* ::Title:: *)
(*DRbar\[LongDash]MSbar scheme conversion*)


Needs["Matchete`"]


Begin["SchemeChange`"]


(* ::Section::Closed:: *)
(*Extracting the Lagrangian terms relevant for tree-level matching*)


GetTreeLevelTerms::usage="";

GetTreeLevelTerms::sanitycheck="Could not correctly identify tree-level terms."

GetTreeLevelTerms[lag_]:=Module[{L,Ltree1,Ltree2,Ltree,Lloop},
	(* only keep tree-level contributions *)
	L = lag/. hbar->0;
	
	(* keep only terms with at most one heavy field *)
	L = Expand[L];
	
	(* get all Lagrangian terms with at most one heavy field *)
	Ltree1 = Plus@@Cases[
		L,
		x_/;Length[Cases[Matchete`PackageScope`RemovePower[x],f_Field/;GetFields[First@f,Heavy],All]]<=1,
		1
	];
	
	(* get kinetic and mass terms of the heavy fields that contribute to the tree-level matching *)
	Ltree2 = Plus@@Cases[
		L,
		x_/;(
		(Length[Cases[Matchete`PackageScope`RemovePower[x],_Field,All]]==2)
		&&
		MatchQ[Cases[Matchete`PackageScope`RemovePower[x],f_Field/;GetFields[First@f,Heavy]:>First[f],All],{y1_,y2_}/;((!FreeQ[Ltree1,y1,All])&&(!FreeQ[Ltree1,y2,All]))]),
		1
	];
	
	(* get all Lagrangian terms relevant for tree-level matching *)
	Ltree = Ltree1 + Ltree2;
	
	(* loop part *)
	Lloop = L-Ltree;
	
	(* sanity check *)
	If[((Ltree+Lloop)-L)=!=0, 
		Message[GetTreeLevelTerms::sanitycheck]; 
		Abort[]
	];
	
	Return@{Ltree,Lloop}
]


(* ::Section::Closed:: *)
(*Helper functions*)


GaugeCouplingsFromRep::usage="returns the gauge coupling corresponding the the gauge group of a given representation.";

GaugeCouplingsFromRep[rep_]:=GetGaugeGroups[Head[rep],Coupling][]


GetFieldGaugeReps::usage="returns all representations of non-Abelian groups under which a given field transforms together with all its charges under Abelian groups.";

GetFieldGaugeReps[\[Psi]_]:=GetFields[\[Psi],Indices]~Join~GetFields[\[Psi],Charges]/.Alternatives@@(Keys@GetGlobalGroups[]~Join~Keys@GetFlavorIndices[])->Nothing


QuadraticCasimir::usage=" returns the the quadratic Casimir operator for a given representation.";

QuadraticCasimir::noCasimir="The quadratic Casimir could not be de termined for `1`.";

QuadraticCasimir[rep_]:=If[MemberQ[Keys@GetGroups[],Head[rep]],
	(* non-Abelian groups *)
	Casimir2[Head[rep]/.GetGroups[],GetRepresentations[rep,DynkinCoefficients]]
	,
	(* Abelian groups *)
	If[MatchQ[rep,(_[n_]/;NumberQ[n])],
		First[rep]^2
		,
		Message[QuadraticCasimir::noCasimir,rep];
		Abort[]
	]
]


SumC2g2::usage="for a given field label this functions returns the quadratic Casimir operator times the square of the associated gauge coupling summed over all reresentations of the field.";

SumC2g2[label_]:=Module[{tmp=GetFieldGaugeReps[label],C2,g},
	C2 = QuadraticCasimir/@ tmp;
	g  = GaugeCouplingsFromRep/@ tmp;
	Plus@@(C2*(g)^2)
]


UnifyOpIndices::usage="for every dummy index contracted within a given operator a Kronecker delta is extracted.";

(* unifies the dummy indices in an operator by introducing new indices and factoring out appropriate Delta functions *)
UnifyOpIndices[oper_Matchete`PackageScope`Operator]:=Module[{op=oper,inds,newInds,deltas,rules},
	inds    = Matchete`PackageScope`FindDummyIndices[op];
	newInds = Table[ind/.First[ind]->Unique[],{ind,inds}];
	deltas  = Times@@(Delta[First@#,Last@#]&/@Transpose@{inds,newInds});
	rules   = (First@#->Last@#)&/@Transpose@{inds,newInds};
	
	(* ensure every index replacement is only performed once *)
	Do[
		op=Matchete`PackageScope`ReplaceFirst[op,rule]
		,
		{rule,rules}
	];
	op*deltas
]


(* ::Section::Closed:: *)
(*Scheme change shift for quartic scalar couplings*)


QuarticScalarShift::usage="takes four scalar fields as input in the form [Bar@\!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\),Bar@\!\(\*SubscriptBox[\(\[Phi]\), \(2\)]\),\!\(\*SubscriptBox[\(\[Phi]\), \(3\)]\),\!\(\*SubscriptBox[\(\[Phi]\), \(4\)]\)] and returns the corresponding finite counterterm originating from the DRbar-MSbar scheme change.";

QuarticScalarShift[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Module[{groups,A,B,a1,b1,c1,d1,e1,f1,a2,b2,c2,d2,e2,f2,ind1,ind2,ind3,ind4,id13$24,id14$23,res,includeDeltas},
	
	(* includes a Delta[a_i,b_i] for indA={a_1,a_2,...} and indB={b_1,b_2,...} if indA and indB have same length, otherwise 0 is returned *)
	includeDeltas[indA_,indB_]:=If[Length[indA]==Length[indB],Times@@(Delta[First[#],Last[#]]&/@Transpose@{indA,indB}),0];
	
	(* markers for when \[Phi]i=\[Phi]j *)
	id13$24=If[(First@\[Phi]1===First@\[Phi]3)&&(First@\[Phi]2===First@\[Phi]4),1,0];
	id14$23=If[(First@\[Phi]1===First@\[Phi]4)&&(First@\[Phi]2===First@\[Phi]3),1,0];
	
	(* list of gauge groups commen between all fields *)
	groups=Intersection[Head/@GetFieldGaugeReps[First@\[Phi]1],
		Head/@GetFieldGaugeReps[First@\[Phi]2],
		Head/@GetFieldGaugeReps[First@\[Phi]3],
		Head/@GetFieldGaugeReps[First@\[Phi]4]
	];
	
	res=-(hbar/4)*Sum[
		(* field indices that belong to the current gauge rep. *)
		a1=FirstCase[\[Phi]1,Index[l_,gr1[_]]:>l,0,All];
		c1=FirstCase[\[Phi]3,Index[l_,gr1[_]]:>l,0,All];
		d1=FirstCase[\[Phi]2,Index[l_,gr1[_]]:>l,0,All];
		f1=FirstCase[\[Phi]4,Index[l_,gr1[_]]:>l,0,All];
		a2=FirstCase[\[Phi]1,Index[l_,gr2[_]]:>l,0,All];
		c2=FirstCase[\[Phi]3,Index[l_,gr2[_]]:>l,0,All];
		d2=FirstCase[\[Phi]2,Index[l_,gr2[_]]:>l,0,All];
		f2=FirstCase[\[Phi]4,Index[l_,gr2[_]]:>l,0,All];
		(* field indices that do not belong to the current gauge reps and that should be contracted by Deltas *)
		ind1=Cases[\[Phi]1,Index[Except[a1|a2],_],All];
		ind2=Cases[\[Phi]2,Index[Except[d1|d2],_],All];
		ind3=Cases[\[Phi]3,Index[Except[c1|c2],_],All];
		ind4=Cases[\[Phi]4,Index[Except[f1|f2],_],All];
		(* include proper power of gauge couplings *)
		GetGaugeGroups[gr1,Coupling][]^2*GetGaugeGroups[gr2,Coupling][]^2*
		(* distinguish Abelian and non-Abelian groups *)
		Switch[{GetGaugeGroups[gr1,Abelian],GetGaugeGroups[gr2,Abelian]},
			(* {Abelian, Abelian} *)
			{True,True},
				(id13$24*includeDeltas[ind1,ind3]*includeDeltas[ind2,ind4]+id14$23*includeDeltas[ind1,ind4]*includeDeltas[ind2,ind3])(
					FieldGenerators[First@\[Phi]1,gr1]FieldGenerators[First@\[Phi]3,gr2]+
					FieldGenerators[First@\[Phi]1,gr2]FieldGenerators[First@\[Phi]3,gr1]
				)*(
					FieldGenerators[First@\[Phi]2,gr1]FieldGenerators[First@\[Phi]4,gr2]+
					FieldGenerators[First@\[Phi]2,gr2]FieldGenerators[First@\[Phi]4,gr1]
				),
			(* {Abelian, non-Abelian} *)
			{True,False},
				id13$24*includeDeltas[ind1,ind3]*includeDeltas[ind2,ind4]*(
					FieldGenerators[First@\[Phi]1,gr1]FieldGenerators[First@\[Phi]3,gr2,{B,a2,c2}]+
					FieldGenerators[First@\[Phi]1,gr2,{B,a2,c2}]FieldGenerators[First@\[Phi]3,gr1]
				)*(
					FieldGenerators[First@\[Phi]2,gr1]FieldGenerators[First@\[Phi]4,gr2,{B,d2,f2}]+
					FieldGenerators[First@\[Phi]2,gr2,{B,d2,f2}]FieldGenerators[First@\[Phi]4,gr1]
				)+
				id14$23*includeDeltas[ind1,ind4]*includeDeltas[ind2,ind3]*(
					FieldGenerators[First@\[Phi]1,gr1]FieldGenerators[First@\[Phi]3,gr2,{B,d2,c2}]+
					FieldGenerators[First@\[Phi]1,gr2,{B,d2,c2}]FieldGenerators[First@\[Phi]3,gr1]
				)*(
					FieldGenerators[First@\[Phi]2,gr1]FieldGenerators[First@\[Phi]4,gr2,{B,a2,f2}]+
					FieldGenerators[First@\[Phi]2,gr2,{B,a2,f2}]FieldGenerators[First@\[Phi]4,gr1]
				),
			(* {non-Abelian, Abelian} *)
			{False,True},
				id13$24*includeDeltas[ind1,ind3]*includeDeltas[ind2,ind4]*(
					FieldGenerators[First@\[Phi]1,gr1,{A,a1,c1}]FieldGenerators[First@\[Phi]3,gr2]+
					FieldGenerators[First@\[Phi]1,gr2]FieldGenerators[First@\[Phi]3,gr1,{A,a1,c1}]
				)*(
					FieldGenerators[First@\[Phi]2,gr1,{A,d1,f1}]FieldGenerators[First@\[Phi]4,gr2]+
					FieldGenerators[First@\[Phi]2,gr2]FieldGenerators[First@\[Phi]4,gr1,{A,d1,f1}]
				)+
				id14$23*includeDeltas[ind1,ind4]*includeDeltas[ind2,ind3]*(
					FieldGenerators[First@\[Phi]1,gr1,{A,d1,c1}]FieldGenerators[First@\[Phi]3,gr2]+
					FieldGenerators[First@\[Phi]1,gr2]FieldGenerators[First@\[Phi]3,gr1,{A,d1,c1}]
				)*(
					FieldGenerators[First@\[Phi]2,gr1,{A,a1,f1}]FieldGenerators[First@\[Phi]4,gr2]+
					FieldGenerators[First@\[Phi]2,gr2]FieldGenerators[First@\[Phi]4,gr1,{A,a1,f1}]
				),
			(* {non-Abelian, non-Abelian} *)
			{False,False},
				If[gr1===gr2,
					(* if both groups are identical *)
					id13$24*includeDeltas[ind1,ind3]*includeDeltas[ind2,ind4]*(
						FieldGenerators[First@\[Phi]1,gr1,{A,a1,b1}]FieldGenerators[First@\[Phi]3,gr1,{B,b1,c1}]+
						FieldGenerators[First@\[Phi]1,gr1,{B,a1,b1}]FieldGenerators[First@\[Phi]3,gr1,{A,b1,c1}]
					)*(
						FieldGenerators[First@\[Phi]2,gr1,{A,d1,e1}]FieldGenerators[First@\[Phi]4,gr1,{B,e1,f1}]+
						FieldGenerators[First@\[Phi]2,gr1,{B,d1,e1}]FieldGenerators[First@\[Phi]4,gr1,{A,e1,f1}]
					)+
					id14$23*includeDeltas[ind1,ind4]*includeDeltas[ind2,ind3]*(
						FieldGenerators[First@\[Phi]1,gr1,{A,d1,b1}]FieldGenerators[First@\[Phi]3,gr1,{B,b1,c1}]+
						FieldGenerators[First@\[Phi]1,gr1,{B,d1,b1}]FieldGenerators[First@\[Phi]3,gr1,{A,b1,c1}]
					)*(
						FieldGenerators[First@\[Phi]2,gr1,{A,a1,e1}]FieldGenerators[First@\[Phi]4,gr1,{B,e1,f1}]+
						FieldGenerators[First@\[Phi]2,gr1,{B,a1,e1}]FieldGenerators[First@\[Phi]4,gr1,{A,e1,f1}]
					)
					,
					(* if both groups are different *)
					id13$24*includeDeltas[ind1,ind3]*includeDeltas[ind2,ind4]*(
						2*FieldGenerators[First@\[Phi]1,gr1,{A,a1,c1}]FieldGenerators[First@\[Phi]3,gr2,{B,a2,c2}]
					)*(
						2*FieldGenerators[First@\[Phi]2,gr1,{A,d1,f1}]FieldGenerators[First@\[Phi]4,gr2,{B,d2,f2}]
					)+
					id14$23*includeDeltas[ind1,ind4]*includeDeltas[ind2,ind3]*(
						2*FieldGenerators[First@\[Phi]1,gr1,{A,d1,c1}]FieldGenerators[First@\[Phi]3,gr2,{B,d2,c2}]
					)*(
						2*FieldGenerators[First@\[Phi]2,gr1,{A,a1,f1}]FieldGenerators[First@\[Phi]4,gr2,{B,a2,f2}]
					)
				]
		]
		,
		{gr1,groups},{gr2,groups}
	]//ContractCGs
]


(* ::Section::Closed:: *)
(*Scheme conversion function: DR-MS*)


Global`DR2MS::noncanonicnorm="The kinetic terms of the gauge bosons are not canonically normalized.";
Global`DR2MS[lagrangian_]:=Module[
	{
		lagTree,lagLoop,
		lag,
		LagGaugeKin,
		LagFermionMass,
		LagYukawa,
		LagQuartic,
		LagNonShifted
		,
		fieldContent,
		gaugeShifts,
		gaugeShiftsInverse={},
		fermionMassShift,
		yukawaShifts,
		fields,
		baredFields,
		quarticShifts
	}
	,
	(* separate terms relevant for tree level matching *)
	{lagTree,lagLoop} = GetTreeLevelTerms[lagrangian];
	
	(* group by operators *)
	lag=CollectOperators[lagTree];
	
	(* get kinetic terms of all gauge bosons and check they are canonically normalized *)
	LagGaugeKin = Matchete`PackageScope`KineticTerms[lag]/._Field->0//Matchete`PackageScope`OperatorToNormalForm;
	If[!Matchete`PackageScope`KineticCanonicalQ[LagGaugeKin],
		Message[Global`DR2MS::noncanonicnorm];
		Abort[]
	];
	LagGaugeKin = Matchete`PackageScope`TermsToList[LagGaugeKin];
	
	(*  Matchete`PackageScope`TermsToList@lag; is problematic here since it expands lag before going to the List *)
	lag=If[Head[lag]===Plus,
		List@@lag,
		{lag}
	];
	
	
	(* get all Fermion mass terms *)
	LagFermionMass=Cases[
		lag/._FieldStrength->0/.Field[__,Except[{}]]->0,
		x_/;((Count[x,Field[_,Fermion,___],All]==2)&&(Count[x,Field[_,Except[Fermion],___],All]==0)),
		1
	];
	
	(* get all scalar-fermion-fermion terms *)
	LagYukawa=Cases[
		lag/._FieldStrength->0,
		x_/;((Count[x,Field[_,Fermion,___],All]==2)&&(Count[Matchete`PackageScope`RemovePower[x],Field[_,Scalar,___],All]==1)),
		1
	];
	
	(* get all scalar quartics *)
	(* NOTE: instead of strating from the present scalar quartics, it might be better to start from the gauge interactions in the availavle sclar kinetic terms! *)
	LagQuartic=Cases[
		lag/._FieldStrength->0,
		x_/;((Count[Matchete`PackageScope`RemovePower[x],Field[_,Scalar,___],All]==4)&&(Count[x,Field[_,Except[Scalar],___],All]==0)),
		1
	];
	
	(* determine all terms that will not be shifted *)
	LagNonShifted=Matchete`PackageScope`LagrangianExpand[(Plus@@lag)-(Plus@@LagGaugeKin+Plus@@LagFermionMass+Plus@@LagYukawa+Plus@@LagQuartic)];
	(*
	Echo[Format[Plus@@LagGaugeKin,NiceForm],"LagGaugeKin"];
	Echo[Format[Plus@@LagFermionMass,NiceForm],"LagFermionMass"];
	Echo[Format[Plus@@LagYukawa,NiceForm],"LagYukawa"];
	Echo[Format[Plus@@LagQuartic,NiceForm],"LagQuartic"];
	Echo[Format[LagNonShifted,NiceForm],"LagNonShifted"];
	*)
	

(* compute scheme changes  *)
	(* proper gauge coupling shifts *)
	If[Length[LagGaugeKin]>0&&LagGaugeKin=!={0},
		gaugeShifts=(FirstCase[#,_Coupling,0,All]->FirstCase[#,_Coupling,0,All](1+hbar/6 SumC2g2[FirstCase[#,FieldStrength[l_,___]:>l,0,All]]))&/@LagGaugeKin;
		(* This shift leads to NON-canonically normalized kinetic terms for the gauge bosons *)
		LagGaugeKin=(First[#]/.Last[#])&/@Transpose@{LagGaugeKin,gaugeShifts};
		(* compute the inverse shift of the MSbar gauge coupling to regain a canonical normalization at the end *)
		gaugeShiftsInverse=gaugeShifts/.hbar->-hbar;
	];

	(* Fermion mass shifts *)
	If[Length[LagFermionMass]>0,
		fermionMassShift=hbar*SumC2g2[FirstCase[#,Field[l_,Fermion,___]:>l,0,All]]&/@LagFermionMass;
		LagFermionMass=(First[#]/.Last[#])&/@Transpose@{
			Matchete`PackageScope`NormalToOperatorForm[LagFermionMass],
			(op_Matchete`PackageScope`Operator:>op(1+#))&/@fermionMassShift
		};
	];

	(* Yukawa mass shifts *)
	If[Length[LagYukawa]>0,
		fieldContent=Join[Cases[#,Field[l_,Fermion,___]:>l,All],Cases[#,Field[l_,Scalar,___]:>l,All]]&/@LagYukawa;
		yukawaShifts=Simplify[hbar/2 (SumC2g2[#[[1]]]+SumC2g2[#[[2]]]-2SumC2g2[#[[3]]])]&/@fieldContent;
		LagYukawa=(First[#]/.Last[#])&/@Transpose@{Matchete`PackageScope`NormalToOperatorForm[LagYukawa],(op_Matchete`PackageScope`Operator:>op(1-#))&/@yukawaShifts};
		LagYukawa=Matchete`PackageScope`OperatorToNormalForm[LagYukawa];
	];
	
	(* shifting scalar quartics *)
	If[Length[LagQuartic]>0,
		LagQuartic=CollectOperators[LagQuartic,Matchete`Simplifications`PackagePrivate`NormalForm->False]/.op_Matchete`PackageScope`Operator:>UnifyOpIndices[op];
		fields=Cases[#,Field[_,Scalar,___],All]&/@LagQuartic;(*all fields*)
		baredFields=Cases[#,Bar[f:Field[_,Scalar,___]]:>f,All]&/@LagQuartic;(*all bared fields*)
		fields=Complement[First@#,Last@#]&/@Transpose@{fields,baredFields};(*all non-bared fields*)
		fields=Join[First@#,Last@#]&/@Transpose@{baredFields,fields};(*sort: {Bar@\[Phi],Bar@\[Phi],\[Phi],\[Phi]}*)
		quarticShifts=QuarticScalarShift[Sequence@@#]&/@fields;
		LagQuartic=(First@#+(Last@#)*FirstCase[First@#,_Matchete`PackageScope`Operator,All])&/@Transpose@{LagQuartic,quarticShifts};
	];
	
	
	(* combine results *)
	Collect[
		ExpandDenominator[
			Contract@RelabelIndices[
				Plus@@LagGaugeKin+Plus@@LagFermionMass+Plus@@LagYukawa+ContractCGs[Plus@@LagQuartic]+LagNonShifted
			]/.gaugeShiftsInverse
		],
		hbar,
		CollectOperators
	]+lagLoop
]


(* ::Subtitle:: *)
(*End of file*)


End[]
