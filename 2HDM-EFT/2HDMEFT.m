(* ::Package:: *)

(* ::Title:: *)
(*2HDMEFT model file*)


(* ::Text:: *)
(*This file provides the a dimensions six EFT basis for 2 Higgs doublet models (2HDMs).*)


(* ::Text:: *)
(*Notice: For using this EFT basis in combination with MapEffectiveCouplings Matchete v0.3.2 or newer is required!*)


(* ::Text:: *)
(*Author: Felix Wilsch*)
(*E-Mail: felix.wilsch@physik.rwth-aachen.de*)
(*Date (created): 17.10.2025*)


(* ::Text:: *)
(*Change log (current version v1):*)
(* - v1: initial release.*)


(* ::Text:: *)
(*This operator basis was first derived in: *)
(*"Two-Higgs Doublet Model Effective Field Theory", *)
(*Radovan Dermisek and Keith Hermanek, *)
(*Phys.Rev.D 110 (2024) 3, 035026, *)
(*[arXiv:2405.20511].*)


(* ::Text:: *)
(*If you use this Matchete model file implementation, please also cite: *)
(*"SUSY meets SMEFT: Complete one-loop matching of the general MSSM",*)
(*Sabine Kraml, Andre Lessa, Suraj Prakash, and Felix Wilsch,*)
(*[arXiv:2506.05201].*)


ParentModel["2HDM"]


(* ::Section:: *)
(*EFT Wilson coefficients*)


(* ::Subsection::Closed:: *)
(*d=5*)


DefineCoupling[c\[Nu]\[Nu]\[Phi]11, 
	Indices->{Flavor, Flavor}, 
	Symmetries-> {SymmetricPermutation[2, 1]},
	NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Nu]\[Nu]\[Phi]\), \((11)\)]\)"
];
DefineCoupling[c\[Nu]\[Nu]\[Phi]22, 
	Indices->{Flavor, Flavor}, 
	Symmetries-> {SymmetricPermutation[2, 1]},
	NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Nu]\[Nu]\[Phi]\), \((22)\)]\)"
];
DefineCoupling[c\[Nu]\[Nu]\[Phi]12, 
	Indices->{Flavor, Flavor}, 
	(*Symmetries-> {SymmetricPermutation[2, 1]},*)
	NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Nu]\[Nu]\[Phi]\), \((12)\)]\)"
];
(* is the last symmetry really absent? *)


(* ::Subsection::Closed:: *)
(*d=6*)


DefineCoupling[ce\[Phi]111, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(e\[Phi]\), \(1\)], \((11)\)]\)"];
DefineCoupling[ce\[Phi]122, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(e\[Phi]\), \(1\)], \((22)\)]\)"];
DefineCoupling[ce\[Phi]121, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(e\[Phi]\), \(1\)], \((21)\)]\)"];
DefineCoupling[ce\[Phi]112, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(e\[Phi]\), \(1\)], \((12)\)]\)"];
DefineCoupling[ce\[Phi]211, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(e\[Phi]\), \(2\)], \((11)\)]\)"];
DefineCoupling[ce\[Phi]222, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(e\[Phi]\), \(2\)], \((22)\)]\)"];
DefineCoupling[ce\[Phi]221, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(e\[Phi]\), \(2\)], \((21)\)]\)"];
DefineCoupling[ce\[Phi]212, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(e\[Phi]\), \(2\)], \((12)\)]\)"];

DefineCoupling[cd\[Phi]111, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(d\[Phi]\), \(1\)], \((11)\)]\)"];
DefineCoupling[cd\[Phi]122, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(d\[Phi]\), \(1\)], \((22)\)]\)"];
DefineCoupling[cd\[Phi]121, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(d\[Phi]\), \(1\)], \((21)\)]\)"];
DefineCoupling[cd\[Phi]112, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(d\[Phi]\), \(1\)], \((12)\)]\)"];
DefineCoupling[cd\[Phi]211, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(d\[Phi]\), \(2\)], \((11)\)]\)"];
DefineCoupling[cd\[Phi]222, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(d\[Phi]\), \(2\)], \((22)\)]\)"];
DefineCoupling[cd\[Phi]221, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(d\[Phi]\), \(2\)], \((21)\)]\)"];
DefineCoupling[cd\[Phi]212, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(d\[Phi]\), \(2\)], \((12)\)]\)"];

DefineCoupling[cu\[Phi]111, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(u\[Phi]\), \(1\)], \((11)\)]\)"];
DefineCoupling[cu\[Phi]122, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(u\[Phi]\), \(1\)], \((22)\)]\)"];
DefineCoupling[cu\[Phi]121, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(u\[Phi]\), \(1\)], \((21)\)]\)"];
DefineCoupling[cu\[Phi]112, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(u\[Phi]\), \(1\)], \((12)\)]\)"];
DefineCoupling[cu\[Phi]211, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(u\[Phi]\), \(2\)], \((11)\)]\)"];
DefineCoupling[cu\[Phi]222, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(u\[Phi]\), \(2\)], \((22)\)]\)"];
DefineCoupling[cu\[Phi]221, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(u\[Phi]\), \(2\)], \((21)\)]\)"];
DefineCoupling[cu\[Phi]212, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(u\[Phi]\), \(2\)], \((12)\)]\)"];


DefineCoupling[ceB\[Phi]1, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(eB\[Phi]\), \(1\)]]\)"];
DefineCoupling[ceW\[Phi]1, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(eW\[Phi]\), \(1\)]]\)"];
DefineCoupling[ceB\[Phi]2, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(eB\[Phi]\), \(2\)]]\)"];
DefineCoupling[ceW\[Phi]2, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(eW\[Phi]\), \(2\)]]\)"];

DefineCoupling[cdB\[Phi]1, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(dB\[Phi]\), \(1\)]]\)"];
DefineCoupling[cdW\[Phi]1, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(dW\[Phi]\), \(1\)]]\)"];
DefineCoupling[cdG\[Phi]1, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(dG\[Phi]\), \(1\)]]\)"];
DefineCoupling[cdB\[Phi]2, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(dB\[Phi]\), \(2\)]]\)"];
DefineCoupling[cdW\[Phi]2, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(dW\[Phi]\), \(2\)]]\)"];
DefineCoupling[cdG\[Phi]2, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(dG\[Phi]\), \(2\)]]\)"];

DefineCoupling[cuB\[Phi]1, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(uB\[Phi]\), \(1\)]]\)"];
DefineCoupling[cuW\[Phi]1, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(uW\[Phi]\), \(1\)]]\)"];
DefineCoupling[cuG\[Phi]1, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(uG\[Phi]\), \(1\)]]\)"];
DefineCoupling[cuB\[Phi]2, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(uB\[Phi]\), \(2\)]]\)"];
DefineCoupling[cuW\[Phi]2, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(uW\[Phi]\), \(2\)]]\)"];
DefineCoupling[cuG\[Phi]2, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), SubscriptBox[\(uG\[Phi]\), \(2\)]]\)"];


DefineCoupling[c\[Phi]e11,  Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]e\), \((11)\)]\)"];
DefineCoupling[c\[Phi]e22,  Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]e\), \((22)\)]\)"];
DefineCoupling[c\[Phi]e12,  Indices-> {Flavor, Flavor}, (*SelfConjugate-> {2, 1},*) NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]e\), \((12)\)]\)"];

DefineCoupling[c\[Phi]l111, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]l\), \(\((11)\)[1]\)]\)"];
DefineCoupling[c\[Phi]l221, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]l\), \(\((22)\)[1]\)]\)"];
DefineCoupling[c\[Phi]l121, Indices-> {Flavor, Flavor}, (*SelfConjugate-> {2, 1},*) NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]l\), \(\((12)\)[1]\)]\)"];
DefineCoupling[c\[Phi]l113, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]l\), \(\((11)\)[3]\)]\)"];
DefineCoupling[c\[Phi]l223, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]l\), \(\((22)\)[3]\)]\)"];
DefineCoupling[c\[Phi]l123, Indices-> {Flavor, Flavor}, (*SelfConjugate-> {2, 1},*) NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]l\), \(\((12)\)[3]\)]\)"];

DefineCoupling[c\[Phi]d11,  Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]d\), \((11)\)]\)"];
DefineCoupling[c\[Phi]d22,  Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]d\), \((22)\)]\)"];
DefineCoupling[c\[Phi]d12,  Indices-> {Flavor, Flavor}, (*SelfConjugate-> {2, 1},*) NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]d\), \((12)\)]\)"];

DefineCoupling[c\[Phi]u11,  Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]u\), \((11)\)]\)"];
DefineCoupling[c\[Phi]u22,  Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]u\), \((22)\)]\)"];
DefineCoupling[c\[Phi]u12,  Indices-> {Flavor, Flavor}, (*SelfConjugate-> {2, 1},*) NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]u\), \((12)\)]\)"];

DefineCoupling[c\[Phi]q111, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]q\), \(\((11)\)[1]\)]\)"];
DefineCoupling[c\[Phi]q221, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]q\), \(\((22)\)[1]\)]\)"];
DefineCoupling[c\[Phi]q121, Indices-> {Flavor, Flavor}, (*SelfConjugate-> {2, 1},*) NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]q\), \(\((12)\)[1]\)]\)"];
DefineCoupling[c\[Phi]q113, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]q\), \(\((11)\)[3]\)]\)"];
DefineCoupling[c\[Phi]q223, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]q\), \(\((22)\)[3]\)]\)"];
DefineCoupling[c\[Phi]q123, Indices-> {Flavor, Flavor}, (*SelfConjugate-> {2, 1},*) NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]q\), \(\((12)\)[3]\)]\)"];

DefineCoupling[c\[Phi]ud11, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]ud\), \((11)\)]\)"];
DefineCoupling[c\[Phi]ud22, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]ud\), \((22)\)]\)"];
DefineCoupling[c\[Phi]ud12, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]ud\), \((12)\)]\)"];


DefineCoupling[cG,  SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(G\)]\)"];
DefineCoupling[cGt, SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), OverscriptBox[\(G\), \(~\)]]\)"];
DefineCoupling[cW,  SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(W\)]\)"];
DefineCoupling[cWt, SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), OverscriptBox[\(W\), \(~\)]]\)"];


DefineCoupling[c\[Phi]G11,  SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]G\), \((11)\)]\)"];
DefineCoupling[c\[Phi]G22,  SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]G\), \((22)\)]\)"];
DefineCoupling[c\[Phi]G21,  SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]G\), \((21)\)]\)"];
DefineCoupling[c\[Phi]Gt11, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(G\), \(~\)]\), \((11)\)]\)"];
DefineCoupling[c\[Phi]Gt22, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(G\), \(~\)]\), \((22)\)]\)"];
DefineCoupling[c\[Phi]Gt21, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(G\), \(~\)]\), \((21)\)]\)"];

DefineCoupling[c\[Phi]B11,  SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]B\), \((11)\)]\)"];
DefineCoupling[c\[Phi]B22,  SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]B\), \((22)\)]\)"];
DefineCoupling[c\[Phi]B21,  SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]B\), \((21)\)]\)"];
DefineCoupling[c\[Phi]Bt11, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(B\), \(~\)]\), \((11)\)]\)"];
DefineCoupling[c\[Phi]Bt22, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(B\), \(~\)]\), \((22)\)]\)"];
DefineCoupling[c\[Phi]Bt21, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(B\), \(~\)]\), \((21)\)]\)"];

DefineCoupling[c\[Phi]W11,  SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]W\), \((11)\)]\)"];
DefineCoupling[c\[Phi]W22,  SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]W\), \((22)\)]\)"];
DefineCoupling[c\[Phi]W21,  SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]W\), \((21)\)]\)"];
DefineCoupling[c\[Phi]Wt11, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(W\), \(~\)]\), \((11)\)]\)"];
DefineCoupling[c\[Phi]Wt22, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(W\), \(~\)]\), \((22)\)]\)"];
DefineCoupling[c\[Phi]Wt21, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(W\), \(~\)]\), \((21)\)]\)"];

DefineCoupling[c\[Phi]WB11,  SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]WB\), \((11)\)]\)"];
DefineCoupling[c\[Phi]WB22,  SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]WB\), \((22)\)]\)"];
DefineCoupling[c\[Phi]WB21,  SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]WB\), \((21)\)]\)"];
DefineCoupling[c\[Phi]WtB11, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(W\), \(~\)] B\), \((11)\)]\)"];
DefineCoupling[c\[Phi]WtB22, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(W\), \(~\)] B\), \((22)\)]\)"];
DefineCoupling[c\[Phi]WtB21, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*OverscriptBox[\(W\), \(~\)] B\), \((21)\)]\)"];


DefineCoupling[c\[Phi]PD1111, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*SuperscriptBox[\(\[PartialD]\), \(2\)]\), \(\((11)\) \((11)\)\)]\)"];
DefineCoupling[c\[Phi]PD2222, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*SuperscriptBox[\(\[PartialD]\), \(2\)]\), \(\((22)\) \((22)\)\)]\)"];
DefineCoupling[c\[Phi]PD1122, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*SuperscriptBox[\(\[PartialD]\), \(2\)]\), \(\((11)\) \((22)\)\)]\)"];
DefineCoupling[c\[Phi]PD2121, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*SuperscriptBox[\(\[PartialD]\), \(2\)]\), \(\((21)\) \((21)\)\)]\)"];
DefineCoupling[c\[Phi]PD2112, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*SuperscriptBox[\(\[PartialD]\), \(2\)]\), \(\((21)\) \((12)\)\)]\)"];
DefineCoupling[c\[Phi]PD2111, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*SuperscriptBox[\(\[PartialD]\), \(2\)]\), \(\((21)\) \((11)\)\)]\)"];
DefineCoupling[c\[Phi]PD2122, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi] \*SuperscriptBox[\(\[PartialD]\), \(2\)]\), \(\((21)\) \((22)\)\)]\)"];

DefineCoupling[c\[Phi]D1111, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]D\), \(\((11)\) \((11)\)\)]\)"];
DefineCoupling[c\[Phi]D2222, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]D\), \(\((22)\) \((22)\)\)]\)"];
DefineCoupling[c\[Phi]D1122, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]D\), \(\((11)\) \((22)\)\)]\)"];
DefineCoupling[c\[Phi]D2121, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]D\), \(\((21)\) \((21)\)\)]\)"];
DefineCoupling[c\[Phi]D2112, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]D\), \(\((21)\) \((12)\)\)]\)"];
DefineCoupling[c\[Phi]D2111, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]D\), \(\((21)\) \((11)\)\)]\)"];
DefineCoupling[c\[Phi]D2122, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]D\), \(\((21)\) \((22)\)\)]\)"];


DefineCoupling[c\[Phi]111111, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((11)\) \((11)\) \((11)\)\)]\)"];
DefineCoupling[c\[Phi]111122, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((11)\) \((11)\) \((22)\)\)]\)"];
DefineCoupling[c\[Phi]112222, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((11)\) \((22)\) \((22)\)\)]\)"];
DefineCoupling[c\[Phi]111121, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((11)\) \((11)\) \((21)\)\)]\)"];
DefineCoupling[c\[Phi]222221, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((22)\) \((22)\) \((21)\)\)]\)"];
DefineCoupling[c\[Phi]222222, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((22)\) \((22)\) \((22)\)\)]\)"];

DefineCoupling[c\[Phi]112121, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((11)\) \((21)\) \((21)\)\)]\)"];
DefineCoupling[c\[Phi]112112, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((11)\) \((21)\) \((12)\)\)]\)"];
DefineCoupling[c\[Phi]222121, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((22)\) \((21)\) \((21)\)\)]\)"];
DefineCoupling[c\[Phi]222112, SelfConjugate-> True,  NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((22)\) \((21)\) \((12)\)\)]\)"];
DefineCoupling[c\[Phi]212121, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((21)\) \((21)\) \((21)\)\)]\)"];
DefineCoupling[c\[Phi]212112, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((21)\) \((21)\) \((12)\)\)]\)"];
DefineCoupling[c\[Phi]112221, SelfConjugate-> False, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\[Phi]\), \(\((11)\) \((22)\) \((21)\)\)]\)"];


DefineCoupling[cll,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2]}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(ll\)]\)"];
DefineCoupling[cqq1, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2]}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qq\), \((1)\)]\)"];
DefineCoupling[cqq3, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2]}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qq\), \((3)\)]\)"];
DefineCoupling[clq1, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lq\), \((1)\)]\)"];
DefineCoupling[clq3, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lq\), \((3)\)]\)"];
DefineCoupling[cee,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2], SymmetricPermutation[1, 4, 3, 2]}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(ee\)]\)"];
DefineCoupling[cuu,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2]}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(uu\)]\)"];
DefineCoupling[cdd,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2]}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(dd\)]\)"];
DefineCoupling[ceu,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(eu\)]\)"];
DefineCoupling[ced,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(ed\)]\)"];
DefineCoupling[cud1, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(ud\), \((1)\)]\)"];
DefineCoupling[cud8, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(ud\), \((8)\)]\)"];
DefineCoupling[cle,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(le\)]\)"];
DefineCoupling[clu,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(lu\)]\)"];
DefineCoupling[cld,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(ld\)]\)"];
DefineCoupling[cqe,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(qe\)]\)"];
DefineCoupling[cqu1, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qu\), \((1)\)]\)"];
DefineCoupling[cqu8, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qu\), \((8)\)]\)"];
DefineCoupling[cqd1, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qd\), \((1)\)]\)"];
DefineCoupling[cqd8, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qd\), \((8)\)]\)"];


DefineCoupling[cduq, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(duq\)]\)"];
DefineCoupling[cqqu, Indices-> {Flavor, Flavor, Flavor, Flavor}, Symmetries-> {SymmetricPermutation[2, 1, 3, 4]}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(qqu\)]\)"];
DefineCoupling[cqqq, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(qqq\)]\)"];
DefineCoupling[cduu, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(duu\)]\)"];


DefineCoupling[cledq,  Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(ledq\)]\)"];
DefineCoupling[cquqd1, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quqd\), \((1)\)]\)"];
DefineCoupling[cquqd8, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quqd\), \((8)\)]\)"];
DefineCoupling[clequ1, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequ\), \((1)\)]\)"];
DefineCoupling[clequ3, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequ\), \((3)\)]\)"];


(* ::Section:: *)
(*Lagrangians*)


Module[
	{
		p,r,s,t,i,j,k,m,n,o,J,K,L,\[Alpha],\[Beta],\[Delta],A,C,D,\[Nu],\[Rho],\[Theta],\[Eta],\[Kappa],
		\[Tau]SU2L,fSU2L,\[Epsilon]SU2L,TSU3c,fSU3c,\[Epsilon]SU3c,HermitianCD,
		Lag, Lag5, Lag6\[Psi]2\[Phi]3, Lag6\[Psi]2X\[Phi], Lag6\[Psi]2\[Phi]2D, Lag6X3, Lag6X2\[Phi]2, Lag6\[Phi]4D2, Lag6\[Phi]6, Lag6\[Psi]4, Lag6BNV
	}
	,
	\[Tau]SU2L[Jadj_,ifund_,jfund_]:= 2 CG[gen[SU2L[fund]],{Jadj,ifund,jfund}];
	fSU2L[Iadj_,Jadj_,Kadj_]:=CG[fStruct[SU2L],{Iadj,Jadj,Kadj}];
	\[Epsilon]SU2L[ifund_,jfund_]:= CG[eps[SU2L],{ifund,jfund}];
	TSU3c[Aadj_,\[Alpha]fund_,\[Beta]fund_]:=CG[gen[SU3c[fund]],{Aadj,\[Alpha]fund,\[Beta]fund}];
	fSU3c[Aadj_,Badj_,Cadj_]:=CG[fStruct[SU3c],{Aadj,Badj,Cadj}];
	\[Epsilon]SU3c[\[Alpha]fund_,\[Beta]fund_,\[Gamma]fund_]:= CG[eps[SU3c],{\[Alpha]fund,\[Beta]fund,\[Gamma]fund}];
	HermitianCD[ind_,term1_,term2_]:= I term1 CD[ind,term2] - I CD[ind,term1] term2;

	Lag5 = c\[Nu]\[Nu]\[Phi]11[p,r] Bar@\[Epsilon]SU2L[i,j]Bar@\[Epsilon]SU2L[k,m] \[Phi]1[i] \[Phi]1[k] Bar@CConj@l[j,p]**l[m,r] +
		c\[Nu]\[Nu]\[Phi]22[p,r] Bar@\[Epsilon]SU2L[i,j]Bar@\[Epsilon]SU2L[k,m] \[Phi]2[i] \[Phi]2[k] Bar@CConj@l[j,p]**l[m,r] + 
		c\[Nu]\[Nu]\[Phi]12[p,r] Bar@\[Epsilon]SU2L[i,j]Bar@\[Epsilon]SU2L[k,m] \[Phi]1[i] \[Phi]2[k] Bar@CConj@l[j,p]**l[m,r];
		
	Lag6\[Psi]2\[Phi]3 = ce\[Phi]111[p,r] (Bar@\[Phi]1[i] \[Phi]1[i]) \[Phi]1[j] Bar@l[j,p]**e[r]
		+ ce\[Phi]122[p,r] (Bar@\[Phi]2[i] \[Phi]2[i]) \[Phi]1[j] Bar@l[j,p]**e[r]
		+ ce\[Phi]121[p,r] (Bar@\[Phi]2[i] \[Phi]1[i]) \[Phi]1[j] Bar@l[j,p]**e[r]
		+ ce\[Phi]112[p,r] (Bar@\[Phi]1[i] \[Phi]2[i]) \[Phi]1[j] Bar@l[j,p]**e[r]
		+ ce\[Phi]211[p,r] (Bar@\[Phi]1[i] \[Phi]1[i]) \[Phi]2[j] Bar@l[j,p]**e[r]
		+ ce\[Phi]222[p,r] (Bar@\[Phi]2[i] \[Phi]2[i]) \[Phi]2[j] Bar@l[j,p]**e[r]
		+ ce\[Phi]221[p,r] (Bar@\[Phi]2[i] \[Phi]1[i]) \[Phi]2[j] Bar@l[j,p]**e[r]
		+ ce\[Phi]212[p,r] (Bar@\[Phi]1[i] \[Phi]2[i]) \[Phi]2[j] Bar@l[j,p]**e[r]
		+ cd\[Phi]111[p,r] (Bar@\[Phi]1[i] \[Phi]1[i]) \[Phi]1[j] Bar@q[\[Alpha],j,p]**d[\[Alpha],r]
		+ cd\[Phi]122[p,r] (Bar@\[Phi]2[i] \[Phi]2[i]) \[Phi]1[j] Bar@q[\[Alpha],j,p]**d[\[Alpha],r]
		+ cd\[Phi]121[p,r] (Bar@\[Phi]2[i] \[Phi]1[i]) \[Phi]1[j] Bar@q[\[Alpha],j,p]**d[\[Alpha],r]
		+ cd\[Phi]112[p,r] (Bar@\[Phi]1[i] \[Phi]2[i]) \[Phi]1[j] Bar@q[\[Alpha],j,p]**d[\[Alpha],r]
		+ cd\[Phi]211[p,r] (Bar@\[Phi]1[i] \[Phi]1[i]) \[Phi]2[j] Bar@q[\[Alpha],j,p]**d[\[Alpha],r]
		+ cd\[Phi]222[p,r] (Bar@\[Phi]2[i] \[Phi]2[i]) \[Phi]2[j] Bar@q[\[Alpha],j,p]**d[\[Alpha],r]
		+ cd\[Phi]221[p,r] (Bar@\[Phi]2[i] \[Phi]1[i]) \[Phi]2[j] Bar@q[\[Alpha],j,p]**d[\[Alpha],r]
		+ cd\[Phi]212[p,r] (Bar@\[Phi]1[i] \[Phi]2[i]) \[Phi]2[j] Bar@q[\[Alpha],j,p]**d[\[Alpha],r]
		+ cu\[Phi]111[p,r] (Bar@\[Phi]1[i] \[Phi]1[i]) \[Epsilon]SU2L[j,k] Bar@\[Phi]1[k] Bar@q[\[Alpha],j,p]**u[\[Alpha],r]
		+ cu\[Phi]122[p,r] (Bar@\[Phi]2[i] \[Phi]2[i]) \[Epsilon]SU2L[j,k] Bar@\[Phi]1[k] Bar@q[\[Alpha],j,p]**u[\[Alpha],r]
		+ cu\[Phi]121[p,r] (Bar@\[Phi]2[i] \[Phi]1[i]) \[Epsilon]SU2L[j,k] Bar@\[Phi]1[k] Bar@q[\[Alpha],j,p]**u[\[Alpha],r]
		+ cu\[Phi]112[p,r] (Bar@\[Phi]1[i] \[Phi]2[i]) \[Epsilon]SU2L[j,k] Bar@\[Phi]1[k] Bar@q[\[Alpha],j,p]**u[\[Alpha],r]
		+ cu\[Phi]211[p,r] (Bar@\[Phi]1[i] \[Phi]1[i]) \[Epsilon]SU2L[j,k] Bar@\[Phi]2[k] Bar@q[\[Alpha],j,p]**u[\[Alpha],r]
		+ cu\[Phi]222[p,r] (Bar@\[Phi]2[i] \[Phi]2[i]) \[Epsilon]SU2L[j,k] Bar@\[Phi]2[k] Bar@q[\[Alpha],j,p]**u[\[Alpha],r]
		+ cu\[Phi]221[p,r] (Bar@\[Phi]2[i] \[Phi]1[i]) \[Epsilon]SU2L[j,k] Bar@\[Phi]2[k] Bar@q[\[Alpha],j,p]**u[\[Alpha],r]
		+ cu\[Phi]212[p,r] (Bar@\[Phi]1[i] \[Phi]2[i]) \[Epsilon]SU2L[j,k] Bar@\[Phi]2[k] Bar@q[\[Alpha],j,p]**u[\[Alpha],r];
			
	Lag6\[Psi]2X\[Phi] = ceB\[Phi]1[p,r] Bar@l[i,p]**\[Sigma][\[Nu],\[Rho]]**e[r] \[Phi]1[i] FS[B,\[Nu],\[Rho]]/gY[]
		+ ceW\[Phi]1[p,r] Bar@l[i,p]**\[Sigma][\[Nu],\[Rho]]**e[r] \[Tau]SU2L[J,i,j] \[Phi]1[j] FS[W,\[Nu],\[Rho],J]/gL[]
		+ ceB\[Phi]2[p,r] Bar@l[i,p]**\[Sigma][\[Nu],\[Rho]]**e[r] \[Phi]2[i] FS[B,\[Nu],\[Rho]]/gY[]
		+ ceW\[Phi]2[p,r] Bar@l[i,p]**\[Sigma][\[Nu],\[Rho]]**e[r] \[Tau]SU2L[J,i,j] \[Phi]2[j] FS[W,\[Nu],\[Rho],J]/gL[]
		+ cdB\[Phi]1[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**d[\[Alpha],r] \[Phi]1[i] FS[B,\[Nu],\[Rho]]/gY[]
		+ cdW\[Phi]1[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**d[\[Alpha],r] \[Tau]SU2L[J,i,j] \[Phi]1[j] FS[W,\[Nu],\[Rho],J]/gL[]
		+ cdG\[Phi]1[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**d[\[Beta],r] \[Phi]1[i] TSU3c[A,\[Alpha],\[Beta]]FS[G,\[Nu],\[Rho],A]/gs[]
		+ cdB\[Phi]2[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**d[\[Alpha],r] \[Phi]2[i] FS[B,\[Nu],\[Rho]]/gY[]
		+ cdW\[Phi]2[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**d[\[Alpha],r] \[Tau]SU2L[J,i,j] \[Phi]2[j] FS[W,\[Nu],\[Rho],J]/gL[]
		+ cdG\[Phi]2[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**d[\[Beta],r] \[Phi]2[i] TSU3c[A,\[Alpha],\[Beta]]FS[G,\[Nu],\[Rho],A]/gs[]
		+ cuB\[Phi]1[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**u[\[Alpha],r] \[Epsilon]SU2L[i,j] Bar@\[Phi]1[j] FS[B,\[Nu],\[Rho]]/gY[]
		+ cuW\[Phi]1[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**u[\[Alpha],r] \[Tau]SU2L[J,i,j]\[Epsilon]SU2L[j,k]Bar@\[Phi]1[k] FS[W,\[Nu],\[Rho],J]/gL[]
		+ cuG\[Phi]1[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**u[\[Beta],r]TSU3c[A,\[Alpha],\[Beta]] \[Epsilon]SU2L[i,j] Bar@\[Phi]1[j] FS[G,\[Nu],\[Rho],A]/gs[]
		+ cuB\[Phi]2[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**u[\[Alpha],r] \[Epsilon]SU2L[i,j] Bar@\[Phi]2[j] FS[B,\[Nu],\[Rho]]/gY[]
		+ cuW\[Phi]2[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**u[\[Alpha],r] \[Tau]SU2L[J,i,j]\[Epsilon]SU2L[j,k]Bar@\[Phi]2[k] FS[W,\[Nu],\[Rho],J]/gL[]
		+ cuG\[Phi]2[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**u[\[Beta],r]TSU3c[A,\[Alpha],\[Beta]] \[Epsilon]SU2L[i,j] Bar@\[Phi]2[j] FS[G,\[Nu],\[Rho],A]/gs[];
		
	Lag6\[Psi]2\[Phi]2D = c\[Phi]e11[p,r] HermitianCD[\[Nu], Bar@\[Phi]1[i], \[Phi]1[i]] Bar@e[p]**\[Gamma][\[Nu]]**e[r]
		+ c\[Phi]e22[p,r] HermitianCD[\[Nu], Bar@\[Phi]2[i],\[Phi]2[i]] Bar@e[p]**\[Gamma][\[Nu]]**e[r]
		+ PlusHc[c\[Phi]e12[p,r] HermitianCD[\[Nu], Bar@\[Phi]1[i], \[Phi]2[i]] Bar@e[p]**\[Gamma][\[Nu]]**e[r]]
		+ c\[Phi]l111[p,r] HermitianCD[\[Nu], Bar@\[Phi]1[i],\[Phi]1[i]] Bar@l[j,p]**\[Gamma][\[Nu]]**l[j,r]
		+ c\[Phi]l221[p,r] HermitianCD[\[Nu], Bar@\[Phi]2[i],\[Phi]2[i]] Bar@l[j,p]**\[Gamma][\[Nu]]**l[j,r]
		+ PlusHc[c\[Phi]l121[p,r] HermitianCD[\[Nu], Bar@\[Phi]1[i],\[Phi]2[i]] Bar@l[j,p]**\[Gamma][\[Nu]]**l[j,r]]
		+ c\[Phi]l113[p,r] HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Tau]SU2L[J,i,j] \[Phi]1[j]] \[Tau]SU2L[J,k,m] Bar@l[k,p]**\[Gamma][\[Nu]]**l[m,r]
		+ c\[Phi]l223[p,r] HermitianCD[\[Nu],Bar@\[Phi]2[i],\[Tau]SU2L[J,i,j] \[Phi]2[j]] \[Tau]SU2L[J,k,m] Bar@l[k,p]**\[Gamma][\[Nu]]**l[m,r]
		+ PlusHc[c\[Phi]l123[p,r] HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Tau]SU2L[J,i,j] \[Phi]2[j]] \[Tau]SU2L[J,k,m] Bar@l[k,p]**\[Gamma][\[Nu]]**l[m,r]]
		+ c\[Phi]d11[p,r] HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Phi]1[i]] Bar@d[\[Alpha],p]**\[Gamma][\[Nu]]**d[\[Alpha],r]
		+ c\[Phi]d22[p,r] HermitianCD[\[Nu],Bar@\[Phi]2[i],\[Phi]2[i]] Bar@d[\[Alpha],p]**\[Gamma][\[Nu]]**d[\[Alpha],r]
		+ PlusHc[c\[Phi]d12[p,r] HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Phi]2[i]] Bar@d[\[Alpha],p]**\[Gamma][\[Nu]]**d[\[Alpha],r]]
		+ c\[Phi]u11[p,r] HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Phi]1[i]] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**u[\[Alpha],r]
		+ c\[Phi]u22[p,r] HermitianCD[\[Nu],Bar@\[Phi]2[i],\[Phi]2[i]] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**u[\[Alpha],r]
		+ PlusHc[c\[Phi]u12[p,r] HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Phi]2[i]] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**u[\[Alpha],r]]
		+ c\[Phi]q111[p,r] HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Phi]1[i]] Bar@q[\[Alpha],j,p]**\[Gamma][\[Nu]]**q[\[Alpha],j,r]
		+ c\[Phi]q221[p,r] HermitianCD[\[Nu],Bar@\[Phi]2[i],\[Phi]2[i]] Bar@q[\[Alpha],j,p]**\[Gamma][\[Nu]]**q[\[Alpha],j,r]
		+ PlusHc[c\[Phi]q121[p,r] HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Phi]2[i]] Bar@q[\[Alpha],j,p]**\[Gamma][\[Nu]]**q[\[Alpha],j,r]]
		+ c\[Phi]q113[p,r] HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Tau]SU2L[J,i,j] \[Phi]1[j]] \[Tau]SU2L[J,k,m] Bar@q[\[Alpha],k,p]**\[Gamma][\[Nu]]**q[\[Alpha],m,r]
		+ c\[Phi]q223[p,r] HermitianCD[\[Nu],Bar@\[Phi]2[i],\[Tau]SU2L[J,i,j] \[Phi]2[j]] \[Tau]SU2L[J,k,m] Bar@q[\[Alpha],k,p]**\[Gamma][\[Nu]]**q[\[Alpha],m,r]
		+ PlusHc[c\[Phi]q123[p,r] HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Tau]SU2L[J,i,j] \[Phi]2[j]] \[Tau]SU2L[J,k,m] Bar@q[\[Alpha],k,p]**\[Gamma][\[Nu]]**q[\[Alpha],m,r]]
		+ PlusHc[c\[Phi]ud11[p,r] I Bar@\[Epsilon]SU2L[i,j] \[Phi]1[i] CD[\[Nu],\[Phi]1[j]] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**d[\[Alpha],r]]
		+ PlusHc[c\[Phi]ud22[p,r] I Bar@\[Epsilon]SU2L[i,j] \[Phi]2[i] CD[\[Nu],\[Phi]2[j]] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**d[\[Alpha],r]]
		+ PlusHc[c\[Phi]ud12[p,r] Bar@\[Epsilon]SU2L[i,j] HermitianCD[\[Nu],\[Phi]2[i],\[Phi]1[j]] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**d[\[Alpha],r]];
		
	Lag6X3 = cG[] fSU3c[A,C,D] FS[G, \[Nu], \[Rho], A] FS[G, \[Rho], \[Theta], C] FS[G, \[Theta], \[Nu], D]/gs[]^3
		+ 1/2 cGt[] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]]fSU3c[A,C,D] FS[G, \[Eta], \[Kappa], A] FS[G, \[Rho], \[Theta], C] FS[G, \[Theta], \[Nu], D]/gs[]^3
	    + cW[] fSU2L[J,K,L] FS[W, \[Nu], \[Rho], J] FS[W, \[Rho], \[Theta], K] FS[W, \[Theta], \[Nu], L]/gL[]^3
	    + 1/2 cWt[] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]]fSU2L[J,K,L] FS[W, \[Eta], \[Kappa], J] FS[W, \[Rho], \[Theta], K] FS[W, \[Theta], \[Nu], L]/gL[]^3;
	    
	Lag6X2\[Phi]2 = c\[Phi]G11[] Bar@\[Phi]1[i] \[Phi]1[i] FS[G, \[Nu], \[Rho], A] FS[G, \[Nu], \[Rho], A]/gs[]^2
		+ c\[Phi]G22[] Bar@\[Phi]2[i] \[Phi]2[i] FS[G, \[Nu], \[Rho], A] FS[G, \[Nu], \[Rho], A]/gs[]^2
		+ PlusHc[c\[Phi]G21[] Bar@\[Phi]2[i] \[Phi]1[i] FS[G, \[Nu], \[Rho], A] FS[G, \[Nu], \[Rho], A]/gs[]^2]
		+ 1/2 c\[Phi]Gt11[] Bar@\[Phi]1[i] \[Phi]1[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[G, \[Eta], \[Kappa], A] FS[G, \[Nu], \[Rho], A]/gs[]^2
		+ 1/2 c\[Phi]Gt22[] Bar@\[Phi]2[i] \[Phi]2[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[G, \[Eta], \[Kappa], A] FS[G, \[Nu], \[Rho], A]/gs[]^2
		+ PlusHc[1/2 c\[Phi]Gt21[] Bar@\[Phi]2[i] \[Phi]1[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[G, \[Eta], \[Kappa], A] FS[G, \[Nu], \[Rho], A]/gs[]^2]
		+ c\[Phi]B11[] Bar@\[Phi]1[i] \[Phi]1[i] FS[B, \[Nu], \[Rho]] FS[B, \[Nu], \[Rho]]/gY[]^2
		+ c\[Phi]B22[] Bar@\[Phi]2[i] \[Phi]2[i] FS[B, \[Nu], \[Rho]] FS[B, \[Nu], \[Rho]]/gY[]^2
		+ PlusHc[c\[Phi]B21[] Bar@\[Phi]2[i] \[Phi]1[i] FS[B, \[Nu], \[Rho]] FS[B, \[Nu], \[Rho]]/gY[]^2]
		+ 1/2 c\[Phi]Bt11[] Bar@\[Phi]1[i] \[Phi]1[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[B, \[Eta], \[Kappa]] FS[B, \[Nu], \[Rho]]/gY[]^2
		+ 1/2 c\[Phi]Bt22[] Bar@\[Phi]2[i] \[Phi]2[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[B, \[Eta], \[Kappa]] FS[B, \[Nu], \[Rho]]/gY[]^2
		+ PlusHc[1/2 c\[Phi]Bt21[] Bar@\[Phi]2[i] \[Phi]1[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[B, \[Eta], \[Kappa]] FS[B, \[Nu], \[Rho]]/gY[]^2]
		+ c\[Phi]W11[] Bar@\[Phi]1[i] \[Phi]1[i] FS[W, \[Nu], \[Rho], J] FS[W, \[Nu], \[Rho], J]/gL[]^2
		+ c\[Phi]W22[] Bar@\[Phi]2[i] \[Phi]2[i] FS[W, \[Nu], \[Rho], J] FS[W, \[Nu], \[Rho], J]/gL[]^2
		+ PlusHc[c\[Phi]W21[] Bar@\[Phi]2[i] \[Phi]1[i] FS[W, \[Nu], \[Rho], J] FS[W, \[Nu], \[Rho], J]/gL[]^2]
		+ 1/2 c\[Phi]Wt11[] Bar@\[Phi]1[i] \[Phi]1[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[W, \[Eta], \[Kappa], J] FS[W, \[Nu], \[Rho], J]/gL[]^2
		+ 1/2 c\[Phi]Wt22[] Bar@\[Phi]2[i] \[Phi]2[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[W, \[Eta], \[Kappa], J] FS[W, \[Nu], \[Rho], J]/gL[]^2
		+ PlusHc[1/2 c\[Phi]Wt21[] Bar@\[Phi]2[i] \[Phi]1[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[W, \[Eta], \[Kappa], J] FS[W, \[Nu], \[Rho], J]/gL[]^2]
		+ c\[Phi]WB11[] Bar@\[Phi]1[i] \[Tau]SU2L[J,i,j] \[Phi]1[j] FS[W, \[Nu], \[Rho], J] FS[B, \[Nu], \[Rho]]/(gL[]gY[])
		+ c\[Phi]WB22[] Bar@\[Phi]2[i] \[Tau]SU2L[J,i,j] \[Phi]2[j] FS[W, \[Nu], \[Rho], J] FS[B, \[Nu], \[Rho]]/(gL[]gY[])
		+ PlusHc[c\[Phi]WB21[] Bar@\[Phi]2[i] \[Tau]SU2L[J,i,j] \[Phi]1[j] FS[W, \[Nu], \[Rho], J] FS[B, \[Nu], \[Rho]]/(gL[]gY[])]
		+ 1/2 c\[Phi]WtB11[] Bar@\[Phi]1[i] \[Tau]SU2L[J,i,j] \[Phi]1[j] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[W, \[Eta], \[Kappa], J] FS[B, \[Nu], \[Rho]]/(gL[]gY[])
		+ 1/2 c\[Phi]WtB22[] Bar@\[Phi]2[i] \[Tau]SU2L[J,i,j] \[Phi]2[j] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[W, \[Eta], \[Kappa], J] FS[B, \[Nu], \[Rho]]/(gL[]gY[])
		+ PlusHc[1/2 c\[Phi]WtB21[] Bar@\[Phi]2[i] \[Tau]SU2L[J,i,j] \[Phi]1[j] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[W, \[Eta], \[Kappa], J] FS[B, \[Nu], \[Rho]]/(gL[]gY[])];
	
	Lag6\[Phi]4D2 = c\[Phi]PD1111[] CD[{\[Nu]},Bar@\[Phi]1[i]\[Phi]1[i]]CD[{\[Nu]},Bar@\[Phi]1[j]\[Phi]1[j]]
		+ c\[Phi]PD2222[] CD[{\[Nu]},Bar@\[Phi]2[i]\[Phi]2[i]]CD[{\[Nu]},Bar@\[Phi]2[j]\[Phi]2[j]]
		+ c\[Phi]PD1122[] CD[{\[Nu]},Bar@\[Phi]1[i]\[Phi]1[i]]CD[{\[Nu]},Bar@\[Phi]2[j]\[Phi]2[j]]
		+ PlusHc[c\[Phi]PD2121[] CD[{\[Nu]},Bar@\[Phi]2[i]\[Phi]1[i]]CD[{\[Nu]},Bar@\[Phi]2[j]\[Phi]1[j]]]
		+ c\[Phi]PD2112[] CD[{\[Nu]},Bar@\[Phi]2[i]\[Phi]1[i]]CD[{\[Nu]},Bar@\[Phi]1[j]\[Phi]2[j]]
		+ PlusHc[c\[Phi]PD2111[] CD[{\[Nu]},Bar@\[Phi]2[i]\[Phi]1[i]]CD[{\[Nu]},Bar@\[Phi]1[j]\[Phi]1[j]]]
		+ PlusHc[c\[Phi]PD2122[] CD[{\[Nu]},Bar@\[Phi]2[i]\[Phi]1[i]]CD[{\[Nu]},Bar@\[Phi]2[j]\[Phi]2[j]]]
		- c\[Phi]D1111[]HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Phi]1[i]]HermitianCD[\[Nu],Bar@\[Phi]1[j],\[Phi]1[j]]
		- c\[Phi]D2222[]HermitianCD[\[Nu],Bar@\[Phi]2[i],\[Phi]2[i]]HermitianCD[\[Nu],Bar@\[Phi]2[j],\[Phi]2[j]]
		- c\[Phi]D1122[]HermitianCD[\[Nu],Bar@\[Phi]1[i],\[Phi]1[i]]HermitianCD[\[Nu],Bar@\[Phi]2[j],\[Phi]2[j]]
		- PlusHc[c\[Phi]D2121[]HermitianCD[\[Nu],Bar@\[Phi]2[i],\[Phi]1[i]]HermitianCD[\[Nu],Bar@\[Phi]2[j],\[Phi]1[j]]]
		- c\[Phi]D2112[]HermitianCD[\[Nu],Bar@\[Phi]2[i],\[Phi]1[i]]HermitianCD[\[Nu],Bar@\[Phi]1[j],\[Phi]2[j]]
		- PlusHc[c\[Phi]D2111[]HermitianCD[\[Nu],Bar@\[Phi]2[i],\[Phi]1[i]]HermitianCD[\[Nu],Bar@\[Phi]1[j],\[Phi]1[j]]]
		- PlusHc[c\[Phi]D2122[]HermitianCD[\[Nu],Bar@\[Phi]2[i],\[Phi]1[i]]HermitianCD[\[Nu],Bar@\[Phi]2[j],\[Phi]2[j]]];		
		
	Lag6\[Phi]6 = c\[Phi]111111[] (Bar@\[Phi]1[i]\[Phi]1[i]) (Bar@\[Phi]1[j]\[Phi]1[j]) (Bar@\[Phi]1[k]\[Phi]1[k])
		+ c\[Phi]111122[] (Bar@\[Phi]1[i]\[Phi]1[i]) (Bar@\[Phi]1[j]\[Phi]1[j]) (Bar@\[Phi]2[k]\[Phi]2[k])
		+ c\[Phi]112222[] (Bar@\[Phi]1[i]\[Phi]1[i]) (Bar@\[Phi]2[j]\[Phi]2[j]) (Bar@\[Phi]2[k]\[Phi]2[k])
		+ PlusHc[c\[Phi]111121[] (Bar@\[Phi]1[i]\[Phi]1[i]) (Bar@\[Phi]1[j]\[Phi]1[j]) (Bar@\[Phi]2[k]\[Phi]1[k])]
		+ PlusHc[c\[Phi]222221[] (Bar@\[Phi]2[i]\[Phi]2[i]) (Bar@\[Phi]2[j]\[Phi]2[j]) (Bar@\[Phi]2[k]\[Phi]1[k])]
		+ c\[Phi]222222[] (Bar@\[Phi]2[i]\[Phi]2[i]) (Bar@\[Phi]2[j]\[Phi]2[j]) (Bar@\[Phi]2[k]\[Phi]2[k])
		+ PlusHc[c\[Phi]112121[] (Bar@\[Phi]1[i]\[Phi]1[i]) (Bar@\[Phi]2[j]\[Phi]1[j]) (Bar@\[Phi]2[k]\[Phi]1[k])]
		+ c\[Phi]112112[] (Bar@\[Phi]1[i]\[Phi]1[i]) (Bar@\[Phi]2[j]\[Phi]1[j]) (Bar@\[Phi]1[k]\[Phi]2[k])
		+ PlusHc[c\[Phi]222121[] (Bar@\[Phi]2[i]\[Phi]2[i]) (Bar@\[Phi]2[j]\[Phi]1[j]) (Bar@\[Phi]2[k]\[Phi]1[k])]
		+ c\[Phi]222112[] (Bar@\[Phi]2[i]\[Phi]2[i]) (Bar@\[Phi]2[j]\[Phi]1[j]) (Bar@\[Phi]1[k]\[Phi]2[k])
		+ PlusHc[c\[Phi]212121[] (Bar@\[Phi]2[i]\[Phi]1[i]) (Bar@\[Phi]2[j]\[Phi]1[j]) (Bar@\[Phi]2[k]\[Phi]1[k])]
		+ PlusHc[c\[Phi]212112[] (Bar@\[Phi]2[i]\[Phi]1[i]) (Bar@\[Phi]2[j]\[Phi]1[j]) (Bar@\[Phi]1[k]\[Phi]2[k])]
		+ PlusHc[c\[Phi]112221[] (Bar@\[Phi]1[i]\[Phi]1[i]) (Bar@\[Phi]2[j]\[Phi]2[j]) (Bar@\[Phi]2[k]\[Phi]1[k])];
		
	Lag6\[Psi]4 = cll[p,r,s,t] Bar@l[i,p]**\[Gamma][\[Nu]]**l[i,r] Bar@l[j,s]**\[Gamma][\[Nu]]**l[j,t]
		+ cqq1[p,r,s,t] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Alpha],i,r] Bar@q[\[Beta],j,s]**\[Gamma][\[Nu]]**q[\[Beta],j,t]
		+ cqq3[p,r,s,t] \[Tau]SU2L[J,i,j]Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Alpha],j,r] \[Tau]SU2L[J,k,m]Bar@q[\[Beta],k,s]**\[Gamma][\[Nu]]**q[\[Beta],m,t]
		+ clq1[p,r,s,t] Bar@l[i,p]**\[Gamma][\[Nu]]**l[i,r] Bar@q[\[Alpha],j,s]**\[Gamma][\[Nu]]**q[\[Alpha],j,t]
		+ clq3[p,r,s,t] \[Tau]SU2L[J,i,j]Bar@l[i,p]**\[Gamma][\[Nu]]**l[j,r] \[Tau]SU2L[J,k,m]Bar@q[\[Alpha],k,s]**\[Gamma][\[Nu]]**q[\[Alpha],m,t]
		+ cee[p,r,s,t] Bar@e[p]**\[Gamma][\[Nu]]**e[r] Bar@e[s]**\[Gamma][\[Nu]]**e[t]
		+ cuu[p,r,s,t] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**u[\[Alpha],r] Bar@u[\[Beta],s]**\[Gamma][\[Nu]]**u[\[Beta],t]
		+ cdd[p,r,s,t] Bar@d[\[Alpha],p]**\[Gamma][\[Nu]]**d[\[Alpha],r] Bar@d[\[Beta],s]**\[Gamma][\[Nu]]**d[\[Beta],t]
		+ ceu[p,r,s,t] Bar@e[p]**\[Gamma][\[Nu]]**e[r] Bar@u[\[Alpha],s]**\[Gamma][\[Nu]]**u[\[Alpha],t]
		+ ced[p,r,s,t] Bar@e[p]**\[Gamma][\[Nu]]**e[r] Bar@d[\[Alpha],s]**\[Gamma][\[Nu]]**d[\[Alpha],t]
		+ cud1[p,r,s,t] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**u[\[Alpha],r] Bar@d[\[Beta],s]**\[Gamma][\[Nu]]**d[\[Beta],t]
		+ cud8[p,r,s,t] TSU3c[A,\[Alpha],\[Beta]] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**u[\[Beta],r] TSU3c[A,\[Delta],\[Kappa]] Bar@d[\[Delta],s]**\[Gamma][\[Nu]]**d[\[Kappa],t]
		+ cle[p,r,s,t] Bar@l[i,p]**\[Gamma][\[Nu]]**l[i,r] Bar@e[s]**\[Gamma][\[Nu]]**e[t]
		+ clu[p,r,s,t] Bar@l[i,p]**\[Gamma][\[Nu]]**l[i,r] Bar@u[\[Alpha],s]**\[Gamma][\[Nu]]**u[\[Alpha],t]
		+ cld[p,r,s,t] Bar@l[i,p]**\[Gamma][\[Nu]]**l[i,r] Bar@d[\[Alpha],s]**\[Gamma][\[Nu]]**d[\[Alpha],t]
		+ cqe[p,r,s,t] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Alpha],i,r] Bar@e[s]**\[Gamma][\[Nu]]**e[t]
		+ cqu1[p,r,s,t] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Alpha],i,r] Bar@u[\[Beta],s]**\[Gamma][\[Nu]]**u[\[Beta],t]
		+ cqu8[p,r,s,t] TSU3c[A,\[Alpha],\[Beta]] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Beta],i,r] TSU3c[A,\[Delta],\[Kappa]] Bar@u[\[Delta],s]**\[Gamma][\[Nu]]**u[\[Kappa],t]
		+ cqd1[p,r,s,t] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Alpha],i,r] Bar@d[\[Beta],s]**\[Gamma][\[Nu]]**d[\[Beta],t]
		+ cqd8[p,r,s,t] TSU3c[A,\[Alpha],\[Beta]] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Beta],i,r] TSU3c[A,\[Delta],\[Kappa]] Bar@d[\[Delta],s]**\[Gamma][\[Nu]]**d[\[Kappa],t]
		+ PlusHc[ 
			+ cledq[p,r,s,t] Bar@l[i,p]**e[r] Bar@d[\[Alpha],s]**q[\[Alpha],i,t]
			+ cquqd1[p,r,s,t] Bar@q[\[Alpha],i,p]**u[\[Alpha],r] \[Epsilon]SU2L[i,j] Bar@q[\[Beta],j,s]**d[\[Beta],t]
			+ cquqd8[p,r,s,t] TSU3c[A,\[Alpha],\[Beta]] Bar@q[\[Alpha],i,p]**u[\[Beta],r] \[Epsilon]SU2L[i,j] TSU3c[A,\[Delta],\[Kappa]] Bar@q[\[Delta],j,s]**d[\[Kappa],t]
			+ clequ1[p,r,s,t] Bar@l[i,p]**e[r] \[Epsilon]SU2L[i,j] Bar@q[\[Alpha],j,s]**u[\[Alpha],t]
			+ clequ3[p,r,s,t] Bar@l[i,p]**\[Sigma][\[Nu],\[Rho]]**e[r] \[Epsilon]SU2L[i,j] Bar@q[\[Alpha],j,s]**\[Sigma][\[Nu],\[Rho]]**u[\[Alpha],t]
		];
		
	Lag6BNV = cduq[p,r,s,t] Bar@\[Epsilon]SU3c[\[Alpha],\[Beta],\[Delta]]Bar@\[Epsilon]SU2L[i,j] CConj@Bar@d[\[Alpha],p]**u[\[Beta],r] CConj@Bar@q[\[Delta],i,s]**l[j,t]
		+ cqqu[p,r,s,t] Bar@\[Epsilon]SU3c[\[Alpha],\[Beta],\[Delta]]Bar@\[Epsilon]SU2L[i,j] CConj@Bar@q[\[Alpha],i,p]**q[\[Beta],j,r] CConj@Bar@u[\[Delta],s]**e[t]
		+ cqqq[p,r,s,t] Bar@\[Epsilon]SU3c[\[Alpha],\[Beta],\[Delta]]Bar@\[Epsilon]SU2L[i,j] Bar@\[Epsilon]SU2L[k,m] CConj@Bar@q[\[Alpha],i,p]**q[\[Beta],k,r] CConj@Bar@q[\[Delta],m,s]**l[j,t]
		+ cduu[p,r,s,t] Bar@\[Epsilon]SU3c[\[Alpha],\[Beta],\[Delta]] CConj@Bar@d[\[Alpha],p]**u[\[Beta],r] CConj@Bar@u[\[Delta],s]**e[t];

	Lag= PlusHc[Lag5] + PlusHc[Lag6\[Psi]2\[Phi]3] + PlusHc[Lag6\[Psi]2X\[Phi]] + Lag6\[Psi]2\[Phi]2D + Lag6X3 + Lag6X2\[Phi]2 + Lag6\[Phi]4D2 + Lag6\[Phi]6 + Lag6\[Psi]4 + PlusHc[Lag6BNV];
	
	RelabelIndices@Lag
]
