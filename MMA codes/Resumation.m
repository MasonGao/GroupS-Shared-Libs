(* ::Package:: *)

(* ::Subsection:: *)
(*Anomalous Dimension input*)


betare={\[Beta]0->11/3 Ca-4/3 Tf nf  ,
\[Beta]1->34/3 Ca^2-20/3 Ca Tf nf-4Cf Tf nf ,
\[Beta]2->2857/54 Ca^3+(2Cf^2-205/9 Cf Ca-1415/27 Ca^2)Tf nf+(44/9 Cf+158/27 Ca)Tf^2 nf^2,
\[Beta]3->149753/6+3564 Zeta[3]-(1078361/162+6508/27 Zeta[3])nf+(50065/162+6472/81 Zeta[3])nf^2+1093/729 nf^3};


cuspre={\[Gamma]cusp0->4,
\[Gamma]cusp1->4((67/9-\[Pi]^2/3)Ca-20/9 Tf nf ) ,
\[Gamma]cusp2->4(Ca^2(245/6-(134\[Pi]^2)/27+(11\[Pi]^4)/45+22/3 Zeta[3])+Ca Tf nf(-(418/27)+(40\[Pi]^2)/27-56/3 Zeta[3])+Cf Tf nf(-(55/3)+16Zeta[3])-16/27 Tf^2nf^2)};


ADimReABC=Join[betare,cuspre];


adimre=Join[{Cf->4/3,nf->5,Tf->1/2,Ca->3},ADimReABC];


\[Gamma]cexpand={\[Gamma]cusp->(\[Alpha]/(4\[Pi]))\[Gamma]cusp0+(\[Alpha]/(4\[Pi]))^2 \[Gamma]cusp1+(\[Alpha]/(4\[Pi]))^3 \[Gamma]cusp2};


\[Gamma]qexpand={\[Gamma]q->(\[Alpha]/(4\[Pi]))\[Gamma]q0+(\[Alpha]/(4\[Pi]))^2 \[Gamma]q1};


\[Gamma]gexpand={\[Gamma]g->(\[Alpha]/(4\[Pi]))\[Gamma]g0+(\[Alpha]/(4\[Pi]))^2 \[Gamma]g1};


\[Gamma]qgre={\[Gamma]q0->-3Cf,\[Gamma]q1->Cf^2 (-(3/2)+2\[Pi]^2-24\[Zeta]3)+Cf Ca (-961/54-11 \[Pi]^2/6+26\[Zeta]3)+Cf Tf nf (130/27+2 \[Pi]^2/3),\[Gamma]g0->-\[Beta]0,\[Gamma]g1->Ca^2 (-(692/27)+11/18 \[Pi]^2+2\[Zeta]3)+Ca Tf nf (256/27-2/9 \[Pi]^2)+4Cf Tf nf};


rre={r[\[Nu]_,\[Mu]_]:>1+(\[Beta]0 Log[\[Nu]/\[Mu]] \[Alpha][\[Mu]])/(2 \[Pi])+(\[Beta]1 Log[\[Nu]/\[Mu]] \[Alpha][\[Mu]]^2)/(8 \[Pi]^2)};


(* ::Subsection:: *)
(*SudakovS,  a\[Gamma] *)


OrderVec[o_]:=Switch[o,-1,{1,0,0,0},0,{1,1,0,0},1,{1,1,1,0},2,{1,1,1,1}];


SudakovSO[\[Nu]_,\[Mu]_]:=\[CapitalGamma]0/(4\[Beta]0^2) {((4\[Pi])/\[Alpha][\[Nu]] (1-1/r[\[Nu],\[Mu]]-Log[r[\[Nu],\[Mu]]])),(\[CapitalGamma]1/\[CapitalGamma]0-\[Beta]1/\[Beta]0)(1-r[\[Nu],\[Mu]]+Log[r[\[Nu],\[Mu]]])+\[Beta]1/(2\[Beta]0) Log[r[\[Nu],\[Mu]]]^2,\[Alpha][\[Nu]]/(4\[Pi]) (((\[Beta]1 \[CapitalGamma]1)/(\[Beta]0 \[CapitalGamma]0)-\[Beta]2/\[Beta]0)(1-r[\[Nu],\[Mu]]+r[\[Nu],\[Mu]]Log[r[\[Nu],\[Mu]]])+(\[Beta]1^2/\[Beta]0^2-\[Beta]2/\[Beta]0)(1-r[\[Nu],\[Mu]])Log[r[\[Nu],\[Mu]]]-(\[Beta]1^2/\[Beta]0^2-\[Beta]2/\[Beta]0-(\[Beta]1 \[CapitalGamma]1)/(\[Beta]0 \[CapitalGamma]0)+\[CapitalGamma]2/\[CapitalGamma]0) (1-r[\[Nu],\[Mu]])^2/2),(\[Alpha][\[Nu]]/(4\[Pi]))^2 ( 
((\[Beta]1 \[Beta]2)/\[Beta]0^2-\[Beta]1^3/(2\[Beta]0^3)-\[Beta]3/(2\[Beta]0)+\[Beta]1/\[Beta]0 (\[CapitalGamma]2/\[CapitalGamma]0-\[Beta]2/\[Beta]0+\[Beta]1^2/\[Beta]0^2-(\[Beta]1 \[CapitalGamma]1)/(\[Beta]0 \[CapitalGamma]0)) r[\[Nu],\[Mu]]^2/2)Log[r[\[Nu],\[Mu]]]
+(\[CapitalGamma]3/\[CapitalGamma]0-\[Beta]3/\[Beta]0+(2\[Beta]1 \[Beta]2)/\[Beta]0^2+\[Beta]1^2/\[Beta]0^2 (\[CapitalGamma]1/\[CapitalGamma]0-\[Beta]1/\[Beta]0)-(\[Beta]2 \[CapitalGamma]1)/(\[Beta]0 \[CapitalGamma]0)-(\[Beta]1 \[CapitalGamma]2)/(\[Beta]0 \[CapitalGamma]0)) (1-r[\[Nu],\[Mu]])^2/3
+((3\[Beta]3)/(4\[Beta]0)-\[CapitalGamma]3/(2\[CapitalGamma]0)+\[Beta]1^3/\[Beta]0^3-(3\[Beta]1^2 \[CapitalGamma]1)/(4\[Beta]0^2 \[CapitalGamma]0)+(\[Beta]2 \[CapitalGamma]1)/(\[Beta]0 \[CapitalGamma]0)+(\[Beta]1 \[CapitalGamma]2)/(4\[Beta]0 \[CapitalGamma]0)-(7\[Beta]1 \[Beta]2)/(4\[Beta]0^2))(1-r[\[Nu],\[Mu]])^2
+((\[Beta]1 \[Beta]2)/\[Beta]0^2-\[Beta]3/\[Beta]0-(\[Beta]1^2 \[CapitalGamma]1)/(\[Beta]0^2 \[CapitalGamma]0)+(\[Beta]1 \[CapitalGamma]2)/(\[Beta]0 \[CapitalGamma]0)) (1-r[\[Nu],\[Mu]])/2)}/.{\[CapitalGamma]0->\[Gamma]cusp0,\[CapitalGamma]1->\[Gamma]cusp1,\[CapitalGamma]2->\[Gamma]cusp2};a\[CapitalGamma]O[\[Nu]_,\[Mu]_]:=\[CapitalGamma]0/(2\[Beta]0) {0,Log[r[\[Nu],\[Mu]]],\[Alpha][\[Nu]](\[CapitalGamma]1/\[CapitalGamma]0-\[Beta]1/\[Beta]0) (r[\[Nu],\[Mu]]-1)/(4\[Pi]),(\[CapitalGamma]2/\[CapitalGamma]0-\[Beta]2/\[Beta]1-\[Beta]1/\[Beta]0 (\[CapitalGamma]1/\[CapitalGamma]0-\[Beta]1/\[Beta]0)) \[Alpha][\[Nu]]^2 (r[\[Nu],\[Mu]]^2-1)/(32\[Pi]^2)}/.{\[CapitalGamma]0->\[Gamma]cusp0,\[CapitalGamma]1->\[Gamma]cusp1,\[CapitalGamma]2->\[Gamma]cusp2};a\[Gamma]O[\[Nu]_,\[Mu]_]:=\[CapitalGamma]0/(2\[Beta]0) {0,Log[r[\[Nu],\[Mu]]],(\[CapitalGamma]1/\[CapitalGamma]0-\[Beta]1/\[Beta]0) \[Alpha][\[Nu]](r[\[Nu],\[Mu]]-1)/(4\[Pi]),(\[CapitalGamma]2/\[CapitalGamma]0-\[Beta]2/\[Beta]1-\[Beta]1/\[Beta]0 (\[CapitalGamma]1/\[CapitalGamma]0-\[Beta]1/\[Beta]0)) \[Alpha][\[Nu]]^2 (r[\[Nu],\[Mu]]^2-1)/(32\[Pi]^2)}/.{\[CapitalGamma]0->\[Gamma]0,\[CapitalGamma]1->\[Gamma]1};


SudakovS[\[Nu]_,\[Mu]_,o_]:=Dot[SudakovSO[\[Nu],\[Mu]],OrderVec[o]];
a\[CapitalGamma][\[Nu]_,\[Mu]_,o_]:=Dot[a\[CapitalGamma]O[\[Nu],\[Mu]],OrderVec[o]];
a\[Gamma][\[Nu]_,\[Mu]_,o_]:=Dot[a\[Gamma]O[\[Nu],\[Mu]],OrderVec[o]];
