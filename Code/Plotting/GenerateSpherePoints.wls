#!/usr/bin/env wolframscript

partN=ToExpression@($ScriptCommandLine[[2]]);
targetdensity=ToExpression@($ScriptCommandLine[[3]]);

distance=Compile[{{x, _Real, 1}, {y, _Real, 1}}, Sqrt[Total[(x - y)^2]]];

Rm=SolveValues[{partN/(4*3.1415* R^2) == 1.0, R > 0}, R][[1]];
sph = Rm*SpherePoints[partN];
sphd = Table[Table[distance[sph[[i]], sph[[j]]], {j, 1, partN}], {i,1, partN}];

(* asd = Table[Sort[sphd[[i]]][[2 ;; 7]], {i, 1, Length[sphd]}]; *)
allpairs2 = Flatten[Table[ord = sphd[[j]];
    (*as=Thread[j->Flatten@Table[Position[ord,i],{i,2,7}]];*)
    Map[{j, #} &, Flatten[Table[Position[ord, RankedMin[ord, i]], {i, 2, 7}]]]
    ,
    {j, 1, Length[sphd]}], 1];



nap2 = DeleteDuplicates@(Map[Sort, allpairs2]);

Export[Directory[]<>"/IsocohedronI.csv", nap2 - 1]

Export[Directory[]<>"/IsocohedronP.csv", sph]