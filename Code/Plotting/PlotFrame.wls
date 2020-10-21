#!/usr/bin/env wolframscript
(* ::Package:: *)
(* SetDirectory["/home/dino/Code/DinoFastMD/sim-20-04-16-13:08:06/"]; *)
string=$ScriptCommandLine[[2]];

(* PeriodicAddVector[x_, y_, ll_] := Module[{dx, dx2},dx = x + y;dx2 = Map[If[# > ll, # - ll, If[# < 0, # + ll, #]] &, dx];dx2]
rho0=5;
n=1000; *)
x = Import[string];
d = Length[x[[1]]];
l = Length[x];
If[l==0,Exit[1]];

s = Max[Flatten[x]];

PlotFrame[pos_,sm_,dm_] := Module[{},
Return[
Switch[
dm,
2,
Graphics[Map[Disk[#,0.5]&,pos],PlotRange->{{0,sm},{0,sm}}],
3,
Graphics3D[Transpose[{Table[Blue,{i,1,l}],Map[Sphere[#,0.5]&,pos]}],PlotRange->{{0,sm},{0,sm},{0,sm}}],
_,
Graphics[Disk[]]
]]
]

img = PlotFrame[x,s,d];

dir = DirectoryName[string];

str = FileBaseName[string];



Export[dir<>str<>".jpg",img]
