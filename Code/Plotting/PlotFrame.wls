#!/usr/bin/env wolframscript
(* ::Package:: *)
(* SetDirectory["/home/dino/Code/DinoFastMD/sim-20-04-16-13:08:06/"]; *)
string=$ScriptCommandLine[[2]];
colstring=If[Length[$ScriptCommandLine]>2,$ScriptCommandLine[[3]],"void"];
(* PeriodicAddVector[x_, y_, ll_] := Module[{dx, dx2},dx = x + y;dx2 = Map[If[# > ll, # - ll, If[# < 0, # + ll, #]] &, dx];dx2]
rho0=5;
n=1000; *)
Print[$ScriptCommandLine];

x = Import[string];
d = Length[x[[1]]];
l = Length[x];
If[l==0,Exit[1]];
col=If[Length[$ScriptCommandLine]>2,Map[RGBColor,Import[colstring]],Table[Blue,{i,1,l}]];

s = Max[Flatten[x]];

PlotFrame[pos_,sm_,dm_] := Module[{},
Return[
Switch[
dm,
2,
Graphics[Map[Disk[#,0.5]&,pos],PlotRange->{{0,sm},{0,sm}}],
3,
Graphics3D[Transpose[{col,MapIndexed[Sphere[#1, If[First[#2] < 1000, 0.5, 0.25]] &, pos]}],PlotRange->{{0,sm},{0,sm},{0,sm}},ImageSize->1000],
_,
Graphics[Disk[]]
]]
];

Print["defined everything"];

img = PlotFrame[x,77.8181,d];

Print["Plotted Image"];

dir = DirectoryName[string];

str = FileBaseName[string];



Export[dir<>str<>".jpg",img]
Print["exported successfully"];


Exit[]
