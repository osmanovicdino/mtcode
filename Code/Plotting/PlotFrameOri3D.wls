#!/usr/bin/env wolframscript
string1=$ScriptCommandLine[[2]];

string2=$ScriptCommandLine[[3]];

x = Import[string1];
o = Import[string2]
d = Length[x[[1]]];
l = Length[x];
If[l==0,Exit[1]];

s = Max[Flatten[x]];

PlotFrameOri[strp_,stro_,vecs_,sm_]:=Module[{pos,ori,sur,pars},
pos=Import[strp];
ori=Import[stro];
ori=Map[Transpose[Partition[#,3]]&,ori];
sur=Flatten[Table[{Red,Sphere[pos[[i]]+1/2 ori[[i]].vecs[[j]],0.1]},{i,1,Length[pos]},{j,1,Length[vecs]}],1];
pars=Map[{Blue,Sphere[#,0.4]}&,pos];
Graphics3D[{pars,sur},PlotRange->{{0,sm},{0,sm},{0,sm}}]
]


vecs = {{1,0,0},{-1,0,0}};

img =  PlotFrameOri[string1,string2,vecs,s];

dir = DirectoryName[string1];

str = FileBaseName[string1];

Export[dir<>str<>".x3d",img]