Print["hello"];
f = FileNames["./pos*.csv"];
Print[Length[f]," files found"];
Do[img = Graphics3D[Transpose[{Map[RGBColor, Import["col.csv"]], MapIndexed[Sphere[#1, If[First[#2] < 1000, 0.5, 0.25]] &, Import[f[[i]]]]}], ImageSize -> 1000];
Export[FileBaseName[f[[i]]] <> ".jpg", img];
Clear[img];
Print[i], {i, 1, Length[f]}]