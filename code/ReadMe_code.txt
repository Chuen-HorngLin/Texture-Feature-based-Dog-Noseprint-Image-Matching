
*** Main code for this system

1. ImageYoloCut.py: Cuts the dog noses from original photos (subprograms: all.bat, half.bat, image_scale.exe).

2. dognose_segmentation.m: Extracts dog nose pattern images.

3. dognose_DB: Establishes a database of scale-like block features.

4. dognose_match.m: Matches dog nose images.

5. show_result.m: Visualizes the comparison of scale-like block patterns.


=========================================
*** Training scale-like block weights

1. get_match_db.m: Constructs a dataset for comparing scale-like block features, used for Genetic Algorithm (GA) to find the optimal parameters.

2. dognose_match_GA.m: Genetic algorithm for optimizing parameters in dog nose pattern matching (subprogram: ga_fun.m).
