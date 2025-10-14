# TreeUnimodalityPatternBoost
Pattern boost code used to build counter-examples to the log concavity conjecture for trees.

This repo contains all of the code used during our exploration of the unimodality conjecture for trees. The only code that was written by the authors of this paper are the julia files beginning "problem_erdos", and the python file "PruferToTree." All of the other python and julia files, which can be found in the folder "PatternBoost", were written by the authors of the original PatternBoost paper - François Charton, Jordan S. Ellenberg, Adam Zsolt Wagner, and Geordie Williamson. In certain cases we needed to very lightly modify these files for our purposes, but the vast majority of that code can be found on its original github repo: https://github.com/zawagner22/transformers_math_experiments.

To run this, first edit the top lines in search_fc.jl to select the local optimization that will be used. The options are:

1. problem_erdos_6.0.jl - This is the most standard version.

2. problem_erdos_6.0Big.jl - This is the version you will need to use if the number of vertices is above 68. It uses BigInt to avoid overflow. **This is the only version with this particular protection on the coefficients of the independence polynomial**. All other versions, however, have overflow protection on the score computation.

3. problem_erdos_NoLocal.jl - This is the version which skips local optimization entirely, and only relies on the transformer teaching itself.

4. problem_erdos_NoPath.jl - This is the version which punishes the machine if it produces the path graph. It seems necessary to do this if you are look for breakage at an index away from half the number of vertices.

Once you've selected the local optimization file, be sure to go into the file to change whatever parameters you need to. This includes the length of the prufer codes that are being used (i.e. two less than the number of vertices), where you are looking for log concavity breakage, as well as the total number of swaps to perform and how non-edges are added during optimization. **IMPORTANTLY, the parameter N in all of these files will be the length of the prufer code being considered, NOT the number of vertices.** In other words, N should be two less than the the number of vertices of the tree.

Once all of this has been settled, you just have to run fc_loop.py with your selected parameters. We recommend looking through the code to see all parameters that you can tinker with. In my experience, the following works perfectly well on a GPU:

python3 fc_loop.py --sample-only 100000 --max-steps 8000 --max_epochs 5 --n_tokens 60  --dump_path experiment_output --exp_name test_run

**Note that you must always set n_tokens to the number of vertices in the trees you are exploring.** If you are trying to run this on a local machine, the numerical parameters should be made smaller, and you must add the line --cpu true For instance,

python3 fc_loop.py --cpu true --sample-only 2000 --max-steps 1000 --max_epochs 3 --n_tokens 60 --dump_path experiment_output --exp_name test_run

The file PruferToTree.ipyn is SageMath code that is unnecessary for all of the above computation. It allows you to input any Prufer code, and it will output a drawing of the resulting tree. The file Problem_Erdos_IP.jl is Julia code which allows you to input a prufer code, and it will output its indepednence polynomial, as well as the value of the scoring function at every index of the independence sequence.

Finally, the folder 60_vertex_output contains all of the output files from our application of these methods to the case of trees with 60 vertices.
