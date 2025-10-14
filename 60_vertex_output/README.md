This folder contains the output of our methods to the 60 vertex case. This experiment used problem_erdos_NoPath as its local search method, with 10 swaps, and randomization of added edges (i.e. order_edges = -1). The PatternBoost parameters can be seen in the terminal command present at the very top of train.log.

The most relevant files in this folder are:

- distribution.txt contains the distribution of scores among the final output of the transformer after optimization

- training_distribution.txt contains the distribution of the 50,000 best performing trees in the final output of the transformer after optimization.

- plot_i.png is a histogram for the full output of the (i-1)-st epoch after optimization.

- plot_training_i.png is a histograph for the top 50,000 performers of the (i-1)-st epoch after optimization.

- search_output_i.txt contains the prufer codes for the top 50,000 performers of the (i-1)-st epoch after optimization. Importantly, these prufer codes appear in decreasing order of score, so that one can align the trees of the final output (search_output_11) with training_distribution.txt

- search_output_i-tokenized.txt contains the same content as the previous file, but tokenized.

- train.log this is a log of terminal screen that ran with the experiment. One can see from this each individual epoch's distribution, as well as the time that the entire process took (just under 9 hours).


Other files in this folder which are less important are:

- transformer-output-decoded.txt This is the raw output of the final epoch of the transformer. All 150,000 samples in some order.

- out.txt is the final tokenized output of the transformer after the 10th, and final in this case, epoch.

- model.pt is the final transformer model used by PatternBoost.

- params.pkl is created to store python parameters.


To translate the prufer code outputs to actual trees, one can use the SageMath code PruferToTree.ipynb included in this git repo.

