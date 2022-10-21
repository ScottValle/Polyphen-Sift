# Polyphen-Sift
A data repository for the paper "Comparing and evaluating the performance of pathogenicity predicting algorithms PolyPhen-2 and SIFT"

The original dataset:
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_1.0/2014/clinvar_20141202.vcf.gz

The following filtered datasets can be found in the Dataset directory:
  -BLOSUM62.txt
    The BLOSUM score for all possible amino acid mutations
  -HGVS_2014_VEP_baseline.tsv
    The HGVS ID of 400 mutations, the resulting change in amino acids and in nucleotides.
  -HGVS_2014_benchmark.tsv
    The HGVS ID of 2982 mutations and their impact
  -HGVS_2014_polyphen_scores.tsv
    The HGVS ID of 400 mutations and their polyphen prediction score
  -HGVS_2014_sift_scores.tsv
    The HGVS ID of 400 mutations and their SIFT prediction score

The following scripts can be found in the Scripts directory:
  -skeleton_script_baseline_model.py
  
    Baseline impact predictor of SNPs in VEP format. Uses raw BLOSUM62 matrix from a text file for scoring.
    usage: skeleton_script_baseline_model.py [-h] -o OUT_PATH vep blosum

    Generates baseline model predictions.

    positional arguments:
      vep          a path to the VEP input file
      blosum       a path to the BLOSUM62 input file

    optional arguments:
      -h, --help   show this help message and exit
      -o OUT_PATH  a path to write the output .tsv file with baseline model
                   scores. This arguments is required!

  -skeleton_script_create_roc_plot.py
  
    This script draws ROC plot with/without gradient color and calculates the AUC score.
    Code partially adapted from https://gist.github.com/podshumok/c1d1c9394335d86255b8
    It should be executed by specifying input files and output file path e.g.:
    python3 skeleton_script_create_roc_plot.py -ibench <benchmark_filepath> -ipred <predictor_filepath -o <output_filepath> -color
    -color is an optional argument: ROC plot for a single predictor (SIFT, PolyPhen-2, or BLOSUM62) will show
    threshold scores with gradient color. If the argument is absent, gradient color will not be shown.
    To plot ROC curves for all three predictors in one figure (without gradient color):
    python3 skeleton_script_create_roc_plot.py -ibench <benchmark_filepath> -ipred <sift_scores_filepath> -ipred <polyphen_scores_filepath> -ipred <baseline_scores_filepath>  -o <ROCplot_output_filepath>
    
  -skeleton_script_roc_plot_tsv.py
  
     Draws a ROC plot from three .tsv files.
    It should be executed by specifying three input predictor .tsv files from one of the ClinVar datasets
    with HGVS IDs and scores (three predictors: SIFT, PolyPhen, and BLOSUM62).
    -ititle or --plot_title (optional): with this argument followed by a title string, a title provided by the user
    will be added to plot. Otherwise, the default will be used.
    Execute with:
    python3 skeleton_script_plot_tsv.py -itsv <ROCplot_filename_sift>_xy.tsv -itsv <ROCplot_filename_polyphen>_xy.tsv -itsv <ROCplot_filename_baseline>_xy.tsv -o <ROCplot_filename_all>_xy.tsv>.png
