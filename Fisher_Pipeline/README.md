# Fisher Pipeline
Command:
<br> `python3 Fisher_Pipeline.py` 
<br>      `--p [path to diretory of bamfiles]`
<br>  `--ctrl [path to txt file of control sample names]`
<br>  `--test [path to txt file of test sample names]` 
<br>  `--dust [threshold for dust score, default 2.5]`
<br>   `--t [threshold for filtering control, default 3]`


## Instruction to Run
  1. To run the script, make sure to have a directory with all the bamfiles to analyze. 
  2. Make sure to have a txt file with all the control sample names and another txt file with all the test sample names. Names 
     should be separated by a line within each file. Make sure the sample names match the bamfile names. For instance,
     S01 sample name in the txt file should correspond to S01.bam 
     See `ex_control.txt` and `ex_test.txt` for example.
  3. Example call:
      <br> `python3 Fisher_Pipeline.py --p ex_data/ --ctrl ex_control.txt --test ex_test.txt`
     
## Arguments
   `--dust`: Set threshold for dust score to eliminate low-complexity sequences. Pipeline will exclude RNAs that have a dust score higher or equal to [dust]
   <br>`--t`: Set threshold for filtering controls. Pipeline will exclude RNAs that are present in at least [t] amount of the control samples. 
   <br> In addition, pipeline will automatically filter out all RNAs that are present in only one sample.
## Output
  1. A binary heatmap of all the significant RNAs
  2. A normalized expression heatmap of all the significant RNAs
  3. A bed file of the significant RNAs
