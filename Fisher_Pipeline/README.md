# Fisher Pipeline
command:
<BR> python3 Fisher_Pipeline.py --p [path to diretory of bamfiles] 
<br>                       --ctrl [path to txt file of control sample names]
<br>                       --test [path to txt file of test sample names] 
<br>                       --dust [threshold for dust score]
<br>                       --t [threshold for filtering control]


## Instruction to Run
  1. To run the script, make sure to have a directory with all the bamfiles to analyze. 
  2. Make sure to have a txt file with all the control sample names and another txt file with all the test sample names. Names 
     should be separated by a line within each file. Make sure the sample names match the bamfile names. For instance,
     S01 sample name in the txt file should correspond to S01.bam 
     
## Arguments

## Output
  1. A binary heatmap of all the significant RNAs
  2. A normalized expression heatmap of all the significant RNAs
  3. A bed file for the significant RNAs
