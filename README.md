# absolute-model-data-fit



The repository contains Perl scripts for assessment of absolute fit of evolutionary models -\
test_stage1.pl and test_stage2.pl published in (place for citation).\
Please see the original article for detailed test description and examples.

The test compares observed data (multiple sequence alignment of DNA sequences, "OD") evolved\
under "empirical evolutionary model" (EM) to distributions of parametric replicates generated\
under explicitly postulated full evolutionary models (substitution model + tree combinations),\
termed "simulation models" (SMs).

If the observed alignment contains character states not predicted by DNA substitution models,\
these must be removed before the method is used for model evaluation. In this case, a necessary\
data preparation step might include, in any combination, removal of corresponding sequences and/or\
removal of sites. The preferred trimming scheme should be selected by an analyst prior to analysis.

Each SM should assume a fixed parameter combination and a fully resolved tree with fixed branch\
lengths. When OD contains biological sequences, all SM model parameters should assume values\
registered at the ML optimum under each evolutionary model specification (ML tree + ML substitution\
model parameters for the OD). The cladograms used to constrain ML searches should represent\
competing evolutionary hypotheses for the OD. These hypotheses might include conflicting trees\
obtained in exploratory phylogenetic analyzes under various models based on the OD or its\
partitions. Published hypotheses concerning phylogenetic relationship of taxa in the OD\
can also be included into analysis. In this case one needs to specify each such hypothesis as a set\
of full/partial topological constraints prior to ML search.

When the final set of alternative tree topologies (termed "T set") is assembled, series of ML\
searches should be run based on the observed data, under each distinct full topological constraint\
and optimal substitution model to obtain SM specifications. An optimal substitution model component\
of SMs can be selected using various published criteria. Alternatively, it can also be selected\
using the scripts in the repository. In this case, several sets of ML searches, each assuming a\
distinct sibstitution model (e.g GTR+G,GTR+G+I, etc.) constrained with each constraint from the T set,\
should be conducted. Replicate files containing a fixed number of replicates of the same length as the OD\
should be generated under each SM.

OD and the replicates should be in fasta or sequential phylip format (1 sequence per line).\
A test file with a list of filenames of T set constraints, in newick format, one tree file name per\
line, should be prepared too.

When the above prerequisites are met, the first stage of the test, performed by test_stage1.pl\
script, can be conducted.

If called without arguments, the script outputs instructions for program execution:

test_stage1.pl (v. 1.0) script calculates:

-TTC value                          ("Tvalue_TTC" in the output file),\
-Default GGg-based component of TC1 ("Tvalue_GGg" in the output file),\
-Full Gelfand-Ghosh statistic       ("Tvalue_GGf" in the output file)\
for the test of absolute model-data fit (Goremykin (2019)).

command line synopsis:

-i name of a file containing observed alignment(s) in fasta or sequential phylip format (1 sequence per line)\
-r name of a file containing parametric replicates in fasta or sequential phylip format (1 sequence per line)\
-t list of filenames of model tree files, in newick format, one tree file name per line\
-o output file name. If subsequent analyzes with the test_stage2.pl script should be conducted, please make\
   the output file names compatible with the format that the test_stage2.pl script requires (described when\
   the test_stage2.pl script is called without arguments)\
-s (optional) alphabet type\
   "-s dna":  DNA alphabet\
   "-s prot": Protein alphabet\
   "-s ry":   alphabet which has R and Y characters only\
   Note: The performance of the test has been tested with DNA data alphabet only. Other options are added for\
   exploratory purposes. If the alphabet type is not specified, the script assumes the DNA alphabet by default.\
-p (optional) sets a P factor value (default 10,000)\
-g (optional) enables calculation of full Gelfand and Ghosh statistic
__________

Output format of test_stage1.pl script

in the case of EM representation by alignment of biological sequences, each output file has the following format:

Tvalue_GGf: 959.309695708264\
Tvalue_GGg: 294.473066605363\
Tvalue_TTC: 0.0210386510130374

wherein:\
Tvalue_GGf is full Gelfand and Ghosh statistic (which can be optionally calculated by test_stage1.pl).\
Tvalue_GGg is GGg statistic related to substitution model fit which is a component of Gelfand and Ghosh statistics.\
Tvalue_TTC is the statistic for the test of tree component (TTC)

In the case of EM representation by a distribution of replicates, each output file has the following format:

Tvalue_GGf: 959.309695708264 782.938618148101\
Tvalue_GGg: 294.473066605363 118.101989045201\
Tvalue_TTC: 0.0210386510130374 0.0263469766498732

The leftmost values represent above statistics calculated based on replicate 1, the next values\
represent above statistics calculated based on replicate 2, etc.
______________________________________________________________________________________________________

The second stage of the test is performed by test_stage2.pl script, which should be run\
in the directory containing the files produced by test_stage1.pl script.

If called without arguments, the script outputs instructions for program execution:

test_stage2.pl (v. 1.0) script calculates the test values based on the data produced by test_stage1.pl\
script and provides a summary of the test results.

synopsis: test_stage2.pl argument1 argument2 > result

argument1: a file produced using the following command line:\
"grep Tvalue *.extension > file", wherein *.extension is a common extension of the\
output files produced by the test_stage1.pl script

argument2 (optional): -g\
Specifies the full Gelfand-Ghosh (GG = GGg + GGp) estimator of substitution model fit instead of\
the default GGg function to be used in calculation of TSC. If GG statistic was not calculated by\
test_stage1.pl script, the option should not be used.

Note 1: The script requires test_stage1.pl output file names to be in a certain format.\
If the empirical model is represented by a biological dataset then the file names should be in\
the following format: *simulation_model[MODELNAME]*.extension\
If the empirical model(s) is/are represented by a distribution of replicates then the file names should\
be in the following format: *simulation_model[MODELNAME1].vs.empirical_model[MODELNAME2]*.extension\
wherein "*" is a wildcard character. MODELNAME and MODELNAME1 designate the evolutionary model\
which fit to the observed data should be estimated (SM). MODELNAME2 designate the evolutionary model\
used to generate the observed data (EM).

Each MODELNAME, MODELNAME1 and MODELNAME2 should be 8 characters in length, of which the first\
two are digits (tree code), designating model tree topology and the rest 6 characters designate\
the substitution model used to generate the corresponding replicates, for example: 10GTR+GI

Note 2: In the comparisons involving representation of the observed data by replicates:
- The set of model tree topologies (T set) should be the same for all sets of SMs assuming different 
  substitution model specifications. 
- The set of SMs which shares substitution model specification with EM(s) should be included into analysis. 
- All EMs should share a distinct substitution model component (e.g. GTR+GI). 
- Each of tree codes should designate the same, distinct tree topology. 
- EM model tree topologies should be present in the T set.

__________

Output format of test_stage2.pl script:

If the empirical model (EM) is represented by one multiple sequence alignment, the test_stage2.pl outputs\
a table containing SM-specific T values sorted in ascending order and corresponding SM designations.\
The best fitting SM is the one with the lowest T value.\
An example of the table is provided below:

Models sorted in ascending order of test values:

model: 10GTR__G test value: 7.00545534776865\
model: 16GTR__G test value: 28.7464343080381\
model: 11GTR__G test value: 49.2098367539291\
model: 12GTR__G test value: 103.722078030946\
model: 13GTR__G test value: 107.448315458472\
model: 14GTR__G test value: 107.644576777394\
model: 15GTR__G test value: 114.689475486713\
model: 10HKY__G test value: 120.975080017367\
model: 16HKY__G test value: 542.089439902632\
model: 11HKY__G test value: 913.838186742487\
model: 10GTR__N test value: 1341.15917604558\
model: 12HKY__G test value: 1915.9724144655\
model: 14HKY__G test value: 1948.66018133137\
model: 13HKY__G test value: 1985.93211472103\
model: 15HKY__G test value: 2137.84559966153\
model: 16GTR__N test value: 6032.01897648639\
model: 11GTR__N test value: 10257.8384430192\
model: 13GTR__N test value: 22290.6543232418\
model: 12GTR__N test value: 22342.3701793176\
model: 14GTR__N test value: 24024.3092445402\
model: 15GTR__N test value: 24934.9691206518

Best model is 10GTR__G, test value: 7.00545534776865 
__________

If the empirical model (EM) is represented by replicates, the test_stage2.pl outputs a self-explanatory\
table which has three parts. An example of the table is provided below:

PART 1. Comparison of all simulation models (SMs) to each empirical model (EM)

MS1 values indicate percentage of times when the SM assuming EM tree topology + EM substitution model\
showed better fit to EM-based replicates in comparison to any other SM.

Shown are the mean test values estimated in comparisons of each SM to each EM-based replicate.

EM:10GTR__G SM:10GTR__G test value:7.9536552972                   MS1: N/A\
EM:10GTR__G SM:16GTR__G test value:28.9420718756                  MS1: 100.0\
EM:10GTR__G SM:11GTR__G test value:49.0069959436                  MS1: 100.0\
EM:10GTR__G SM:12GTR__G test value:95.9038204760                  MS1: 100.0\
EM:10GTR__G SM:13GTR__G test value:107.4520615314                 MS1: 100.0\
EM:10GTR__G SM:14GTR__G test value:110.7371941324                 MS1: 100.0\
EM:10GTR__G SM:15GTR__G test value:112.5031444522                 MS1: 100.0\
EM:10GTR__G SM:10HKY__G test value:136.4757790546                 MS1: 100.0\
EM:10GTR__G SM:16HKY__G test value:562.2581473864                 MS1: 100.0\
EM:10GTR__G SM:11HKY__G test value:919.8166452323                 MS1: 100.0\
EM:10GTR__G SM:10GTR__N test value:1331.8883894357                MS1: 100.0\
EM:10GTR__G SM:12HKY__G test value:1814.5569774804                MS1: 100.0\
EM:10GTR__G SM:13HKY__G test value:2037.6335797516                MS1: 100.0\
EM:10GTR__G SM:14HKY__G test value:2054.4223711852                MS1: 100.0\
EM:10GTR__G SM:15HKY__G test value:2140.2748792433                MS1: 100.0\
EM:10GTR__G SM:16GTR__N test value:6167.9382428951                MS1: 100.0\
EM:10GTR__G SM:11GTR__N test value:10141.4400140610               MS1: 100.0\
EM:10GTR__G SM:12GTR__N test value:20835.4701148395               MS1: 100.0\
EM:10GTR__G SM:13GTR__N test value:22429.3271299974               MS1: 100.0\
EM:10GTR__G SM:15GTR__N test value:24421.7637117566               MS1: 100.0\
EM:10GTR__G SM:14GTR__N test value:24514.6916720335               MS1: 100.0

EM:11GTR__G SM:11GTR__G test value:10.9878731993                  MS1: N/A\
EM:11GTR__G SM:16GTR__G test value:27.6538934487                  MS1: 100.0\
EM:11GTR__G SM:12GTR__G test value:47.0590774977                  MS1: 100.0\
EM:11GTR__G SM:10GTR__G test value:59.7485197459                  MS1: 100.0\
EM:11GTR__G SM:15GTR__G test value:61.2800660671                  MS1: 100.0\
EM:11GTR__G SM:13GTR__G test value:69.1213271538                  MS1: 100.0\
EM:11GTR__G SM:14GTR__G test value:84.5152367685                  MS1: 100.0\
EM:11GTR__G SM:11HKY__G test value:165.9773344128                 MS1: 100.0\
EM:11GTR__G SM:16HKY__G test value:408.7807462001                 MS1: 100.0\
EM:11GTR__G SM:12HKY__G test value:729.4230790075                 MS1: 100.0\
EM:11GTR__G SM:10HKY__G test value:899.1315594127                 MS1: 100.0\
EM:11GTR__G SM:15HKY__G test value:935.1123931735                 MS1: 100.0\
EM:11GTR__G SM:13HKY__G test value:1062.7946341026                MS1: 100.0\
EM:11GTR__G SM:14HKY__G test value:1260.7747463272                MS1: 100.0\
EM:11GTR__G SM:11GTR__N test value:1656.8597144909                MS1: 100.0\
EM:11GTR__G SM:16GTR__N test value:3860.5662158918                MS1: 100.0\
EM:11GTR__G SM:12GTR__N test value:9171.2265222537                MS1: 100.0\
EM:11GTR__G SM:10GTR__N test value:9750.9502802821                MS1: 100.0\
EM:11GTR__G SM:15GTR__N test value:11285.8252479354               MS1: 100.0\
EM:11GTR__G SM:13GTR__N test value:11621.6905898025               MS1: 100.0\
EM:11GTR__G SM:14GTR__N test value:14377.3317896123               MS1: 100.0

PART 2. Comparison of fit of simulation models sharing common substitution model scheme\
but assuming different model tree topologies

MS2 values indicate percentage of times when the SM assuming EM tree topology and a certain\
substitution model scheme, (termed 'preferred model') showed better fit to EM-based replicates\
in comparison to any other SM assuming the same substitution model scheme as the preferred model.

Shown are the mean test values estimated in comparisons of each SM to each EM-based replicate.

EM:10GTR__G SM:10GTR__G test value:7.9536552972                   MS2: N/A\
EM:10GTR__G SM:16GTR__G test value:28.9420718756                  MS2: 100.0\
EM:10GTR__G SM:11GTR__G test value:49.0069959436                  MS2: 100.0\
EM:10GTR__G SM:12GTR__G test value:95.9038204760                  MS2: 100.0\
EM:10GTR__G SM:13GTR__G test value:107.4520615314                 MS2: 100.0\
EM:10GTR__G SM:14GTR__G test value:110.7371941324                 MS2: 100.0\
EM:10GTR__G SM:15GTR__G test value:112.5031444522                 MS2: 100.0

EM:10GTR__G SM:10GTR__N test value:1331.8883894357                MS2: N/A\
EM:10GTR__G SM:16GTR__N test value:6167.9382428951                MS2: 100.0\
EM:10GTR__G SM:11GTR__N test value:10141.4400140610               MS2: 100.0\
EM:10GTR__G SM:12GTR__N test value:20835.4701148395               MS2: 100.0\
EM:10GTR__G SM:13GTR__N test value:22429.3271299974               MS2: 100.0\
EM:10GTR__G SM:15GTR__N test value:24421.7637117566               MS2: 100.0\
EM:10GTR__G SM:14GTR__N test value:24514.6916720335               MS2: 100.0

EM:10GTR__G SM:10HKY__G test value:136.4757790546                 MS2: N/A\
EM:10GTR__G SM:16HKY__G test value:562.2581473864                 MS2: 100.0\
EM:10GTR__G SM:11HKY__G test value:919.8166452323                 MS2: 100.0\
EM:10GTR__G SM:12HKY__G test value:1814.5569774804                MS2: 100.0\
EM:10GTR__G SM:13HKY__G test value:2037.6335797516                MS2: 100.0\
EM:10GTR__G SM:14HKY__G test value:2054.4223711852                MS2: 100.0\
EM:10GTR__G SM:15HKY__G test value:2140.2748792433                MS2: 100.0

EM:11GTR__G SM:11GTR__G test value:10.9878731993                  MS2: N/A\
EM:11GTR__G SM:16GTR__G test value:27.6538934487                  MS2: 100.0\
EM:11GTR__G SM:12GTR__G test value:47.0590774977                  MS2: 100.0\
EM:11GTR__G SM:10GTR__G test value:59.7485197459                  MS2: 100.0\
EM:11GTR__G SM:15GTR__G test value:61.2800660671                  MS2: 100.0\
EM:11GTR__G SM:13GTR__G test value:69.1213271538                  MS2: 100.0\
EM:11GTR__G SM:14GTR__G test value:84.5152367685                  MS2: 100.0

EM:11GTR__G SM:11GTR__N test value:1656.8597144909                MS2: N/A\
EM:11GTR__G SM:16GTR__N test value:3860.5662158918                MS2: 100.0\
EM:11GTR__G SM:12GTR__N test value:9171.2265222537                MS2: 100.0\
EM:11GTR__G SM:10GTR__N test value:9750.9502802821                MS2: 100.0\
EM:11GTR__G SM:15GTR__N test value:11285.8252479354               MS2: 100.0\
EM:11GTR__G SM:13GTR__N test value:11621.6905898025               MS2: 100.0\
EM:11GTR__G SM:14GTR__N test value:14377.3317896123               MS2: 100.0

EM:11GTR__G SM:11HKY__G test value:165.9773344128                 MS2: N/A\
EM:11GTR__G SM:16HKY__G test value:408.7807462001                 MS2: 100.0\
EM:11GTR__G SM:12HKY__G test value:729.4230790075                 MS2: 100.0\
EM:11GTR__G SM:10HKY__G test value:899.1315594127                 MS2: 100.0\
EM:11GTR__G SM:15HKY__G test value:935.1123931735                 MS2: 100.0\
EM:11GTR__G SM:13HKY__G test value:1062.7946341026                MS2: 100.0\
EM:11GTR__G SM:14HKY__G test value:1260.7747463272                MS2: 100.0

PART 3. Summary of the results

-Subtable 1 shows success (1) or failure (0) of identification of a preferred SM\
- sharing model tree topology with the EM - by the mean test values in comparison to all SMs\
which share a common substitution model scheme (shown above the subtable) with the preferred SM.

-Subtable 2 shows the mean test values for the preferred SMs assuming substitution\
model components shown above the subtable.

-Subtable 3 shows the worst model separation (MS2) values between the preferred SM\
and any other SM which assumes the same substitution model scheme (shown above the subtable)\
as the preferred SM, but different model tree topology.\
See PART2 of the report (above) for the full list of MS2 values.
Subtable 1\
|               | GTR__G  | HKY__G  | GTR__N  |
| ------------- |:------- |:------- |:------- |
| EM:10GTR__G   |    1    |    1    | 1       |
| EM:11GTR__G   |    1    |    1    | 1       |


Subtable 2\

|               | GTR__G  | HKY__G  | GTR__N  |
| ------------- |:------- |:------- |:------- |
| EM:10GTR__G   |    8    | 136     | 1332    |
| EM:11GTR__G   |    11   | 166     | 1657    |

Subtable 3\

|               | GTR__G  | HKY__G  | GTR__N  |
| ------------- |:------- |:------- |:------- |
| EM:10GTR__G   | 100	  | 100     | 100     |
| EM:11GTR__G   | 100	  | 100     | 100     |



The mean over values in above subtable is: 100.00

