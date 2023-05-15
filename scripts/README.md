# How to run
## Estimate Top N Genefamilies/Pathways abundance
The top N genefamilies/pathways abundance across the samples can be estimated by using the `extract_topNabundance.py`.
```
python3 extract_topNabundance.py --help

usage: extract_topNabundance.py [-h] -i INPUT -o OUTPUT [-st STATS] [-gs GSTRING] [-ps PSTRING] [-qc QC] [-m METADATA] [-t TOPN] [--version]

optional arguments:
  -h,  --help             Show this help message and exit
  
  -i,  --input      Gene Families/Path(way) Abundance summary file across multiple samples
  -o,  --output    Prefix for output file
  -st, --stats     A tsv file containing sample names and total number of reads per sample.
  -gs, --gstring   String suffix to trim from sample names in genefamilies abundance summary file (default: _Abundance-RPKs)
  -ps, --pstring   String suffix to trim from sample names in PathAbundance summary file (default: _Abundance)
  -qc, --qc        Save the qc report
  -m, --metadata   Metadata about the samples in csv format with colum name: 'SampleId' containing the sample names
  -t, --topN       Top N genefamilies (default N = 20)
  --version        Show program's version number and exit
```
### Estimate top N genefamilies abundance
To get the top N (for e.g. N = 5) genefamilies from the Genefamilies summary file obtained from Humann2/3 pipeline.
* By default the script trim **_Abundance-RPKs** suffix from samples names from genefamilies summary file. If there is something else then use `-gs` flag to change. 
```
python3 extract_topNabundance.py -i ./test/input/Sample_Genefamilies.tsv -o top5gf -st stats.csv -t 5

Outputs:
top5gf-rel.csv: Top5 Genefamilies relative abundance across the samples
top5gf-abs.csv: Top5 Genefamilies absolute abundance across the samples
qc-report.csv: QC file --> Samplenames, # Total_reads, # Unmapped reads, # Mapped reads
```
Note:
* The stats.csv file is the file containing total number of reads.
```
For e.g. Sample Name, Total no of reads
Sample1	1000
Sample2	950
Sample3	900
``` 
* The flag `-qc` and `-st` should only be used when using Genefamilies summary as input. But these are optional.

### Estimate top N pathways abundance
To get the top N (for e.g. N = 5) pathways from the Pathabundance summary file obtained from Humann2/3 pipeline.
* By default the script trim **_Abundance** suffix from samples names from PathAbundance summary file. If there is something else then use `-ps` flag to change.
```
python3 extract_topNabundance.py -i ./test/input/Sample_Pathabundance.tsv -o top5pa -t 5

Outputs:
top5pa-rel.csv: Top5 Pathways relative abundance across the samples
top5gf-abs.csv: Top5 Pathways absolute abundance across the samples
```

### Plot the top N genefamiles/Pathways relative abundance
To make a bubble plot showing the relative abundance of Top N genefamiles/Pathways.
```
Rscript bubble_plot.R -h

Usage: bubble_plot.R [options]

Options:
	-i, --input Genefamilies/Pathways/taxa Relative abundance file
	-d, --datatype  Input datatype:For GeneFamiles -> GF, Pathways -> PW, Taxa -> T, [default=GF]
	-m, --metadata  If INPUT file has metadata in it.[default= FALSE]
	-q, --qcfile    QC file with SampleName, Total_Reads, Unmapped, Mapped
	-o, --output    Prefix for output files [defualt=bubble_plot]
	-h, --help      Show this help message and exit
```
An example:
To plot the top 5 genefamilies relative abundance (i.e., top5gf-rel.csv) and the qc-report (i.e., qc-report.csv) generated above
* Note Change the datatype `-d` as per the input data.
* `-q` is optional.
```
Rscript bubble_plot.R -i top5-rel.csv -d GF -q qc-repot.csv -o top5gf_plot

Output:
top5gf_plot.svg: Top 5 Genefamilies abundance plot.
qc_plot.svg: Its default output file genereate showing the barplots for mapped vs unmapped reads per sample if "-q" is used.
```

### Adding metadata
`humann_barplot` script required summarized humann files (genefamilies/pathway abundance)
with the metadadata at the end of the those files. This is helpful if the grouping of samples as per the metadata is inteded as final plot.

Use the script to add the metadata of samples to summary files from `humann_join_tables`.

This script is not part of the main pipeline.

```
# The input shoul be provided using the command line.

usage: python3 addmetada.py

Humann summary file (.tsv) name: filename or /path_to_file/name
Metadata (in .csv format) filename: filename or /path_to_file/name
Output file name prefix: some_string
```
