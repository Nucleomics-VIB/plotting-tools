[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)
## plotting-tools

Create plots from numbers.

### **fpkm2heatmap.R** 

**[fpkm2heatmap.R](fpkm2heatmap/fpkm2heatmap.R)** creates a heatmap from **expXXXX-RNAseqCounts.xlsx** delivered with a RNASeq Nucleomics Core project (<a href="fpkm2heatmap/data/expXXXX-RNAseqCounts.xlsx">fake example with only 10 gene products</a>). An arbitrary example is shown next to ilustrate a typical result. Although care was taken, some minor edits may be required to adapt to changes in column order or eliminate 'toxic' data rows.

<img src="fpkm2heatmap/data/test-signature-heatmap.png?raw=true" alt="fpkm2heatmap example" style="width: 300px;"/>

**[fpkm2heatmap-cmd.R](fpkm2heatmap/fpkm2heatmap-cmd.R)** does the same but from terminal with command arguments

```{bash}
fpkm2heatmap-cmd.R -h
Usage: ./fpkm2heatmap-cmd.R [options] file

Options:
	-i INFILE, --infile=INFILE
		input file name [REQUIRED]

	-s SIGNATURE, --signature=SIGNATURE
		signature file name (if absent, all data will be plotted)

	-o OUTFILE, --outfile=OUTFILE
		base name for output [default: full_heatmap]

	-f OUTFORMAT, --outformat=OUTFORMAT
		file format for output 1:PNG, 2:PDF [default: 1]

	-h, --help
		Show this help message and exit
```

### **lendist.R** 

**[lendist.R](lendist/lendist.R)** creates two plots from one or more fasta assembly files. The plots report the cumulative sum of the contigs/chromosomes present in each fasta to illustrate assembly improvement.

<img src="lendist/data/assembly_sizes.png?raw=true" alt="lendist example" style="width: 300px;"/>
<img src="lendist/data/assembly_sizes-log.png?raw=true" alt="lendist example" style="width: 300px;"/>


*[[back-to-top](#top)]*  

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)

Licensed under the **GNU General Public License v3.0** — see [`LICENSE`](LICENSE) for the
full text, and [`NOTICE.md`](NOTICE.md) for attribution and the relicensing history.

Attribute the work to *VIB Nucleomics Core*, and licence derived work under the same terms.
Relicensed on 2026-08-27; copies obtained before that date remain available to their
holders under the previous CC BY-SA 3.0 terms.
