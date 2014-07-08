##### Contents
[Overview] (#overview)  
[Copy right] (#copyright)  
[How to cite AlignGraph?] (#cite)  
[Short Manual] (#manual)  
[Eval-AlignGraph] (#eval)
[FAQs] (#faq)
[Erratum] (#error)

<a name="overview"/>
### Overview
AlignGraph is a software that extends and joins contigs or scaffolds by reassembling them with help provided by a reference genome of a closely related organism.

<a name="copyright"/>
###Copy right
AlignGraph is under the [Artistic License 2.0](http://opensource.org/licenses/Artistic-2.0).

<a name="cite"/>
### How to cite AlignGraph?
If you use AlignGraph, please cite the following paper:  
Bao E, Jiang T, Girke T (2014) AlignGraph: algorithm for secondary de novo genome assembly guided by closely related references. Bioinformatics: [epub](http://www.hubmed.org/display.cgi?uids=24932000).

<a name="manual"/>
### Short manual
1. System requirements

   AlignGraph is suitable for 32-bit or 64-bit machines with Linux operating systems. At least 4GB of system memory is recommended for assembling larger data sets.

2. Installation

   Aligners [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [BLAT](http://genome.ucsc.edu/FAQ/FAQblat.html) are required to run AlignGraph.  
   * To use Bowtie2 and BLAT, put them to your $PATH: `export PATH=PATH2BOWTIE2:$PATH` and `export PATH=PATH2BLAT:$PATH`.
   * The downloaded AlignGraph.cpp file can be compiled with command `g++ -o AlignGraph AlignGraph.cpp -lpthread`.

3. Inputs
   * Paired-end DNA reads in FASTA format.
   * De novo contigs or scaffolds assembled by any de novo DNA-Seq assembler (Velvet, ABySS, ALLPATHS-LG, SOAPdenovo, etc.).
   * Reference genome from a closely related species.

4. Using AlignGraph

   ```
   AlignGraph --read1 reads_1.fa --read2 reads_2.fa --contig contigs.fa --genome genome.fa --distanceLow distanceLow --distanceHigh distancehigh --extendedContig extendedContigs.fa --remainingContig remainingContigs.fa [--kMer k --insertVariation insertVariation --coverage coverage --noAlignment --part p --fastMap --ratioCheck --iterativeMap --misassemblyRemoval]
   ```

   Inputs:  
   --read1 is the the first pair of PE DNA reads in fasta format.  
   --read2 is the the second pair of PE DNA reads in fasta format.  
   --contig is the initial contigs/scaffolds in fasta format.  
   --genome is the reference genome in fasta format.  
   --distanceLow is the lower bound of alignment distance between the first and second pairs of PE DNA reads (recommended: max{insert length - 1000, single read length}).  
   --distanceHigh is the upper bound of alignment distance between the first and second pairs of PE DNA reads (recommended: insert length + 1000).  
   Outputs:  
   --extendedContig is the extended contig/scaffold file in fasta format.  
   --remainingContig is the not extended initial contig/scaffold file in fasta format.  
   Options:  
   --kMer is the k-mer size (default: 5).  
   --insertVariation is the standard variation of insert length (default: 100).  
   --coverage is the minimum coverage to keep a path in de Bruijn graph (default: 20).  
   --noAlignment skips the initial time-consuming alignment step, if all the alignment files have been provided in tmp directory (default: none).  
   --part is the number of parts a chromosome is divided into when it is loaded to reduce memory requirement (default: 1).  
   --fastMap makes BLAT alignment faster to avoid super long time waiting on some data but may lower a little sensitivity of AlignGraph (default: none).  
   --ratioCheck checks read alignment ratio to the reference beforehand and warns if the ratio is too low; may take a little more time (default: none).  
   --iterativeMap aligns reads to one chromosome and then another rather than directly to the genome, which increases sensitivity while loses precision (default: none).  
   --misassemblyRemoval detects and then breaks at or removes misassembed regions (default: none).  

5. Outputs
   * Extended contigs or scaffolds in FASTA format. The format of the specification for each extended contig or scaffold (the string following the '>' of FASTA file) is: `AlignGraphX @ chromosomeID : contig/scaffoldID ; contig/scaffoldID ; contig/scaffoldID ... : partY`, where chromosomeID is the specification of the reference chromosome used to generate the extended contig or scaffold, X is a number starting from 0 to identify the extended contig or scaffold for each reference chromosome, and contig/scaffoldIDs are the specifications of the extendable contigs or scaffolds. If misassemblyRemoval is specified, partY shows the Y-th subcontig or subscaffold of the misassembled contig or scaffold split at misassemblies.
   * Remaining contigs or scaffolds not extended in FASTA format.

6. Example commandline

   Given PE reads files reads_1.fa and reads_2.fa with single read length 100 bp and insert length 500 bp, --distanceLow could be max{500 - 1000, 100} = 100 and --distanceHigh could be 500 + 1000 = 1500, so the simplest commandline with pre-assembled contigs file contigs.fa and reference genome genome.fa should be:

  ```
  AlignGraph --read1 reads_1.fa --read2 reads_2.fa --contig contigs.fa --genome genome.fa --distanceLow 100 --distanceHigh 1500 --extendedContig extendedContigs.fa --remainingContig remainingContigs.fa
  ```

<a name="eval"/>
### Eval-AlignGraph
Eval-AlignGraph is the evaluation tool distributed with AlignGraph to generate statistics of the contigs or scaffolds. By default the contigs or scaffolds are aligned to the target genome by BLAT without the -fastMap option, but it can be enabled to avoid super long time waiting on some data with a little sensitivity loss by changing the value of FASTMAP macro from 0 to 1 in the source code.

<a name="faq"/>
### FAQs
1. How can I input multiple libraries with different insert lengths?

   Suppose you have one library x1.fa/y1.fa with insert length I1, and another library x2.fa/y2.fa with insert length I2, then you can simply combine x1.fa and x2.fa into x.fa, combine x2.fa and y2.fa into y.fa, and then input x.fa and y.fa. When you specify --insertLow, let insert length be min{I1, I2}; for --insertHigh, let insert length be max{I1, I2}.

2. Can I use mate pair libraries?

   You can but it is not recommended, since good alignment cannot be guaranteed with the very short reads and the very large insert length.

3. How many threads are used for Bowtie2?

   8 threads are used. Currently users cannot make changes to this, since this is a moderate choice for either single CPU machines (overhead for parallelization would not be too large) or multiple CPU machines. Another reason is, the bottleneck for the runtime is usually from BLAT, no matter how many threads there are for Bowtie2.

<a name="error"/>
### Erratum
   There is a small error in page i322. See [erratum](http://biocluster.ucr.edu/~ebao/erratum.pdf) for more details.

