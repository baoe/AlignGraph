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
   AlignGraph --read1 reads_1.fa --read2 reads_2.fa --contig contigs.fa --genome genome.fa --distanceLow distanceLow --distanceHigh distancehigh --extendedContig extendedContigs.fa --remainingContig remainingContigs.fa [--kMer k --insertVariation insertVariation --coverage coverage --part p --fastMap --ratioCheck --iterativeMap --misassemblyRemoval --resume]
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
   --part is the number of parts a chromosome is divided into when it is loaded to reduce memory requirement (default: 1).  
   --fastMap makes BLAT alignment faster to avoid super long time waiting on some data but may lower a little sensitivity of AlignGraph (default: none).  
   --ratioCheck checks read alignment ratio to the reference beforehand and warns if the ratio is too low; may take a little more time (default: none).  
   --iterativeMap aligns reads to one chromosome and then another rather than directly to the genome, which increases sensitivity while loses precision (default: none).  
   --misassemblyRemoval detects and then breaks at or removes misassembed regions (default: none).  
   --resume resumes the previous unfinished running from several checkpoints (default: none).  

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

   Suppose you have one library x1.fa/y1.fa with insert length I1, and another library x2.fa/y2.fa with insert length I2, then you can simply combine x1.fa and x2.fa into x.fa, combine y1.fa and y2.fa into y.fa, and then input x.fa and y.fa. When you specify --insertLow, let insert length be min{I1, I2}; for --insertHigh, let insert length be max{I1, I2}.

2. Can I use mate pair libraries?

   You can but it is not recommended, since good alignment cannot be guaranteed with the very short reads and the very large insert length.

3. How many threads are used for Bowtie2?

   8 threads are used. Currently users cannot make changes to this, since this is a moderate choice for either single CPU machines (overhead for parallelization would not be too large) or multiple CPU machines. Another reason is, the bottleneck for the runtime is usually from BLAT, no matter how many threads there are for Bowtie2.

4. Why is there rare or no extension made by AlignGraph?

   How much extensions AlignGraph can make is mainly dependent on factors like how close the reference genome and the target genome are, and how well the pre-assembly worked. Therefore, it is possible there is rare or no extension, either because the reference genome is not so similar to the target genome, or because the upstream assemblies are already good enough for the current version of AlignGraph. We are currently working on improving AlignGraph's performance, so that more extensions can be made with a relatively different reference genome, but this may take some time. Besides these, you may want to check the bowtie_doc.txt and blat_doc.txt files to make sure Bowtie2 and BLAT were properly called by AlignGraph. 

5. Why could not my run with AlighGraph finish after a long time?

   This could be due to the runtime bottleneck of AlignGraph including the BLAT alignment from contigs to reference genome and the sequential processing after the alignment. For the first thing, a suggestion besides specifying the -fastMap option is to combine the reference sequences if the reference genome contains not the complete chromosomes but long contigs/scaffolds. Like FAQ 4 above, it may also take some more time to publish an accelerated version of AlignGraph. Here is a tip for you to estimate how much more time it may still need to finish: if there are outputs on screen like steps (1), (2), (3)..., it means AlignGraph is processing the reference chromosomes one by one and you will know the progress; if not, it means AlignGraph is still doing the alignment and you have to wait for more time.

<a name="error"/>
### Erratum
   There is a small error in page i322. See [erratum](http://biocluster.ucr.edu/~ebao/erratum.pdf) for more details.

