### LATEST NEWS
NUCMER is supported to make the contig alignment to reference genome. Though less sensitive and accurate, it is much faster than BLAT and can effectively resolve the runtime bottleneck of BLAT.

### Overview
AlignGraph is a software that extends and joins contigs or scaffolds by reassembling them with help provided by a reference genome of a closely related organism.

### Copy right
AlignGraph is under the [Artistic License 2.0](http://opensource.org/licenses/Artistic-2.0).

### How to cite AlignGraph?
If you use AlignGraph, please cite the following paper:  
Bao E, Jiang T, Girke T (2014) AlignGraph: algorithm for secondary de novo genome assembly guided by closely related references. Bioinformatics: [epub](http://www.hubmed.org/display.cgi?uids=24932000).

### Short manual
1. System requirements

   AlignGraph is suitable for 32-bit or 64-bit machines with Linux operating systems. At least 4GB of system memory is recommended for assembling larger data sets.

2. Installation

   Aligners [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [BLAT](http://genome.ucsc.edu/FAQ/FAQblat.html)/[PBLAT](http://icebert.github.io/pblat/)/[NUCMER](http://mummer.sourceforge.net/)(BLAT's version should be v34 or below) are required to run AlignGraph.  
   * To use Bowtie2 and BLAT/PBLAT/NUCMER, put them to your $PATH: `export PATH=PATH2BOWTIE2:$PATH` and `export PATH=PATH2BLAT/PBLAT/NUCMER:$PATH`.
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
   --fastMap calls NUCMER to make fast but less sensitive and accurate contig alignment instead of BLAT (default: none).  
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

### Eval-AlignGraph
Eval-AlignGraph is the evaluation tool distributed with AlignGraph to generate statistics of the contigs or scaffolds. By default the contigs or scaffolds are aligned to the target genome by BLAT without the -fastMap option, but it can be enabled to avoid super long time waiting on some data with a little sensitivity loss by changing the value of FASTMAP macro from 0 to 1 in the source code.

### FAQs
1. How can I input multiple libraries with different insert lengths?

   Suppose you have one library x1.fa/y1.fa with insert length I1, and another library x2.fa/y2.fa with insert length I2, then you can simply combine x1.fa and x2.fa into x.fa, combine y1.fa and y2.fa into y.fa, and then input x.fa and y.fa. When you specify --distanceLow, let insert length be min{I1, I2}; for --distanceHigh, let insert length be max{I1, I2}.

2. Can I use mate pair libraries?

   You can but it is not recommended, since good alignment cannot be guaranteed with the very short reads and the very large insert length.

3. Which aligner do I need to put to $PATH: BLAT, PBLAT or NUCMER?

   BLAT and PBLAT are more sensitive and accurate than NUCMER, and PBLAT is the parallelized version of BLAT much faster, so PBLAT should be the first choice. However, the current PBLAT is not stable, so it is highly recommended to put both PBLAT and BLAT to your $PATH, so that AlignGraph can call BLAT instead when PBLAT fails somewhere. If PBLAT and BLAT are too slow to finish, you would have to switch to NUCMER and put it to your $PATH (remember to also specify the -fastMap option in command).

4. Why do I get the error "BLAT CALL FAILED!" even if I have put BLAT to my $PATH?

   The current version of BLAT (v35) is not compatible with AlignGraph, so you would have to use an earlier version (preferably v34) to avoid this error.

5. Why is there rare or no extension made by AlignGraph?

   How much extensions AlignGraph can make is mainly dependent on factors like how close the reference genome and the target genome are, and how well the pre-assembly worked. Therefore, it is possible there is rare or no extension, either because the reference genome is not so similar to the target genome, or because the upstream assemblies are already good enough for the current version of AlignGraph. We are currently working on improving AlignGraph's performance, so that more extensions can be made with a relatively different reference genome, but this may take some time.

6. How many threads are used for Bowtie2 and PBLAT?

   8 threads are used. Currently users cannot make changes to this, since this is a moderate choice for either single CPU machines (overhead for parallelization would not be too large) or multiple CPU machines.

### Erratum
   There is a small error in page i322. See [erratum](http://biocluster.ucr.edu/~ebao/erratum.pdf) for more details.

