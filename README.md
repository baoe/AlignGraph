##### Contents
[Overview] (#overview)  
[Copy right] (#copyright)  
[How to cite AlignGraph?] (#cite)  
[Short Manual] (#manual)  

<a name="overview"/>
### Overview
AlignGraph is a software that extends and joins contigs or scaffolds by reassembling them with help provided by a reference genome of a closely related organism.

<a name="copyright"/>
###Copy right
AlignGraph is under the [Artistic License 2.0](http://opensource.org/licenses/Artistic-2.0).

<a name="cite"/>
### How to cite AlignGraph?
If you use AlignGraph, please cite the following paper:  
AlignGraph: algorithm for secondary de novo genome assembly guided by closely related references (submitted).

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
   AlignGraph --read1 reads_1.fa --read2 reads_2.fa --contig contigs.fa --genome genome.fa --distanceLow distanceLow --distanceHigh distancehigh --extendedContig extendedContigs.fa --remainingContig remainingContigs.fa [--kMer k --insertVariation insertVariation --coverage coverage --noAlignment --part p --fastMap]
   ```

   Inputs:  
   --read1 is the the first pair of PE DNA reads in fasta format.  
   --read2 is the the second pair of PE DNA reads in fasta format.  
   --contig is the initial contigs/scaffolds in fasta format.  
   --genome is the reference genome in fasta format.  
   --distanceLow is the lower bound of alignment distance between the first and second pairs of PE DNA reads.  
   --distanceHigh is the upper bound of alignment distance between the first and second pairs of PE DNA reads.  
   Outputs:  
   --extendedContig is the extended contig/scaffold file in fasta format.  
   --remainingContig is the not extended initial contig/scaffold file in fasta format.  
   Options:  
   --kMer is the k-mer size (default: 5).  
   --insertVariation is the standard variation of insert length (default: 100).  
   --coverage is the minimum coverage to keep a path in de Bruijn graph (default: 20).  
   --noAlignment skips the initial time-consuming alignment step, if all the alignment files have been provided in tmp directory (default: none).  
   --part is the number of parts a chromosome is divided into when it is loaded to reduce memory requirement (default: 1).  
   --fastMap makes BLAT alignment faster but may lower the performance of AlignGraph. Useful for large genomes (default: none).  

5. Outputs
   * Extended contigs or scaffolds in FASTA format.
   * Remaining contigs or scaffolds not extended in FASTA format.
