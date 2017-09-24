/*
 * Copyright 2013 Brian Tjaden
 *
 * This file is part of Rockhopper.
 *
 * Rockhopper is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * Rockhopper is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * (in the file gpl.txt) along with Rockhopper.  
 * If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.ArrayList;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import javax.swing.JTextArea;
import javax.swing.SwingWorker;

/**
 * The Rockhopper application reads in RNA-seq data, aligns the
 * sequencing reads to a genome, and then analyzes the data to
 * identify the expression of genes in each condition, differential
 * gene expression between conditions, transcript boundaries, and
 * operons. Rockhopper can also assemble transcripts de novo
 * without a reference genome.
 */
public class Rockhopper {

    /*****************************************
     **********   CLASS VARIABLES   **********
     *****************************************/

    public static final String version = "2.03";
    public static JTextArea output;  // Output to GUI or, if null, to System.out
    public static SwingWorker worker;  // Update progress bar if using GUI.
    private static int stagesProcessed;  // Number of files (stages) processed. Needed for GUI progress bar.
    public static ArrayList<Integer> genomeSizes;  // Sizes of each 1-indexed genome in nucleotides
    public static int genomeSize;  // Size of all genomes combined
    public static int numConditions;  // Number of experimental conditions
    private static ArrayList<String> genome_DIRs;  // PARAM: directories of geneome files
    public static String output_DIR = "Rockhopper_Results/";  // PARAM: write output files to this directory
    public static String browser_DIR = "genomeBrowserFiles/";  // Write browser files to this sub-directory
    private static ArrayList<String> conditionFiles;  // PARAM: list of seq-read files for each condition
    public static boolean computeExpression = true;  // PARAM: compute differential expression
    public static boolean computeOperons = true;     // PARAM: compute operons
    public static boolean computeTranscripts = true;  // PARAM: compute transcript boundaries
    private static PrintWriter summaryWriter = null;  // For outputting summary file
    private static String summaryFile = "summary.txt";
    private static String expressionFile = "transcripts.txt";
    private static String operonGenePairFile = "operonGenePairs.txt";
    private static String operonMergedFile = "operons.txt";
    public static boolean unstranded = false;  // Is RNA-seq data strand specific or ambiguous?
    public static double transcriptSensitivity = 0.5;
    public static boolean verbose = false;
    private static String[] labels;
    private static boolean time = false;
    private static boolean isDeNovo = false;  // Reference based (false) or de novo (true) assembly



    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private ArrayList<Genome> genomes = new ArrayList<Genome>();
    private ArrayList<Transcripts> transcripts = new ArrayList<Transcripts>();
    private ArrayList<Condition> conditions;  // RNA-seq data for each condition



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    /**
     * Constructs a new Rockhopper object.
     */
    public Rockhopper() {

	stagesProcessed = 0;
	if (isDeNovo) {  // De novo assembly
	    Assembler a = new Assembler();
	    return;
	}

	// Read in genome
	Rockhopper.genomeSizes = new ArrayList<Integer>();
	Rockhopper.genomeSize = 0;
	for (String genome_DIR : genome_DIRs) {
	    Genome genome = new Genome(genome_DIR);
	    this.genomes.add(genome);
	    Rockhopper.genomeSizes.add(genome.size());
	    Rockhopper.genomeSize += genome.size();
	}

	// Align sequencing reads to genome and compute basic expression information
	this.conditions = new ArrayList<Condition>();
	Peregrine.sequenceFiles = new String[genomes.size()];
	Peregrine.annotationsPlus = new String[genomes.size()][];
	Peregrine.annotationsMinus = new String[genomes.size()][];
	for (int z=0; z<genomes.size(); z++) {
	    String genome_DIR = genome_DIRs.get(z);
	    File dir = new File(genome_DIR);
	    String[] files = dir.list();
	    for (int i=0; i<files.length; i++) {
		if (files[i].endsWith(".fna")) {
		    Peregrine.sequenceFiles[z] = genome_DIR + files[i];
		    Peregrine.annotationsPlus[z] = genomes.get(z).getAnnotations('+');
		    Peregrine.annotationsMinus[z] = genomes.get(z).getAnnotations('-');
		}
	    }
	}
	Peregrine.numSequences = genomes.size();
	Peregrine.outputBrowserFile = false;
	output("\n");
	for (int i=0; i<conditionFiles.size(); i++) {  // For each condition
	    String[] readFiles = conditionFiles.get(i).split(",");
	    Condition C = new Condition();
	    for (int j=0; j<readFiles.length; j++) {  // For each replicate
		String[] parse_readFiles = readFiles[j].split("%");
		if (parse_readFiles.length == 1) {  // Using single-end reads (except maybe SAM/BAM)
		    String file1 = parse_readFiles[0];
		    output("Aligning sequencing reads from file:\t" + file1.split(separator())[file1.split(separator()).length-1] + "\n");
		    Peregrine.isPairedEnd = false;
		    Peregrine.pairedEndFile = null;
		    Peregrine.readsFile = file1;
		    Peregrine p = new Peregrine();
		    updateProgress();
		    C.addReplicate(new Replicate(p.getListOfCompressedFileNames(file1), file1, unstranded));
		} else {  // Using paired-end reads
		    String file1 = parse_readFiles[0];
		    String file2 = parse_readFiles[1];
		    output("Aligning sequencing reads from files:\n");
		    output("\t" + file1.split(separator())[file1.split(separator()).length-1] + "\n");
		    output("\t" + file2.split(separator())[file2.split(separator()).length-1] + "\n");
		    Peregrine.isPairedEnd = true;
		    Peregrine.readsFile = file1;
		    Peregrine.pairedEndFile = file2;
		    Peregrine p = new Peregrine();
		    updateProgress();
		    C.addReplicate(new Replicate(p.getListOfCompressedFileNames(file1), file1, unstranded));
		}
	    }
	    C.setMinDiffExpressionLevel();
	    conditions.add(C);
	}
	Rockhopper.numConditions = conditions.size();
	computeGeneExpression();

	// Output mapped reads to WIG files that can be loaded by genome browser
	double trackRange = 0.0;
	int totalReplicates = 0;
	for (int i=0; i<conditions.size(); i++) {
	    for (int j=0; j<conditions.get(i).numReplicates(); j++) {
		trackRange += conditions.get(i).getReplicate(j).getAvgReads();
		totalReplicates++;
	    }
	}
	trackRange = (2*trackRange) / totalReplicates;  // Range for display in genome browser
	for (int i=0; i<conditionFiles.size(); i++) {  // For each condition
	    String[] readFiles = conditionFiles.get(i).split(",");
	    for (int j=0; j<readFiles.length; j++) {  // For each replicate
		Replicate r = conditions.get(i).getReplicate(j);
		double multiplier = 100000.0 / (double)r.getUpperQuartile();
		for (int z=0; z<genomes.size(); z++) {
		    int[] coordinates_plus = new int[Rockhopper.genomeSizes.get(z)];
		    int[] coordinates_minus = new int[Rockhopper.genomeSizes.get(z)];
		    for (int k=0; k<Rockhopper.genomeSizes.get(z); k++) {
			coordinates_plus[k] = (int)(multiplier * r.getReads(z, k, '+'));
			coordinates_minus[k] = (int)(multiplier * r.getReads(z, k, '-'));
		    }
		    Peregrine.outputBrowserFile(readFiles[j].split("%")[0], Peregrine.sequenceFiles[z], output_DIR, browser_DIR, r.getCompressedFileName(z), coordinates_plus, coordinates_minus, (int)trackRange);
		}
	    }
	}
	Peregrine.releaseMemory();

	output("Analyzing transcripts..." + "\n");

	// Compute transcript boundaries
	if (computeTranscripts) {
	    for (int z=0; z<genomes.size(); z++) {
		Genome genome = genomes.get(z);
		if (genome.numGenes() > 0) {
		    Transcripts t = new Transcripts(z, genome, conditions, unstranded);
		    t.identifyUTRs();
		    t.identifyRNAs();
		    transcripts.add(t);
		}
	    }
	}

	// Compute gene differential expression
	if (computeExpression) {
	    Gene.setLowessVariances(genomes, conditions);  // Compute lowess variance for each gene
	    for (Genome genome : genomes) {
		if (genome.numGenes() > 0) computeDifferentialExpression(genome);
	    }
	    Gene.correctPvalues(genomes, conditions);  // Compute q-values, Benjamini-Hochberg correction
	}

	for (int z=0; z<genomes.size(); z++) runRockhopper(z);
	output("\nSummary of results written to file:\t" + output_DIR + summaryFile + "\n");
	for (int z=0; z<genomes.size(); z++) {
	    output("Expression of transcripts written to file:\t" + output_DIR + genomes.get(z).getID() + "_" + expressionFile + "\n");
	    if (computeOperons) {
		//output("Pairs of operon genes written to file:\t" + output_DIR + genome.getID() + "_" + operonGenePairFile + "\n");
		output("Operons written to file:\t\t" + output_DIR + genomes.get(z).getID() + "_" + operonMergedFile + "\n");
	    }
	}
	output("\nFINISHED.\n\n");
	summaryWriter.close();
	updateProgress();
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    public void runRockhopper(int z) {
	Genome genome = genomes.get(z);

	output("\n");
	output("************************************************************\n");
	output(genome.getName() + "\n");
	output("************************************************************\n");
	output("\n");

	// Compute transcript boundaries
	if (computeTranscripts && (genome.numGenes() > 0)) {
	    Transcripts t = transcripts.get(z);
	    output("Computing transcript boundaries...\n");
	    output("\tNumber of 5'UTRs:\t\t\t" + t.getNum5UTRs() + "\n");
	    output("\tNumber of 3'UTRs:\t\t\t" + t.getNum3UTRs() + "\n");
	    output("\tNumber of predicted RNAs:\t\t\t" + (t.getNumSenseRNAs() + t.getNumAntisenseRNAs()) + "\n");
	    if (!unstranded) {
		output("\tNumber of predicted RNAs (not antisense):\t\t" + t.getNumSenseRNAs() + "\n");
		output("\tNumber of predicted RNAs (antisense):\t\t" + t.getNumAntisenseRNAs() + "\n");
	    }
	    output("\n");
	    outputUTRsForBrowser(output_DIR + browser_DIR + genome.getID() + "_UTRs.wig", genome.getFormalGenomeName(), genome.getGenes(), genome.size());
	    outputRNAsForBrowser(output_DIR + browser_DIR + genome.getID() + "_ncRNAs.wig", genome.getFormalGenomeName(), genome.getGenes(), genome.size());
	}

	// Compute gene differential expression
	if (computeExpression && (genome.numGenes() > 0)) {
	    output("Computing differential gene expression...\n");
	    output("\tNumber of differentially expressed protein coding genes:\t" + getNumOfDifferentiallyExpressedGenes(genome, 0.01) + "\n\n");
	    outputDifferentiallyExpressedGenesForBrowser(output_DIR + browser_DIR + genome.getID() + "_diffExpressedGenes.wig", genome.getFormalGenomeName(), genome.getGenes(), genome.size());
	}

	// Output transcript file
	try {
	    PrintWriter writer = new PrintWriter(new File(output_DIR + genome.getID() + "_" + expressionFile));
	    writer.println(genome.genesToString(conditions, labels));
            writer.close();
	} catch (FileNotFoundException e) {
	    output("\nError - could not open file " + genome.getID() + "_" + expressionFile + "\n\n");
	}

	// Compute operons
	if (computeOperons && (genome.numGenes() > 0)) {
	    output("Computing likely operons...\n");
	    Operons ops = new Operons(genome.getCodingGenes(), genome.getGenes());
	    output("\tNumber of gene-pairs predicted to be part of the same operon:\t" + ops.getNumOperonGenePairs(genome.getCodingGenes()) + "\n");
	    //ops.outputGenePairOperons(output_DIR + genome.getID() + "_" + operonGenePairFile, genome.getCodingGenes());
	    int numMergedOperons = ops.outputMergedOperons(output_DIR + genome.getID() + "_" + operonMergedFile, output_DIR + browser_DIR + genome.getID() + "_operons.wig", genome.getFormalGenomeName(), genome.getCodingGenes(), genome.size());
	    output("\tNumber of predicted multi-gene operons:\t\t" + numMergedOperons + "\n\n");
	}
    }

    /**
     * Returns a list of Genome objects.
     */
    public ArrayList<Genome> getGenomes() {
	return genomes;
    }



    /**************************************************
     **********   PRIVATE INSTANCE METHODS   **********
     **************************************************/

    /**
     * For each gene, compute its raw read counts, normalized read
     * counts, mean, variance, and RPKM.
     */
    private void computeGeneExpression() {
	int totalReplicates = 0;
	long sumUpperQuartile = 0;
	for (int j=0; j<conditions.size(); j++) {
	    for (int k=0; k<conditions.get(j).numReplicates(); k++) {
		Replicate r = conditions.get(j).getReplicate(k);
		r.setMinExpression(transcriptSensitivity);
		ArrayList<Long> currentGeneReads = new ArrayList<Long>();  // For normalization
		int zeroCountGenes = 0;  // For normalization, keep track of genes with zero reads
		int totalGenes = 0;
		for (int z=0; z<genomes.size(); z++) {
		    Genome genome = genomes.get(z);
		    for (int i=0; i<genome.numGenes(); i++) {
			Gene g = genome.getGene(i);

			long readsForGene = r.getReadsInRange(z, g.getMinCoordinate(), g.getMaxCoordinate(), g.getStrand());
			if (unstranded)  // Strand ambiguous. Count reads on both strands.
			    readsForGene = r.getReadsInRange(z, g.getMinCoordinate(), g.getMaxCoordinate(), '+') + r.getReadsInRange(z, g.getMinCoordinate(), g.getMaxCoordinate(), '-');
			g.setRawCount(j, k, readsForGene);  // Set raw counts for gene
			if (r.getAvgLengthReads() == 0) g.setRawCount_reads(j, k, 0);
			else g.setRawCount_reads(j, k, readsForGene/r.getAvgLengthReads());
			currentGeneReads.add(readsForGene);  // Populate list for normalization
			if (readsForGene == 0) zeroCountGenes++;  // Genes with zero reads for normalization
		    }
		    totalGenes += genome.numGenes();
		}
		long upperQuartile = Misc.select_Long(currentGeneReads, (int)(zeroCountGenes+0.75*(totalGenes-zeroCountGenes)));
		r.setUpperQuartile(upperQuartile);
		sumUpperQuartile += upperQuartile;
		totalReplicates++;
	    }
	}
	Condition.setAvgUpperQuartile(sumUpperQuartile / totalReplicates);

	// Set normalized read counts for each gene
	for (int j=0; j<conditions.size(); j++) {
	    for (int k=0; k<conditions.get(j).numReplicates(); k++) {
		for (int z=0; z<genomes.size(); z++) {
		    Genome genome = genomes.get(z);
		    for (int i=0; i<genome.numGenes(); i++) {  // Set normalized counts
			genome.getGene(i).setNormalizedCount(j, k, 100000.0, conditions.get(j).getReplicate(k).getUpperQuartile());
		    }
		}
	    }
	}

	// Compute mean and RPKM for each gene in each condition
	for (int z=0; z<genomes.size(); z++) {
	    Genome genome = genomes.get(z);
	    for (int i=0; i<genome.numGenes(); i++) genome.getGene(i).computeExpression(conditions);
	}

	identifySimilarConditions();  // We require means to identify similar conditions

	// Compute variance for each gene in each condition
	for (int z=0; z<genomes.size(); z++) {
	    Genome genome = genomes.get(z);
	    for (int i=0; i<genome.numGenes(); i++) genome.getGene(i).computeVariance(conditions);
	}
    }

    private void computeDifferentialExpression(Genome genome) {
	// Compute differential expression for each gene in pairs of conditions
	for (int i=0; i<genome.numGenes(); i++) genome.getGene(i).computeDifferentialExpression();
    }

    /**
     * Returns the number of significantly differentially expressed
     * protein-coding genes.
     * The significance threshold is specified by the parameter.
     */
    private int getNumOfDifferentiallyExpressedGenes(Genome genome, double significance) {
	int numDiffExpressedGenes = 0;
	for (int i=0; i<genome.numGenes(); i++) {
	    if (genome.getGene(i).isDifferntiallyExpressedORF(significance))
		numDiffExpressedGenes++;
	}
	return numDiffExpressedGenes;
    }

    /**
     * For each condition, determines the other condition that is most similar,
     * i.e., its "partner". If there are no replicate experiments then the
     * partner of each condition is used as a surrogate replicate.
     * Similarity between two conditions is measured by the Pearson correlation 
     * coefficient of normalized gene expression.
     * This method sets the "partner" of each condition.
     */
    private void identifySimilarConditions() {
	if (conditions.size() == 1) {  // Special case. Only 1 condition. It partners itself.
	    conditions.get(0).setPartner(0);
	    return;
	}

	double[][] correlationMatrix = new double[conditions.size()][conditions.size()];
	for (int i=0; i<conditions.size(); i++) {
	    for (int j=0; j<conditions.size(); j++) {
		ArrayList<Long> e1 = new ArrayList<Long>();  // Gene expression values 1
		ArrayList<Long> e2 = new ArrayList<Long>();  // Gene expression values 2
		for (Genome genome : genomes) {
		    for (int k=0; k<genome.numGenes(); k++) {
			int geneLength = genome.getGene(k).getMaxCoordinate() - genome.getGene(k).getMinCoordinate() + 1;
			e1.add(genome.getGene(k).getMean(i)/geneLength);
			e2.add(genome.getGene(k).getMean(j)/geneLength);
		    }
		}
		correlationMatrix[i][j] = Misc.correlation(e1, e2);
	    }
	}

	// Determine partner for each condition
	int partner = -1;
	for (int i=0; i<conditions.size(); i++) {
	    if (i == 0) partner = 1;
	    else partner = 0;
	    for (int j=0; j<conditions.size(); j++) {
		if ((i != j) && (correlationMatrix[i][j] > correlationMatrix[i][partner]))
		    partner = j;
	    }
	    conditions.get(i).setPartner(partner);
	}
    }



    /***********************************************
     **********   PRIVATE CLASS METHODS   **********
     ***********************************************/

    /**
     * Create WIG file of predicted UTRs that
     * can be loaded into a genome browser.
     */
    private static void outputUTRsForBrowser(String outputFile, String genomeName, ArrayList<Gene> genes, int size) {
	// Determine coordinates of UTRs
	int[] geneCoordinates = new int[size];
	for (int i=0; i<genes.size(); i++) {
	    Gene g = genes.get(i);

	    // UTRs
	    if (g.isORF() && (g.getStartT() > 0) && (g.getStrand() == '+'))  // 5'UTR on plus strand
		for (int j=g.getStartT(); j<g.getStart(); j++) geneCoordinates[j] = 1;
	    if (g.isORF() && (g.getStartT() > 0) && (g.getStrand() == '-'))  // 5'UTR on minus strand
		for (int j=g.getStart()+1; j<=g.getStartT(); j++) geneCoordinates[j] = -1;
	    if (g.isORF() && (g.getStopT() > 0) && (g.getStrand() == '+'))  // 3'UTR on plus strand
		for (int j=g.getStop()+1; j<=g.getStopT(); j++) geneCoordinates[j] = 1;
	    if (g.isORF() && (g.getStopT() > 0) && (g.getStrand() == '-'))  // 3'UTR on minus strand
		for (int j=g.getStopT(); j<g.getStop(); j++) geneCoordinates[j] = -1;
	}

	// Output differentially expressed genes to genome browser file
	try {
	    PrintWriter writer = new PrintWriter(new File(outputFile));
	    writer.println("track name=" + "\"" + "UTRs" + "\"" + " color=255,0,255 altColor=255,0,255 graphType=bar viewLimits=-1:1");
	    writer.println("fixedStep chrom=" + genomeName + " start=1 step=1");
	    for (int j=1; j<geneCoordinates.length; j++) writer.println(geneCoordinates[j]);
	    writer.close();
	} catch (FileNotFoundException e) {
	    output("\nError - could not open file " + outputFile + "\n\n");
	}
    }

    /**
     * Create WIG file of predicted RNAs that
     * can be loaded into a genome browser.
     */
    private static void outputRNAsForBrowser(String outputFile, String genomeName, ArrayList<Gene> genes, int size) {
	// Determine coordinates of ncRNAs
	int[] geneCoordinates = new int[size];
	for (int i=0; i<genes.size(); i++) {
	    Gene g = genes.get(i);

	    // ncRNAs
	    if (!g.isORF() && (g.getName().equals("predicted RNA"))) {  // novel RNA
		int value = 1;
		if (g.getStrand() == '-') value = -1;
		for (int j=g.getFirst(); j<=g.getLast(); j++) geneCoordinates[j] = value;
	    }
	}

	// Output differentially expressed genes to genome browser file
	try {
	    PrintWriter writer = new PrintWriter(new File(outputFile));
	    writer.println("track name=" + "\"" + "Novel RNAs" + "\"" + " color=0,255,0 altColor=0,255,0 graphType=bar viewLimits=-1:1");
	    writer.println("fixedStep chrom=" + genomeName + " start=1 step=1");
	    for (int j=1; j<geneCoordinates.length; j++) writer.println(geneCoordinates[j]);
	    writer.close();
	} catch (FileNotFoundException e) {
	    output("\nError - could not open file " + outputFile + "\n\n");
	}
    }

    /**
     * Create WIG file of differentially expressed genes that
     * can be loaded into a genome browser.
     */
    private static void outputDifferentiallyExpressedGenesForBrowser(String outputFile, String genomeName, ArrayList<Gene> genes, int size) {
	// Determine coordinates of differentially expressed genes
	int[] geneCoordinates = new int[size];
	for (int i=0; i<genes.size(); i++) {
	    Gene g = genes.get(i);
	    double qValue = g.getMinQvalue();
	    int value = 0;
	    if (qValue == 0.0) value = 300;  // Special case. We cannot take log of zero.
	    else value = (int)(-Math.log10(qValue));
	    if (g.getStrand() == '-') value = -value;
	    for (int j=g.getFirst(); j<=g.getLast(); j++) geneCoordinates[j] = value;
	}

	// Output differentially expressed genes to genome browser file
	try {
	    PrintWriter writer = new PrintWriter(new File(outputFile));
	    writer.println("track name=" + "\"" + "Differentially expressed genes" + "\"" + " color=0,255,0 altColor=0,255,0 graphType=bar viewLimits=-10:10");
	    writer.println("fixedStep chrom=" + genomeName + " start=1 step=1");
	    for (int j=1; j<geneCoordinates.length; j++) writer.println(geneCoordinates[j]);
	    writer.close();
	} catch (FileNotFoundException e) {
	    output("\nError - could not open file " + outputFile + "\n\n");
	}
    }

    /**
     * Returns the system-dependent file separator that
     * can be used as a RegEx.
     */
    private static String separator() {
	if (File.separatorChar == '\\') return "\\\\";
	else return File.separator;
    }



    /**********************************************
     **********   PUBLIC CLASS METHODS   **********
     **********************************************/

    // Used for updating the GUI progress bar.
    public static void updateProgress() {
	if (worker != null) {
	    double totalStages = 0.0;
	    for (String files : conditionFiles) totalStages += files.split(",").length;
	    if (isDeNovo) {
		for (String files : conditionFiles) totalStages += files.split(",").length;
		for (String files : conditionFiles) totalStages += files.split(",").length;
	    }
	    totalStages += 1.0;  // For analyzing transcripts and computing expression
	    stagesProcessed++;
	    worker.firePropertyChange("progress", worker.getProgress(), (int)(100.0*stagesProcessed/totalStages));
	}
    }

    public static void commandLineArguments(String[] args) {
	if (args.length < 1) {
	    output("\n*************************************************\n");
	    output("**********   Rockhopper version " + Rockhopper.version + "   **********\n");
	    output("*************************************************\n");
	    output("\nThe Rockhopper application has the following required command line arguments.\n");

	    output("\nREQUIRED ARGUMENTS\n\n");
	    output("\t" + "exp1A.fastq,exp1B.fastq,exp1C.fastq  exp2A.fastq,exp2B.fastq" + "\t" + "a comma separated list of sequencing files (in FASTQ, QSEQ, FASTA, SAM, or BAM format) for replicate experiments, one list per experimental condition (mate-pair files should be delimited by '%')\n");

	    output("\nREFERENCE BASED ASSEMBLY VS. DE NONO ASSEMBLY:\n");
	    output("IF THE -g OPTION IS USED THEN ROCKHOPPER ALIGNS READS TO ONE OR MORE REFERENCE GENOMES,\n");
	    output("OTHERWISE, ROCKHOPPER PERFORMS DE NOVO TRANSCRIPT ASSEMBLY.\n\n");
	    output("\t" + "-g <DIR1,DIR2>" + "\t" + "a comma separated list of directories, each containing a genome file (*.fna), gene file (*.ptt), and rna file (*.rnt)\n");

	    output("\nOPTIONAL ARGUMENTS FOR EITHER REFERENCE BASED ASSEMBLY OR DE NOVO ASSEMBLY\n\n");
            output("\t" + "-c <boolean>" + "\t" + "reverse complement single-end reads (default is false)\n");
            output("\t" + "-ff/fr/rf/rr" + "\t" + "orientation of two mate reads for paired-end read, f=forward and r=reverse_complement (default is fr)\n");
            output("\t" + "-d <integer>" + "\t" + "maximum number of bases between mate pairs for paired-end reads (default is 500)\n");
            output("\t" + "-a <boolean>" + "\t" + "identify 1 alignment (true) or identify all optimal alignments (false), (default is true)\n");
            output("\t" + "-p <integer>" + "\t" + "number of processors (default is self-identification of processors)\n");
	    output("\t" + "-e <boolean>" + "\t" + "compute differential expression for transcripts in pairs of experimental conditions (default is true)\n");
	    output("\t" + "-s <boolean>" + "\t" + "RNA-seq experiments are strand specific (true) or strand ambiguous (false), (default is true)\n");
	    output("\t" + "-L <comma separated list>" + "\t" + "labels for each condition\n");
	    output("\t" + "-o <DIR>" + "\t" + "directory where output files are written (default is Rockhopper_Results/)\n");
	    output("\t" + "-v <boolean>" + "\t" + "verbose output including raw/normalized counts aligning to each gene (default is false)\n");
            output("\t" + "-SAM       " + "\t" + "output a SAM format file\n");
            output("\t" + "-TIME      " + "\t" + "output time taken to execute program\n");

	    output("\nOPTIONAL ARGUMENTS FOR REFERENCE BASED ASSEMBLY ONLY\n\n");
            output("\t" + "-m <number>" + "\t" + "allowed mismatches as percent of read length (default is 0.15)\n");
	    output("\t" + "-l <number>" + "\t" + "minimum seed as percent of read length (default is 0.33)\n");
	    output("\t" + "-y <boolean>" + "\t" + "compute operons (default is true)\n");
	    output("\t" + "-t <boolean>" + "\t" + "identify transcript boundaries including UTRs and ncRNAs (default is true)\n");
	    output("\t" + "-z <number>" + "\t" + "minimum expression of UTRs and ncRNAs, a number in range [0.0, 1.0] (default is 0.5)\n");

	    output("\nOPTIONAL ARGUMENTS FOR DE NOVO ASSEMBLY ONLY\n\n");
	    output("\t" + "-k <integer>" + "\t" + "size of k-mer, range of values is 15 to 31 (default is 25)\n");
            output("\t" + "-j <integer>" + "\t" + "minimum length required to use a sequencing read after trimming/processing (default is 35)\n");
            output("\t" + "-n <integer>" + "\t" + "size of k-mer hashtable is ~ 2^n (default is 25). HINT: should normally be 25 or, if more memory is available, 26. WARNING: if increased above 25 then more than 1.2M of memory must be allocated\n");
            output("\t" + "-b <integer>" + "\t" + "minimum number of full length reads required to map to a de novo assembled trancript (default is 20)\n");
	    output("\t" + "-u <integer>" + "\t" + "minimum length of de novo assembled transcripts (default is 2*k)\n");
	    output("\t" + "-w <integer>" + "\t" + "minimum count of k-mer to use it to seed a new de novo assembled transcript (default is 50)\n");
	    output("\t" + "-x <integer>" + "\t" + "minimum count of k-mer to use it to extend an existing de novo assembled transcript (default is 5)\n");

	    output("\nEXAMPLE EXECUTION: REFERENCE BASED ASSEMBLY WITH SINGLE-END READS\n");
	    output("\njava Rockhopper <options> -g genome_DIR1,genome_DIR2 aerobic_replicate1.fastq,aerobic_replicate2.fastq anaerobic_replicate1.fastq,anaerobic_replicate2.fastq\n");
	    output("\nEXAMPLE EXECUTION: REFERENCE BASED ASSEMBLY WITH PAIRED-END READS\n");
	    output("\njava Rockhopper <options> -g genome_DIR1,genome_DIR2 aerobic_replicate1_pairedend1.fastq%aerobic_replicate1_pairedend2.fastq,aerobic_replicate2_pairedend1.fastq%aerobic_replicate2_pairedend2.fastq anaerobic_replicate1_pairedend1.fastq%anaerobic_replicate1_pairedend2.fastq,anaerobic_replicate2_pairedend1.fastq%anaerobic_replicate2_pairedend2.fastq\n");
	    output("\nEXAMPLE EXECUTION: DE NOVO ASSEMBLY WITH SINGLE-END READS\n");
	    output("\njava Rockhopper <options> aerobic_replicate1.fastq,aerobic_replicate2.fastq anaerobic_replicate1.fastq,anaerobic_replicate2.fastq\n");
	    output("\nEXAMPLE EXECUTION: DE NOVO ASSEMBLY WITH PAIRED-END READS\n");
	    output("\njava Rockhopper <options> aerobic_replicate1_pairedend1.fastq%aerobic_replicate1_pairedend2.fastq,aerobic_replicate2_pairedend1.fastq%aerobic_replicate2_pairedend2.fastq anaerobic_replicate1_pairedend1.fastq%anaerobic_replicate1_pairedend2.fastq,anaerobic_replicate2_pairedend1.fastq%anaerobic_replicate2_pairedend2.fastq\n");

	    output("\n");
	    System.exit(0);
	}

	// Initially set the number of threads
	Peregrine.numThreads = Runtime.getRuntime().availableProcessors();
	if (Peregrine.numThreads > 4) Peregrine.numThreads *= 0.75;
	Assembler.numThreads = Runtime.getRuntime().availableProcessors();
	if (Assembler.numThreads > 4) Assembler.numThreads = (int)Math.min(Assembler.numThreads * 0.75, 8.0);

	int i = 0;
	conditionFiles = new ArrayList<String>();
	genome_DIRs = null;
	while (i < args.length) {
	    if (args[i].equals("-g")) {
		if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
		    output("Error - command line argument -g must be followed by the name of one or more directories.\n");
		    System.exit(0);
		}
		String[] parse_dirs = args[i+1].split(",");
		genome_DIRs = new ArrayList<String>();
		for (String s : parse_dirs) {
		    if (s.endsWith(File.separator)) genome_DIRs.add(s);
		    else genome_DIRs.add(s + File.separatorChar);
		}
		isDeNovo = false;
		i += 2;
            } else if (args[i].startsWith("-a")) {
                if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -a must be followed by a boolean.\n");
                    System.exit(0);
                }
                Peregrine.stopAfterOneHit = Boolean.valueOf(args[i+1]);
                Assembler.stopAfterOneHit = Boolean.valueOf(args[i+1]);
                i += 2;
            } else if (args[i].startsWith("-p")) {
                if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -p must be followed by an integer.\n");
                    System.exit(0);
                }
                Peregrine.numThreads = Integer.parseInt(args[i+1]);
                Assembler.numThreads = Integer.parseInt(args[i+1]);
                i += 2;
            } else if (args[i].startsWith("-e")) {
                if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -e must be followed by a boolean.\n");
                    System.exit(0);
                }
		computeExpression = Boolean.valueOf(args[i+1]);
		Assembler.computeExpression = Boolean.valueOf(args[i+1]);
                i += 2;
            } else if (args[i].startsWith("-s")) {
                if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -s must be followed by a boolean.\n");
                    System.exit(0);
                }
		unstranded = !Boolean.valueOf(args[i+1]);
		Assembler.unstranded = !Boolean.valueOf(args[i+1]);
                i += 2;
            } else if (args[i].startsWith("-d")) {
                if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -d must be followed by an integer.\n");
                    System.exit(0);
                }
                Peregrine.maxPairedEndLength = Integer.parseInt(args[i+1]);
                Assembler.maxPairedEndLength = Integer.parseInt(args[i+1]);
		SamOps.maxPairedEndLength = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].equals("-L")) {
		if (i == args.length-1) {
		    output("Error - command line argument -L must be followed by a comma separated list of names for the conditions.\n");
		    System.exit(0);
		}
		labels = args[i+1].split(",");
		Assembler.labels = args[i+1].split(",");
		i += 2;
	    } else if (args[i].equals("-o")) {
		if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
		    output("Error - command line argument -o must be followed the name of a directory.\n");
		    System.exit(0);
		}
		if (args[i+1].endsWith("/")) output_DIR = args[i+1]; 
		else output_DIR = args[i+1] + "/"; 
		Peregrine.outputDIR = output_DIR;
		Assembler.output_DIR = output_DIR;
		i += 2;
            } else if (args[i].startsWith("-v")) {
                if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -v must be followed by a boolean.\n");
                    System.exit(0);
                }
		verbose = Boolean.valueOf(args[i+1]);
		Assembler.verbose = Boolean.valueOf(args[i+1]);
                i += 2;
            } else if (args[i].startsWith("-m")) {
                if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -m must be followed by a decimal number.\n");
                    System.exit(0);
                }
                Peregrine.percentMismatches = Double.parseDouble(args[i+1]);
                i += 2;
            } else if (args[i].startsWith("-l")) {
                if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -l must be followed by a decimal number.\n");
                    System.exit(0);
                }
                Peregrine.percentSeedLength = Double.parseDouble(args[i+1]);
                i += 2;
            } else if (args[i].startsWith("-y")) {
                if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -y must be followed by a boolean.\n");
                    System.exit(0);
                }
		computeOperons = Boolean.valueOf(args[i+1]);
                i += 2;
            } else if (args[i].startsWith("-t")) {
                if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -t must be followed by a boolean.\n");
                    System.exit(0);
                }
		computeTranscripts = Boolean.valueOf(args[i+1]);
                i += 2;
	    } else if (args[i].equals("-z")) {
		if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
		    output("Error - command line argument -z must be followed by a number in the range [0.0,1.0].\n");
		    System.exit(0);
		}
		transcriptSensitivity = Double.parseDouble(args[i+1]);
		i += 2;
	    } else if (args[i].startsWith("-k")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -k must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.k = Integer.parseInt(args[i+1]);
		if (Assembler.minTranscriptLength == 0) Assembler.minTranscriptLength = 2*Assembler.k;
                i += 2;
	    } else if (args[i].startsWith("-j")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -j must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.minReadLength = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-n")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -n must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.CAPACITY_POWER = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-b")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -b must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.MIN_READS_MAPPING = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-u")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -u must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.minTranscriptLength = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-w")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -w must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.minSeedExpression = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-x")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -x must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.minExpression = Integer.parseInt(args[i+1]);
                i += 2;
            } else if (args[i].startsWith("-c")) {
                if ((i == args.length-1) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -c must be followed by a boolean.\n");
                    System.exit(0);
                }
                Peregrine.singleEndOrientationReverseComplement = Boolean.valueOf(args[i+1]);
                Assembler.singleEndOrientationReverseComplement = Boolean.valueOf(args[i+1]);
                i += 2;
	    } else if (args[i].equals("-ff")) {
		Peregrine.pairedEndOrientation = "ff";
		Assembler.pairedEndOrientation = "ff";
		i++;
	    } else if (args[i].equals("-fr")) {
		Peregrine.pairedEndOrientation = "fr";
		Assembler.pairedEndOrientation = "fr";
		i++;
	    } else if (args[i].equals("-rf")) {
		Peregrine.pairedEndOrientation = "rf";
		Assembler.pairedEndOrientation = "rf";
		i++;
	    } else if (args[i].equals("-rr")) {
		Peregrine.pairedEndOrientation = "rr";
		Assembler.pairedEndOrientation = "rr";
		i++;
	    } else if (args[i].equals("-TIME")) {
		time = true;
		Assembler.time = true;
		i++;
	    } else if (args[i].equals("-SAM")) {
		Peregrine.outputSAM = true;
		Assembler.outputSAM = true;
		i++;
	    } else {
		conditionFiles.add(args[i]);
		i++;
	    }
	}

	// Handle erroneous command line arguments
	if (genome_DIRs == null) {
	    isDeNovo = true;
	} else {
	    for (int j=0; j<genome_DIRs.size(); j++) {
		if (!(new File(genome_DIRs.get(j)).exists())) {
		    output("Error - directory " + genome_DIRs.get(j) + " does not exist.\n");
		    System.exit(0);
		}
	    }
	}
	if (conditionFiles.size() == 0) {
	    output("Error - sequencing reads files (in FASTQ, QSEQ, FASTA, SAM, or BAM format) are required as command line arguments.\n");
	    System.exit(0);
	}
	Assembler.conditionFiles = conditionFiles;
	Assembler.expressionFile = expressionFile;
	Assembler.output_DIR = output_DIR;
	Peregrine.outputDIR = output_DIR;
	for (int j=0; j<conditionFiles.size(); j++) {
	    String[] files = conditionFiles.get(j).split(",");
	    for (int k=0; k<files.length; k++) {
		String[] pairedEnd_files = files[k].split("%");
		if (!(new File(pairedEnd_files[0]).exists())) {
		    output("Error - file " + pairedEnd_files[0] + " does not exist.\n");
		    System.exit(0);
		}
		if ((pairedEnd_files.length > 1) && (!(new File(pairedEnd_files[1]).exists()))) {
		    output("Error - file " + pairedEnd_files[1] + " does not exist.\n");
		    System.exit(0);
		}
		if (pairedEnd_files.length > 1) Assembler.unstranded = false;  // Assume paired-end reads are strand specific
	    }
	}

	// If not set by command line arguments, set de novo assembled transcript length
	if (Assembler.minTranscriptLength == 0) Assembler.minTranscriptLength = 2*Assembler.k;

	// Output directory does not exist. Create it.
	if (!(new File(output_DIR).exists())) {
	    if (!(new File(output_DIR).mkdir())) {
		output("Error - could not create directory " + output_DIR + "." + "\n");
		System.exit(0);
	    }
	}
	Peregrine.browserDIR = Rockhopper.browser_DIR;
	if (!(new File(output_DIR + browser_DIR).exists())) {
	    if (!(new File(output_DIR + browser_DIR).mkdir())) {
		output("Error - could not create directory " + (output_DIR+browser_DIR) + "." + "\n");
		System.exit(0);
	    }
	}
	try {
	    summaryWriter = new PrintWriter(new File(output_DIR + summaryFile));
	    Peregrine.summaryWriter = summaryWriter;
	    Assembler.summaryWriter = summaryWriter;
	} catch (FileNotFoundException e) {
	    output("Error - could not create file " + (output_DIR+summaryFile) + "\n");
	    System.exit(0);
	}

    }

    public static void output(String s) {
	if (summaryWriter != null) summaryWriter.print(s);
	if (output == null) System.out.print(s);
	else {
	    output.append(s);
	    output.setCaretPosition(output.getDocument().getLength());
	}
    }



    /*************************************
     **********   MAIN METHOD   **********
     *************************************/

    /**
     * The Main method, when invoked with the appropriate command line
     * arguments, executes the Rockhopper application for the purpose
     * of analyzing RNA-seq data.
     */
    public static void main(String[] args) {
	commandLineArguments(args);
	long runningTime = System.currentTimeMillis();
	Rockhopper rh = new Rockhopper();
	runningTime = System.currentTimeMillis() - runningTime;
	if (time) {  // Output time taken to execute program
            System.out.println("Execution time:\t" + (runningTime/60000) + " minutes " + ((runningTime%60000)/1000) + " seconds" + "\n");
	}
    }

}

