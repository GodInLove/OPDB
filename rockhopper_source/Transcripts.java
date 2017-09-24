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
import java.util.HashMap;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;

/**
 * The Transcript class identifies UTRs of genes as well
 * ncRNA transcripts based on RNA-seq data.
 */
public class Transcripts {

    /*****************************************
     **********   CLASS VARIABLES   **********
     *****************************************/

    //private static int WINDOW = 10;  // Number of nucleotides in sliding window



    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private Genome genome;
    private int z;  // Index of genome in list of all genomes
    private ArrayList<Condition> conditions;
    private boolean unstranded;        // Is RNA-seq data strand specific or ambiguous?
    private int num5UTRs = 0;          // Number of 5'UTRs
    private int num3UTRs = 0;          // Number of 3'UTRs
    private int numSenseRNAs = 0;      // Number of sense ncRNAs
    private int numAntisenseRNAs = 0;  // Number of antisense ncRNAs
    //public long[] distribution = new long[201];  // Distribution of expressed gene reads



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    public Transcripts(int z, Genome genome, ArrayList<Condition> conditions, boolean unstranded) {
	this.z = z;
	this.genome = genome;
	this.conditions = conditions;
	this.unstranded = unstranded;
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    /**
     * Returns the number of 5'UTRs.
     */ 
    public int getNum5UTRs() {
	return this.num5UTRs;
    }

    /**
     * Returns the number of 3'UTRs.
     */ 
    public int getNum3UTRs() {
	return this.num3UTRs;
    }

    /**
     * Returns the number of sense ncRNAs.
     */ 
    public int getNumSenseRNAs() {
	return this.numSenseRNAs;
    }

    /**
     * Returns the number of antisense ncRNAs.
     */ 
    public int getNumAntisenseRNAs() {
	return this.numAntisenseRNAs;
    }

    /**
     * For each gene, identify its 5'UTR and 3'UTR based on
     * the expression data.
     */
    public void identifyUTRs() {
	for (int x=0; x<genome.numGenes()+1; x++) {  // For each IG region
	    Gene g1 = null;
	    Gene g2 = null;

	    // Identify downstream gene and IG stop coordinate
	    if (x == genome.numGenes()) g2 = new Gene((genome.size()-1) + ".." + (genome.size()-1) + "\t" + genome.getGene(x-1).getStrand() + "\t" + "0" + "\t" + "-" + "\t" + "???" + "\t" + "???" + "\t" + "-" + "\t" + "-" + "\t" + "???", "ORF");
	    else g2 = genome.getGene(x);
	    int stop = Math.min(g2.getStart(), g2.getStop()) - 1;
	    if (!g2.isORF()) stop = g2.getMinCoordinate() - 1;

	    // Identify upstream gene and IG start coordinate
	    int start = 1;
	    int upstreamGeneIndex = x-1;
	    while (upstreamGeneIndex >= 0) {
		if (unstranded) {  // Strand ambiguous
		    start = Math.max(genome.getGene(upstreamGeneIndex).getStart(), genome.getGene(upstreamGeneIndex).getStop()) + 1;
		    if (!genome.getGene(upstreamGeneIndex).isORF()) start = genome.getGene(upstreamGeneIndex).getMaxCoordinate() + 1;
		    break;
		} else {  // Strand specific
		    if (g2.getStrand() == genome.getGene(upstreamGeneIndex).getStrand()) {
			start = Math.max(genome.getGene(upstreamGeneIndex).getStart(), genome.getGene(upstreamGeneIndex).getStop()) + 1;
			if (!genome.getGene(upstreamGeneIndex).isORF()) start = genome.getGene(upstreamGeneIndex).getMaxCoordinate() + 1;
			break;
		    }
		    upstreamGeneIndex--;
		}
	    }
	    if (upstreamGeneIndex < 0) g1 = new Gene(1 + ".." + 1 + "\t" + g2.getStrand() + "\t" + "0" + "\t" + "-" + "\t" + "???" + "\t" + "???" + "\t" + "-" + "\t" + "-" + "\t" + "???", "ORF");
	    else g1 = genome.getGene(upstreamGeneIndex);

	    int IG_length = stop - start + 1;
	    if (IG_length <= 0) continue;  // No action necessary if there is no IG region.

	    int frontUTR_length_merged = -1;
	    int backUTR_length_merged = -1;
	    for (int i=0; i<conditions.size(); i++) {
		for (int j=0; j<conditions.get(i).numReplicates(); j++) {

		    Replicate r = conditions.get(i).getReplicate(j);

		    // Get UTR lengths
		    int frontUTR_length = getUTR_length(g1, start, stop, r, true);
		    int backUTR_length = getUTR_length(g2, start, stop, r, false);

		    // Distinguish overlapping UTRs
		    if ((frontUTR_length > 0) && (backUTR_length > 0) &&
			(Math.max(g1.getStart(), g1.getStop()) + frontUTR_length >= Math.min(g2.getStart(), g2.getStop()) - backUTR_length)) {

			// Determine mean of genes' expression
			double mean1 = r.getMeanOfRange(z, g1.getStart(), g1.getStop(), g1.getStrand());
			double mean2 = r.getMeanOfRange(z, g2.getStart(), g2.getStop(), g2.getStrand());
			if (unstranded) {
			    mean1 = r.getMeanOfRange(z, g1.getStart(), g1.getStop(), '?');
			    mean2 = r.getMeanOfRange(z, g2.getStart(), g2.getStop(), '?');
			}

			int overlapStart = Math.min(g2.getStart(), g2.getStop()) - backUTR_length;
			int overlapStop = Math.max(g1.getStart(), g1.getStop()) + frontUTR_length;
			//while ((overlapStart < Math.min(g2.getStart(), g2.getStop())) && (getPoissonPDF(r.getReads(overlapStart, g1.getStrand()), mean1) > getPoissonPDF(r.getReads(overlapStart, g2.getStrand()), mean2))) overlapStart++;
			if (!unstranded) {  // Strand specific
			    while ((overlapStart <= overlapStop) && (Math.abs(r.getReads(z, overlapStart, g1.getStrand()) - mean1) < Math.abs(r.getReads(z, overlapStart, g2.getStrand()) - mean2))) overlapStart++;
			    //while ((overlapStop > Math.max(g1.getStart(), g1.getStop())) && (getPoissonPDF(r.getReads(overlapStop, g1.getStrand()), mean1) < getPoissonPDF(r.getReads(overlapStop, g2.getStrand()), mean2))) overlapStop--;
			    while ((overlapStop >= overlapStart) && (Math.abs(r.getReads(z, overlapStop, g1.getStrand()) - mean1) > Math.abs(r.getReads(z, overlapStop, g2.getStrand()) - mean2))) overlapStop--;
			} else {  // Strand ambiguous
			    while ((overlapStart <= overlapStop) && (Math.abs(r.getReads(z, overlapStart, '?') - mean1) < Math.abs(r.getReads(z, overlapStart, '?') - mean2))) overlapStart++;
			    while ((overlapStop >= overlapStart) && (Math.abs(r.getReads(z, overlapStop, '?') - mean1) > Math.abs(r.getReads(z, overlapStop, '?') - mean2))) overlapStop--;
			}
			if (overlapStop - overlapStart < -1)
			    Rockhopper.output("Error - overlapping UTR region has invalid size.\n");
			else if (overlapStop - overlapStart == -1) {
			    // Do nothing. There is no overlap.
			} else {  // We still have overlapping UTR regions
			    double mean_of_overlap = r.getMeanOfRange(z, overlapStart, overlapStop, g1.getStrand());
			    if (unstranded) mean_of_overlap = r.getMeanOfRange(z, overlapStart, overlapStop, '?');
			    if (Math.abs(mean_of_overlap - mean1) <= Math.abs(mean_of_overlap - mean2))
				overlapStart = overlapStop + 1;
			    else overlapStop = overlapStart - 1;
			    //double prob1 = 1.0;
			    //double prob2 = 1.0;
			    //for (int y=overlapStart; y<=overlapStop; y++) {
			    //prob1 *= getPoissonPDF(r.getReads(y, g1.getStrand()), mean1);
			    //prob2 *= getPoissonPDF(r.getReads(y, g2.getStrand()), mean2);
			    //}
			    //if (prob1 >= prob2) overlapStart = overlapStop + 1;
			    //else overlapStop = overlapStart - 1;
			}
			frontUTR_length = overlapStart - (Math.max(g1.getStart(), g1.getStop()) + 1);
			backUTR_length = Math.min(g2.getStart(), g2.getStop()) - 1 - overlapStop;
		    }

		    // Merge UTRs from different experiments
		    //frontUTR_length_merged = Math.max(frontUTR_length_merged, frontUTR_length);
		    if (frontUTR_length_merged == -1) frontUTR_length_merged = frontUTR_length;
		    else if (frontUTR_length_merged == 0) frontUTR_length_merged = Math.max(frontUTR_length_merged, frontUTR_length);
		    else frontUTR_length_merged = Math.min(frontUTR_length_merged, frontUTR_length);
		    //backUTR_length_merged = Math.max(backUTR_length_merged, backUTR_length);
		    if (backUTR_length_merged == -1) backUTR_length_merged = backUTR_length;
		    else if (backUTR_length_merged == 0) backUTR_length_merged = Math.max(backUTR_length_merged, backUTR_length);
		    else backUTR_length_merged = Math.min(backUTR_length_merged, backUTR_length);
		}
	    }

	    //  Update transcription start/stop of genes (if genes are expressed)
	    if (frontUTR_length_merged >= 0) {  // Gene is expressed and has flanking IG region
		if (g1.getStrand() == '+') g1.setStopT(g1.getStop() + frontUTR_length_merged);
		if (g1.getStrand() == '-') g1.setStartT(g1.getStart() + frontUTR_length_merged);
	    }
	    if (backUTR_length_merged >= 0) {  // Gene is expressed and has flanking IG region
		if (g2.getStrand() == '+') g2.setStartT(g2.getStart() - backUTR_length_merged);
		if (g2.getStrand() == '-') g2.setStopT(g2.getStop() - backUTR_length_merged);
	    }
	}

	// Keep track of the number of 5'UTRs and 3'UTRs
	for (int x=0; x<genome.numGenes(); x++) {
	    Gene g = genome.getGene(x);
	    if (g.isORF() && (g.getStartT() > 0) && (g.getStrand() == '+') && (g.getStart()-g.getStartT() > 0))
		num5UTRs++;
	    if (g.isORF() && (g.getStartT() > 0) && (g.getStrand() == '-') && (g.getStartT()-g.getStart() > 0))
		num5UTRs++;
	    if (g.isORF() && (g.getStopT() > 0) && (g.getStrand() == '+') && (g.getStopT()-g.getStop() > 0))
		num3UTRs++;
	    if (g.isORF() && (g.getStopT() > 0) && (g.getStrand() == '-') && (g.getStop()-g.getStopT() > 0))
		num3UTRs++;
	}
    }

    /**
     * For each IG, identify ncRNAs based on
     * the expression data.
     */
    public void identifyRNAs() {

	// Generate coordinate map of genes
	String[] genesPlus = new String[genome.size()+1];
	String[] genesMinus = new String[genome.size()+1];
	for (int i=0; i<genesPlus.length; i++) {
	    genesPlus[i] = "";
	    genesMinus[i] = "";
	}
	for (int i=0; i<genome.numGenes(); i++) {
	    Gene g = genome.getGene(i);
	    int geneStart = Math.min(g.getStart(), g.getStop());
	    int geneStop = Math.max(g.getStart(), g.getStop());
	    if (!g.isORF()) {
		geneStart = g.getMinCoordinate();
		geneStop = g.getMaxCoordinate();
	    }
	    String[] genes = genesPlus;  // Plus strand
	    if (g.getStrand() == '-') genes = genesMinus;  // Minus strand
	    for (int j=geneStart; j<=geneStop; j++) genes[j] = g.getName();
	}

	ArrayList<Gene> rnaGenes = new ArrayList<Gene>();  // List of new predicted RNAs
	for (int x=0; x<genome.numGenes()+1; x++) {  // For each IG region
	    Gene g1 = null;
	    Gene g2 = null;

	    // Identify downstream gene and IG stop coordinate
	    if (x == genome.numGenes()) g2 = new Gene((genome.size()-1) + ".." + (genome.size()-1) + "\t" + genome.getGene(x-1).getStrand() + "\t" + "0" + "\t" + "-" + "\t" + "???" + "\t" + "???" + "\t" + "-" + "\t" + "-" + "\t" + "???", "ORF");
	    else g2 = genome.getGene(x);
	    if (!g2.isORF()) continue;  // Do not predict novel ncRNAs near known RNA genes
	    int stop = g2.getMinCoordinate() - 1;

	    // Identify upstream gene and IG start coordinate
	    int start = 1;
	    int upstreamGeneIndex = x-1;
	    while (upstreamGeneIndex >= 0) {
		if (unstranded) {  // Strand ambiguous
		    start = genome.getGene(upstreamGeneIndex).getMaxCoordinate() + 1;
		    break;
		} else {  // Strand specific
		    if (g2.getStrand() == genome.getGene(upstreamGeneIndex).getStrand()) {
			start = genome.getGene(upstreamGeneIndex).getMaxCoordinate() + 1;
			break;
		    }
		    upstreamGeneIndex--;
		}
	    }
	    if (upstreamGeneIndex < 0) g1 = new Gene(1 + ".." + 1 + "\t" + g2.getStrand() + "\t" + "0" + "\t" + "-" + "\t" + "???" + "\t" + "???" + "\t" + "-" + "\t" + "-" + "\t" + "???", "ORF");
	    else g1 = genome.getGene(upstreamGeneIndex);
	    if (!g1.isORF()) continue;  // Do not predict novel ncRNAs near known RNA genes

	    // Determine strand
	    char strand = '?';
	    if (unstranded) strand = '?';  // Strand ambiguous
	    else if ((g1.getStrand() == '+') || (g2.getStrand() == '+')) strand = '+';  // Plus strand
	    else if ((g1.getStrand() == '-') || (g2.getStrand() == '-')) strand = '-';  // Minus strand
	    else strand = '?';

	    // No action necessary if there is no IG region.
	    if (stop - start + 1 <= 0) {
		//output_IG_to_file(g1, g2, strand, new ArrayList<RNA>(), "IGs/Strep/");
		continue;
	    }

	    ArrayList<RNA> rnas = new ArrayList<RNA>();
	    for (int i=0; i<conditions.size(); i++) {
		for (int j=0; j<conditions.get(i).numReplicates(); j++) {
		    // Find novel transcripts on both strands separately
		    identifyNovelTranscriptsInIG(rnas, conditions.get(i).getReplicate(j), start, stop, strand);
		}
	    }

	    // Merge novel RNA transcripts (from different experiments and nearby novel RNAs)
	    ArrayList<RNA> merged_RNAs = new ArrayList<RNA>();
	    if (unstranded) {  // Strand ambiguous
		merge_RNAs(g1, g2, '?', rnas, merged_RNAs);
	    } else {  // Strand specific
		merge_RNAs(g1, g2, '+', rnas, merged_RNAs);
		merge_RNAs(g1, g2, '-', rnas, merged_RNAs);
	    }
	    //output_IG_to_file(g1, g2, strand, merged_RNAs, "IGs/Strep/");

	    /*
	    if (merged_RNAs.size() > 0) {
		int IG_start = Math.max(g1.getStart(), g1.getStop()) + 1;
		int IG_stop = Math.min(g2.getStart(), g2.getStop()) - 1;
		if (!g1.isORF()) IG_start = g1.getMaxCoordinate() + 1;
		if (!g2.isORF()) IG_stop = g2.getMinCoordinate() - 1;
		Rockhopper.output(IG_start + "_" + IG_stop + "\t" + g1.getName() + "\t" + g1.getStrand() + "\t" + g2.getName() + "\t" + g2.getStrand() + "\t" + rnas.size() + "\t" + merged_RNAs.size() + "\n");
	    }
	    */

	    // Convert RNAs to Genes
	    for (int i=0; i<merged_RNAs.size(); i++) {
		String product = getAntisenseAnnotation(merged_RNAs.get(i), genesPlus, genesMinus);
		if (product.indexOf("antisense") >= 0) this.numAntisenseRNAs++;
		rnaGenes.add(new Gene(merged_RNAs.get(i).start + ".." + merged_RNAs.get(i).stop + "\t" + merged_RNAs.get(i).strand + "\t" + "0" + "\t" + "-" + "\t" + "-" + "\t" + "predicted RNA" + "\t" + "-" + "\t" + "-" + "\t" + getAntisenseAnnotation(merged_RNAs.get(i), genesPlus, genesMinus), "RNA"));
	    }
	}
	this.numSenseRNAs = rnaGenes.size() - this.numAntisenseRNAs;

	// Compute expression of each predicted RNA
	for (int i=0; i<rnaGenes.size(); i++) {
	    Gene g = rnaGenes.get(i);
	    for (int j=0; j<conditions.size(); j++) {
		for (int k=0; k<conditions.get(j).numReplicates(); k++) {
		    Replicate r = conditions.get(j).getReplicate(k);
		    long readsForGene = r.getReadsInRange(z, g.getMinCoordinate(), g.getMaxCoordinate(), g.getStrand());
		    g.setRawCount(j, k, readsForGene);  // Set raw counts for gene
		    if (r.getAvgLengthReads() == 0) g.setRawCount_reads(j, k, 0);
		    else g.setRawCount_reads(j, k, readsForGene/r.getAvgLengthReads());
		    //g.setNormalizedCount(j, k, Condition.getAvgUpperQuartile(), r.getUpperQuartile());
		    g.setNormalizedCount(j, k, 100000.0, r.getUpperQuartile());
		}
	    }

	    // Compute mean and RPKM for each gene in each condition
	    g.computeExpression(conditions);

	    // Compute variance for each gene in each condition
	    g.computeVariance(conditions);
	}

	// Add transcripts to list of genes
	genome.addPredictedRNAs(rnaGenes);
    }



    /**************************************************
     **********   PRIVATE INSTANCE METHODS   **********
     **************************************************/

    private int getUTR_length(Gene g, int start, int stop, Replicate r, boolean isFront) {

	// If we do not have a real gene but merely a place holder, do not compute UTR.
	if (g.getStart() == g.getStop()) return -1;

	// If we have an RNA gene rather than an ORF, do not compute UTR.
	if (!g.isORF()) return -1;

	// Determine strand
	char strand = g.getStrand();
	if (unstranded) strand = '?';

	// Determine mean of gene's expression
	double mean = r.getMeanOfRange(z, g.getStart(), g.getStop(), strand);
	if (mean < r.getMinExpressionUTR()) return -1;  // Gene is not expressed

	// Compute expression distribution
	/*
	double thresh = 1.0;
	if (mean >= thresh * r.getMinExpressionUTR()) {
	    int start2 = g.getStart();
	    int stop2 = g.getStop();
	    if (stop2 < start2) {  // Swap
		int temp = start2;
		start2 = stop2;
		stop2 = temp;
	    }
	    for (int i=start2; i<=stop2; i++) {
		int x = (int)Math.round(20000*((r.getMeanOfRange(z, i, i, strand) - mean) / stdev));
		if (Math.abs(x) >= distribution.length/2) x = (distribution.length/2)*(x/Math.abs(x));
		distribution[x+distribution.length/2]++;
	    }
	}
	*/
	
	//int[] IG = new int[stop-start+1+WINDOW];
	//if (isFront) {  // Front UTR
	//   for (int i=start-WINDOW/2; i<=Math.min(stop+WINDOW/2, genome.size()-1); i++) IG[i-start+WINDOW/2] = r.getReads(i, g.getStrand());
	//} else {  // Back UTR
	//    for (int i=stop+WINDOW/2; i>=Math.max(start-WINDOW/2,1); i--) IG[stop+WINDOW/2-i] = r.getReads(i, g.getStrand());
	//}
	int[] IG = new int[stop-start+1];
	if (isFront) {  // Front UTR
	    for (int i=start; i<=Math.min(stop, genome.size()-1); i++) IG[i-start] = r.getReads(z, i, strand);
	} else {  // Back UTR
	    for (int i=stop; i>=Math.max(start, 1); i--) IG[stop-i] = r.getReads(z, i, strand);
	}

	//for (int i=0; i<IG.length-WINDOW; i++) {
	for (int i=0; i<IG.length; i++) {

	    double SCALE = 1.5;
	    if (IG[i] == 0) return i;
	    if (IG[i] >= SCALE*mean) continue;
	    if (IG[i] >= SCALE*r.getMinExpressionUTR()) continue;
	    if (SCALE*r.getBackgroundProb(IG[i]) > getPoissonPDF(IG[i], mean)) return i;
	    else continue;
	    /*
	    // Quick check to see if we can easily classify the window as IG or UTR
	    int numZeros = 0;
	    int numExpressed = 0;
	    for (int j=i; j<i+WINDOW; j++) {
		if (IG[j] == 0) numZeros++;
		if (IG[j] >= mean) numExpressed++;
	    }
	    if (numZeros >= WINDOW/2) return i;
	    if (numExpressed >= WINDOW/2) continue;

	    // More extensive check if probability is closer to IG or UTR
	    double probIG = 1.0;
	    double probUTR = 1.0;
	    for (int j=i; j<i+WINDOW; j++) {
		probIG *= r.getBackgroundProb(IG[j]);
		probUTR *= getPoissonPDF(IG[j], mean);
	    }
	    if (probIG >= probUTR) return i;
	    */
	}
	//return IG.length-WINDOW;
	return IG.length;
    }

    /**
     * For a given IG region, identifies novel transcripts on the specified strand.
     * Each novel transcript is added to the list of "rnas".
     */
    private void identifyNovelTranscriptsInIG(ArrayList<RNA> rnas, Replicate r, int start, int stop, char strand) {

	// Search for a transcript seed, i.e., n consecutive nucleotides above some threshold.
	// Then we extend this seed in both directions.
	//int n = WINDOW/2;
	double THRESHOLD = r.getMinExpressionRNA() * 2.0;
	int MIN_RNA_LENGTH = 10;

	//int[] IG = new int[stop-start+1+WINDOW];
	//int startReadsForIG = Math.max(start-WINDOW/2, 1);
	//for (int i=startReadsForIG; i<=Math.min(stop+WINDOW/2, genome.size()-1); i++) IG[i-startReadsForIG] = r.getReads(i, strand);
	int[] IG = new int[stop-start+1];
	int startReadsForIG = Math.max(start, 1);
	for (int i=startReadsForIG; i<=Math.min(stop, genome.size()-1); i++) IG[i-startReadsForIG] = r.getReads(z, i, strand);

	//int i = WINDOW/2;
	int i = 0;
	int startT = -1;
	int stopT = -1;
	//while (i < IG.length-WINDOW-n+1) {
	while (i < IG.length) {
	    if ((IG[i] >= THRESHOLD) && (startT == -1)) {  // Start of new possible seed
		startT = i;
		stopT = i;
	    } else if ((IG[i] >= THRESHOLD) && (startT >= 0)) {  // Within possible seed
		stopT = i;
	    } else if ((IG[i] < THRESHOLD) && (startT == -1)) {  // Within non-seed
		// Do nothing
	    } else if ((IG[i] < THRESHOLD) && (startT >= 0)) {  // Just ended possible seed
		//if (stopT - startT + 1 >= n) {  // We have a seed
		if (stopT - startT + 1 >= MIN_RNA_LENGTH) {  // We have a seed

		    // Extend the seed downstream
		    //boolean done = (stopT == IG.length-WINDOW);
		    boolean done = (stopT == IG.length-1);
		    while (!done) {
			double mean = r.getMeanOfRange(z, startT, stopT, strand);
			if (mean < THRESHOLD) break;

			if ((IG[stopT+1] >= THRESHOLD) || (r.getBackgroundProb(IG[stopT+1]) <= getPoissonPDF(IG[stopT+1], mean))) stopT++;
			else break;
			done = (stopT == IG.length-1);
			/*
			double probIG = 1.0;
			double probRNA = 1.0;
			for (int j=stopT+1-WINDOW/2+1; j<Math.min(stopT+1+WINDOW/2+1,IG.length); j++) {
			    probIG *= r.getBackgroundProb(IG[j]);
			    probRNA *= getPoissonPDF(IG[j], mean);
			}
			if (probIG >= probRNA) break;
			stopT++;
			if (stopT == IG.length-1-WINDOW) done = true;
			*/
		    }

		    // Extend the seed upstream
		    //done = (startT == WINDOW/2);
		    done = (startT == 0);
		    while (!done) {
			double mean = r.getMeanOfRange(z, startT, stopT, strand);
			if (mean < THRESHOLD) break;

			if ((IG[startT-1] >= THRESHOLD) || (r.getBackgroundProb(IG[startT-1]) <= getPoissonPDF(IG[startT-1], mean))) startT--;
			else break;
			done = (startT == 0);
			/*
			double probIG = 1.0;
			double probRNA = 1.0;
			for (int j=startT-1+WINDOW/2-1; j>Math.max(startT-1-WINDOW/2-1,-1); j--) {
			    probIG *= r.getBackgroundProb(IG[j]);
			    probRNA *= getPoissonPDF(IG[j], mean);
			}
			if (probIG >= probRNA) break;
			startT--;
			if (startT == WINDOW/2) done = true;
			*/
		    }

		    // Add transcript to list
		    rnas.add(new RNA(start + startT, start + stopT, strand));

		    i = stopT;
		    startT = -1;
		} else {  // Seed is too short, so it is not a seed.
		    startT = -1;
		}
	    } else {
		// Impossible to reach this case.
	    }
	    i++;
	}

	// Handle case where seed extends all the way to end of IG region
	if (startT >= 0) {  // Ended possible seed
	    //if (stopT - startT + 1 >= n) {  // We have a seed
	    if (stopT - startT + 1 >= MIN_RNA_LENGTH) {  // We have a seed

		// Extend the seed downstream
		//boolean done = (stopT == IG.length-WINDOW);
		boolean done = (stopT == IG.length-1);
		while (!done) {
		    double mean = r.getMeanOfRange(z, startT, stopT, strand);
		    if (mean < THRESHOLD) break;

		    if ((IG[stopT+1] >= THRESHOLD) || (r.getBackgroundProb(IG[stopT+1]) <= getPoissonPDF(IG[stopT+1], mean))) stopT++;
		    else break;
		    done = (stopT == IG.length-1);
		    /*
		    double probIG = 1.0;
		    double probRNA = 1.0;
		    for (int j=stopT+1-WINDOW/2+1; j<Math.min(stopT+1+WINDOW/2+1,IG.length); j++) {
			probIG *= r.getBackgroundProb(IG[j]);
			probRNA *= getPoissonPDF(IG[j], mean);
		    }
		    if (probIG >= probRNA) break;
		    stopT++;
		    if (stopT == IG.length-1-WINDOW) done = true;
		    */
		}

		// Extend the seed upstream
		//done = (startT == WINDOW/2);
		done = (startT == 0);
		while (!done) {
		    double mean = r.getMeanOfRange(z, startT, stopT, strand);
		    if (mean < THRESHOLD) break;

		    if ((IG[startT-1] >= THRESHOLD) || (r.getBackgroundProb(IG[startT-1]) <= getPoissonPDF(IG[startT-1], mean))) startT--;
		    else break;
		    done = (startT == 0);
		    /*		    
		    double probIG = 1.0;
		    double probRNA = 1.0;
		    for (int j=startT-1+WINDOW/2-1; j>Math.max(startT-1-WINDOW/2-1,-1); j--) {
			probIG *= r.getBackgroundProb(IG[j]);
			probRNA *= getPoissonPDF(IG[j], mean);
		    }
		    if (probIG >= probRNA) break;
		    startT--;
		    if (startT == WINDOW/2) done = true;
		    */
		}

		// Add transcript to list
		rnas.add(new RNA(start + startT, start + stopT, strand));
	    }
	}
    }

    /**
     * Given a list of RNAs, merges overlapping RNAs from different experiments
     * and RNAs that are close in proximity. Adds the merged RNAs to the 
     * list "merged_RNAs".
     */ 
    private void merge_RNAs(Gene g1, Gene g2, char strand, ArrayList<RNA> rnas, ArrayList<RNA> merged_RNAs) {
	if (rnas.size() == 0) return;

	// Get coordinates of genes flanking the IG region
	int start1 = g1.getStart();
	int stop1 = g1.getStop();
	if (!g1.isORF()) {
	    start1 = g1.getStartT();
	    stop1 = g1.getStopT();
	}
	if (start1 > stop1) {  // Swap
	    int temp = start1;
	    start1 = stop1;
	    stop1 = temp;
	}
	int start2 = g2.getStart();
	int stop2 = g2.getStop();
	if (!g2.isORF()) {
	    start2 = g2.getStartT();
	    stop2 = g2.getStopT();
	}
	if (start2 > stop2) {  // Swap
	    int temp = start2;
	    start2 = stop2;
	    stop2 = temp;
	}

	int IG_start = stop1 + 1;
	int IG_stop = start2 - 1;
	int IG_length = IG_stop - IG_start + 1;
	if (IG_length <= 0) return;  // IG region has length 0. Do not output it.

	// For each nucleotide in an IG region, indicates if a RNA is predicted
	// for the nucleotide or not.
	boolean[] IG = new boolean[IG_length];
	for (int i=0; i<rnas.size(); i++) {
	    if ((rnas.get(i).strand == strand) || unstranded) {
		for (int j=rnas.get(i).start; j<=rnas.get(i).stop; j++) {
		    if (j < IG_start) return;
		    IG[j-IG_start] = true;
		}
	    }
	}

	// Determine set of merged RNAs
	int PROXIMITY = 50;  // RNAs within this many NTs of each other are merged
	int start = -1;
	int stop = -1;
	for (int i=0; i<IG.length; i++) {
	    if (IG[i] && (start == -1)) {  // Start of RNA
		start = i;
		stop = i;
	    } else if (IG[i] && (start >= 0)) {  // Within RNA
		stop = i;
	    } else if (!IG[i] && (start == -1)) {  // Within non-RNA
		// do nothing
	    } else if (!IG[i] && (start >= 0)) {  // RNA just ended
		for (int j=i; j<Math.min(i+PROXIMITY, IG.length); j++) {
		    if (IG[j]) stop = j;
		}
		if (stop > i-1) i = stop;  // Continue the loop from the end of the RNA
		else {
		    // Ignore predicted RNAs that abut annotated RNAs
		    if (((g1.isORF()) || (start > 0)) && ((g2.isORF()) || (stop < IG.length-1)))
			merged_RNAs.add(new RNA(IG_start + start, IG_start + stop, strand));
		    start = -1;
		    stop = -1;
		}
	    } else {
		Rockhopper.output("Error - this case should be unreachable!\n");
	    }
	}
	// Check if IG region ends with RNA
	if (start >= 0) {
	    int previousStop = stop;
	    for (int j=previousStop+1; j<Math.min(previousStop+1+PROXIMITY, IG.length); j++) {
		if (IG[j]) stop = j;
	    }
	    // Ignore predicted RNAs that abut annotated RNAs
	    if (((g1.isORF()) || (start > 0)) && ((g2.isORF()) || (stop < IG.length-1)))
		merged_RNAs.add(new RNA(IG_start + start, IG_start + stop, strand));
	}
    }

    /**
     * Given a predicted RNA, determines if it is antisense to any annotated genes.
     * Returns "-" if it is not antisense to any genes.
     * Otherwise, returns "antisense" followed by a list of genes it is antisense to.
     */
    private String getAntisenseAnnotation(RNA rna, String[] genesPlus, String[] genesMinus) {
	HashMap<String, String> genes = new HashMap<String, String>();
	for (int i=rna.start; i<=rna.stop; i++) {
	    if (unstranded || (rna.strand == '+')) genes.put(genesMinus[i], "");
	    if (unstranded || (rna.strand == '-')) genes.put(genesPlus[i], "");
	}

	String s = "";
	for (String key : genes.keySet()) {
	    if (key.length() > 0) s += " " + key;
	}

	if (s.length() == 0) return "-";
	else return "antisense:" + s;
    }

    /**
     * Output reads in IG region along with UTR and RNA annotation of region.
     */
    private void output_IG_to_file(Gene g1, Gene g2, char strand, ArrayList<RNA> rnas, String DIR) {

	// Get coordinates of genes flanking the IG region
	int start1 = g1.getStart();
	int stop1 = g1.getStop();
	if (!g1.isORF()) {
	    start1 = g1.getStartT();
	    stop1 = g1.getStopT();
	}
	if (start1 > stop1) {  // Swap
	    int temp = start1;
	    start1 = stop1;
	    stop1 = temp;
	}
	int start2 = g2.getStart();
	int stop2 = g2.getStop();
	if (!g2.isORF()) {
	    start2 = g2.getStartT();
	    stop2 = g2.getStopT();
	}
	if (start2 > stop2) {  // Swap
	    int temp = start2;
	    start2 = stop2;
	    stop2 = temp;
	}

	int IG_start = stop1 + 1;
	int IG_stop = start2 - 1;
	int IG_length = IG_stop - IG_start + 1;
	if (IG_length <= 0) return;  // IG region has length 0. Do not output it.

	int mean1 = 0;
	int mean2 = 0;
	for (int x=0; x<conditions.size(); x++)
	    for (int y=0; y<conditions.get(x).numReplicates(); y++) {
		if (!unstranded) {  // Strand specific
		    mean1 = Math.max(mean1, (int)conditions.get(x).getReplicate(y).getMeanOfRange(z, start1, stop1, g1.getStrand()));
		    mean2 = Math.max(mean2, (int)conditions.get(x).getReplicate(y).getMeanOfRange(z, start2, stop2, g2.getStrand()));
		} else {  // Starnd ambiguous
		    mean1 = Math.max(mean1, (int)conditions.get(x).getReplicate(y).getMeanOfRange(z, start1, stop1, '?'));
		    mean2 = Math.max(mean2, (int)conditions.get(x).getReplicate(y).getMeanOfRange(z, start2, stop2, '?'));
		}
	    }
	int[] IG = new int[IG_length];
	for (int i=IG_start; i<=IG_stop; i++) {
	    int maxRead = 0;
	    for (int x=0; x<conditions.size(); x++)
		for (int y=0; y<conditions.get(x).numReplicates(); y++)
		    maxRead = Math.max(maxRead, conditions.get(x).getReplicate(y).getReads(z, i, strand));
	    IG[i-IG_start] = maxRead;
	}
	StringBuilder annotation = getAnnotationOfIG(IG, g1, g2, strand, rnas);

	try {

	    // Create directory (if it doesn't exist) to write files to
	    if (DIR.charAt(DIR.length()-1) != '/') DIR += '/';
	    File dir = new File(DIR);
	    if (!dir.exists()) dir.mkdir();

	    String fileName = IG_start + "_" + IG_stop + ".txt";
	    PrintWriter writer = new PrintWriter(new File(DIR + fileName));
	    writer.println(g1.getName() + "\t" + g1.getStrand() + "\t" + mean1);
	    writer.println(g2.getName() + "\t" + g2.getStrand() + "\t" + mean2);
	    writer.println("\n");
	    writer.println(annotation.toString());
	    writer.close();
	} catch (FileNotFoundException ex) {
	    Rockhopper.output("Error - could not output file.\n");
	    System.exit(0);
	}
    }

    /**
     * Output reads in IG region and UTR of specified gene.
     */
    /*
    private void output_IG_to_file(Gene g1, Gene g2, int frontUTR_length, int backUTR_length, int start, int stop, Replicate r, String DIR) {

	int mean1 = (int)r.getMeanOfRange(g1.getStart(), g1.getStop(), g1.getStrand());
	int mean2 = (int)r.getMeanOfRange(g2.getStart(), g2.getStop(), g2.getStrand());
	int[] IG1 = new int[stop-start+1];
	int[] IG2 = null;
	if ((g1.getStrand() != '+') && (g1.getStrand() != '-')) {  // Placeholder gene
	    for (int i=start; i<=stop; i++) IG1[i-start] = r.getReads(i, g2.getStrand());
	} else if ((g2.getStrand() != '+') && (g2.getStrand() != '-')) {  // Placeholder gene
	    for (int i=start; i<=stop; i++) IG1[i-start] = r.getReads(i, g1.getStrand());
	} else if (g1.getStrand() == g2.getStrand()) {  // Same strand
	    for (int i=start; i<=stop; i++) IG1[i-start] = r.getReads(i, g1.getStrand());
	} else {  // Different strands
	    for (int i=start; i<=stop; i++) IG1[i-start] = r.getReads(i, g1.getStrand());
	    IG2 = new int[stop-start+1];
	    for (int i=start; i<=stop; i++) IG2[i-start] = r.getReads(i, g2.getStrand());
	}
	StringBuilder annotation1 = getAnnotationOfIG(IG1, g1, g2, frontUTR_length, backUTR_length);
	StringBuilder annotation2 = getAnnotationOfIG(IG2, g1, g2, frontUTR_length, backUTR_length);

	try {

	    // Create directory (if it doesn't exist) to write files to
	    if (DIR.charAt(DIR.length()-1) != '/') DIR += '/';
	    File dir = new File(DIR);
	    if (!dir.exists()) dir.mkdir();
	    dir = new File(DIR + r.getName() + "/");
	    if (!dir.exists()) dir.mkdir();
	    DIR += r.getName() + "/";

	    String fileName = start + "_" + stop + ".txt";
	    PrintWriter writer = new PrintWriter(new File(DIR + fileName));
	    writer.println(g1.getName() + "\t" + g1.getStrand() + "\t" + mean1);
	    writer.println(g2.getName() + "\t" + g2.getStrand() + "\t" + mean2);
	    writer.println("\n");
	    writer.println(annotation1.toString());
	    if (annotation2 != null) writer.println("\n" + annotation2.toString());
	    writer.close();
	} catch (FileNotFoundException ex) {
	    Rockhopper.output("Error - could not output file.\n");
	    System.exit(0);
	}
    }
    */

    /**
     * Return a StringBuilder representation of an IG region consisting 
     * of reads and an annotation (5'UTR and 3'UTR and RNA) of those reads.
     */
    private StringBuilder getAnnotationOfIG(int[] IG, Gene g1, Gene g2, char strand, ArrayList<RNA> rnas) {
	int FASTA = 25;
	StringBuilder annotation = new StringBuilder(IG.length);
	for (int i=0; i<IG.length; i++) annotation.append(".");

	// Front UTR
	if ((g1.getStrand() == strand) || (this.unstranded)) {
	    int frontUTR_count = g1.getMaxCoordinate() - Math.max(g1.getStart(), g1.getStop());
	    if (!g1.isORF()) frontUTR_count = 0;
	    char UTR_char = '3';  // Plus strand
	    if (g1.getStrand() == '-') UTR_char = '5';  // Minus strand
	    for (int i=0; i<frontUTR_count; i++) annotation.setCharAt(i, UTR_char);
	}

	// Back UTR
	if ((g2.getStrand() == strand) || (this.unstranded)) {
	    int backUTR_count = Math.min(g2.getStart(), g2.getStop()) - g2.getMinCoordinate();
	    if (!g2.isORF()) backUTR_count = 0;
	    char UTR_char = '5';  // Plus strand
	    if (g2.getStrand() == '-') UTR_char = '3';  // Minus strand
	    for (int i=0; i<backUTR_count; i++) annotation.setCharAt(annotation.length()-1-i, UTR_char);
	}

	// RNAs
	int IG_start = Math.max(g1.getStart(), g1.getStop()) + 1;
	if (!g1.isORF()) IG_start = g1.getMaxCoordinate() + 1;
	int IG_stop = Math.min(g2.getStart(), g2.getStop()) - 1;
	if (!g2.isORF()) IG_stop = g2.getMinCoordinate() - 1;
	for (int i=0; i<rnas.size(); i++) {
	    RNA rna = rnas.get(i);
	    if ((rna.strand == strand) || (this.unstranded)) {
		for (int j=rna.start; j<=rna.stop; j++)
		    annotation.setCharAt(j - IG_start, 'R');
	    }
	}

	StringBuilder IG_annotation = new StringBuilder(IG.length);
	int i = 0;
	for (i=0; i<IG.length; i++) {
	    String num = "" + IG[i];
	    if (IG[i] < 10) num  = " " + num + " ";
	    else if (IG[i] < 100) num = num + " ";
	    else num = num.charAt(0) + "e" + (num.length()-1);
	    IG_annotation.append(" " + num);
	    if ((i+1) % FASTA == 0) {
		IG_annotation.append("\n");
		for (int j=0; j<FASTA; j++) {
		    IG_annotation.append(" " + " " + annotation.charAt((i+1) - FASTA + j) + " ");
		}
		IG_annotation.append("\n\n");
	    }
	}
	if (i % FASTA != 0) {
	    IG_annotation.append("\n");
	    for (int j=0; j<(i%FASTA); j++) {
		IG_annotation.append(" " + " " + annotation.charAt(i - (i%FASTA) + j) + " ");
	    }
	    IG_annotation.append("\n\n");
	}
	return IG_annotation;
    }

    /**
     * Return a StringBuilder representation of an IG region consisting 
     * of reads and an annotation (5'UTR and 3'UTR) of those reads.
     */
    /*
    private StringBuilder getAnnotationOfIG(int[] IG, Gene g1, Gene g2, int frontUTR_length, int backUTR_length) {
	int FASTA = 25;
	if (IG == null) return null;
	StringBuilder annotation = new StringBuilder();
	int frontUTR_count = 0;
	int backUTR_count = 0;
	int i = 0;
	for (i=0; i<IG.length; i++) {
	    String num = "" + IG[i];
	    if (IG[i] < 10) num  = " " + num + " ";
	    else if (IG[i] < 100) num = num + " ";
	    else num = num.charAt(0) + "e" + (num.length()-1);
	    annotation.append(" " + num);
	    if ((i+1) % FASTA == 0) {
		annotation.append("\n");
		for (int j=0; j<FASTA; j++) {
		    if (frontUTR_count < frontUTR_length) {
			if (g1.getStrand() == '+') annotation.append(" " + " 3 ");
			else annotation.append(" " + " 5 ");
			frontUTR_count++;
		    }
		    else if (i+1-FASTA+j >= IG.length - backUTR_length) {
			if (g2.getStrand() == '+') annotation.append(" " + " 5 ");
			else annotation.append(" " + " 3 ");
		    } else annotation.append(" " + " . ");
		}
		annotation.append("\n\n");
	    }
	}
	if ((i+1) % FASTA != 0) annotation.append("\n");
	for (int j=0; j<(IG.length%FASTA); j++) {
	    if (frontUTR_count < frontUTR_length) {
		if (g1.getStrand() == '+') annotation.append(" " + " 3 ");
		else annotation.append(" " + " 5 ");
		frontUTR_count++;
	    }
	    else if (IG.length - (IG.length%FASTA) + j >= IG.length - backUTR_length) {
		if (g2.getStrand() == '+') annotation.append(" " + " 5 ");
		else annotation.append(" " + " 3 ");
	    } else annotation.append(" " + " . ");
	}
	annotation.append("\n\n");
	return annotation;
    }
    */



    /***********************************************
     **********   PRIVATE CLASS METHODS   **********
     ***********************************************/

    /**
     * Returns the PDF value at x for the specified Poisson distribution.
     */
    private static double getPoissonPDF(int x, double lambda) {
	double result = Math.exp(-lambda);
	int k = x;
	while (k >= 1) {
	    result *= lambda/k;
	    k--;
	}
	return result;
    }

}





/***********************************
 **********   RNA CLASS   **********
 ***********************************/

class RNA {

    public int start;
    public int stop;
    public char strand;

    public RNA(int start, int stop, char strand) {
	this.start = start;
	this.stop = stop;
	this.strand = strand;
    }

    public String toString() {
	return start + "\t" + stop + "\t" + strand;
    }

}

