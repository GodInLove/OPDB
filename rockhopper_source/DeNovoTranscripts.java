/*
 * Copyright 2014 Brian Tjaden
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
import java.util.Random;
import java.util.concurrent.atomic.AtomicIntegerArray;

/**
 * An instance of the DeNovoTranscripts class represents
 * a collection of de novo assembled transcripts.
 */
public class DeNovoTranscripts {

    /********************************************
     **********   Instance Variables   **********
     ********************************************/

    private ArrayList<DeNovoTranscript> transcripts;
    private ArrayList<ArrayList<Long>> upperQuartiles;
    private ArrayList<Integer> minDiffExpressionLevels;  // Min diff expression level per condition



    /*****************************************
     **********   Class Variables   **********
     *****************************************/

    public static int[] avgLengthOfReads;  // Avg length of mapping reads in each file
    public static int[] numReads;          // Total reads in each file
    public static int[] numMappingReads;   // Reads mapping to transcripts in each file
    public static ArrayList<ArrayList<Long>> totalReads;  // Total number of nt reads in each replicate
    private static Random rand = new Random();  // Random number generator



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    public DeNovoTranscripts(DeNovoIndex bwtIndex) {
	DeNovoTranscripts.avgLengthOfReads = bwtIndex.getAvgLengthOfReads();
	DeNovoTranscripts.numReads = bwtIndex.getNumReads();
	DeNovoTranscripts.numMappingReads = bwtIndex.getNumMappingReads();
	DeNovoTranscripts.totalReads = bwtIndex.getTotalReads();
	transcripts = determineTranscripts(bwtIndex.sequence, bwtIndex.readCounts);
	computeExpression();
    }

    /**
     * Used when reading transcripts in from compressed file.
     */
    public DeNovoTranscripts(ArrayList<DeNovoTranscript> transcripts) {
	this.transcripts = transcripts;
    }



    /*************************************************
     **********   Public Instance Methods   **********
     *************************************************/

    public int getNumTranscripts() {
	return transcripts.size();
    }

    public DeNovoTranscript getTranscript(int i) {
	return transcripts.get(i);
    }

    /**
     * Returns the average transcript length.
     */
    public int getAverageTranscriptLength() {
	long length = 0;
	if (transcripts.size() == 0) return 0;
	for (DeNovoTranscript transcript : transcripts) length += transcript.length();
	return (int)(length/transcripts.size());
    }

    /**
     * Returns the median transcript length.
     */
    public int getMedianTranscriptLength() {
	ArrayList<Long> a = new ArrayList<Long>(transcripts.size());
	for (DeNovoTranscript transcript : transcripts) a.add((long)transcript.length());
	return (int)Misc.select_Long(a, a.size()/2);
    }

    /**
     * Returns the sum of all transcript lengths.
     */
    public long getTotalAssembledBases() {
	long sum = 0;
	for (DeNovoTranscript transcript : transcripts) sum += transcript.length();
	return sum;
    }

    public String getTranscriptSequences() {
	StringBuilder sb = new StringBuilder();
	for (int z=0; z<transcripts.size(); z++) sb.append(transcripts.get(z).sequence() + "\n");
	return sb.toString();
    }

    public void computeDifferentialExpression() {
	computeVarianceAndLowess();
	for (DeNovoTranscript transcript : transcripts) transcript.computeDifferentialExpression();
	correctPvalues();
    }

    /**
     * The background probabilty decreases as the number of reads increases.
     * We use this fact to estimate a minimum level of expression necessary
     * for us to be able to compute a p-value of differential expression.
     * For a given Condition, we ensure that all Replicates meet the 
     * expression threshold.
     * Currently, the expression threshold is set (based on anecdote and
     * experience) to 0.005.
     */
    public void setMinDiffExpressionLevels(DeNovoIndex bwtIndex) {
	double THRESHOLD = 0.005;
	String allSequences = bwtIndex.sequence;
	AtomicIntegerArray[] readCounts = bwtIndex.readCounts;
	int backgroundLength = 2;
	int linearIndex = 0;
	minDiffExpressionLevels = new ArrayList<Integer>(Assembler.conditionFiles.size());
	for (int i=0; i<Assembler.conditionFiles.size(); i++) {  // For each condition
	    String[] files = Assembler.conditionFiles.get(i).split(",");
	    int[][] background = new int[files.length][backgroundLength];
	    double[] backgroundParameter = new double[files.length];  // Parameter of geometric distribution of background reads. One per replicate in current condition.
	    for (int j=0; j<files.length; j++) {  // For each replicate

		// Populate background for current replicate
		for (int z=0; z<allSequences.length()-1; z++) {  // We go to length-1 to ignore final '$'
		    if (readCounts[linearIndex].get(z) < background[j].length) background[j][readCounts[linearIndex].get(z)]++;
		    else background[j][background[j].length-1]++;
		}

		// Compute background parameter for current replicate
		for (int k=0; k<backgroundLength-1; k++) backgroundParameter[j] += background[j][k];
		backgroundParameter[j] /= (3.0*1.0*allSequences.length());  // We use 3.0*1.0 since there is only one strand

		linearIndex++;
	    }

	    // Compute minDiffExpressionLevel for current condition
	    int minDiffExpressionLevel = 0;
	    while (true) {
		boolean foundThreshold = true;
		for (int j=0; j<files.length; j++) {
		    double backgroundProb = Math.pow(1.0 - backgroundParameter[j], minDiffExpressionLevel) * backgroundParameter[j];  // Based on geometric distribution
		    if (backgroundProb >= THRESHOLD)
			foundThreshold = false;
		}
		if (foundThreshold) break;
		minDiffExpressionLevel++;
	    }
	    minDiffExpressionLevels.add(minDiffExpressionLevel);
	}
    }

    public String toString(String[] labels) {

	// Include Header
	StringBuilder sb = new StringBuilder();
	sb.append("Sequence" + "\t" + "Length");
	for (int i=0; i<Assembler.conditionFiles.size(); i++) {
	    String conditionName = "" + (i+1);
	    if ((labels != null) && (labels.length == Assembler.conditionFiles.size())) conditionName = labels[i];
	    if (Assembler.verbose) {  // Verbose output
		String[] files = Assembler.conditionFiles.get(i).split(",");
		if (files.length == 1) {
		    sb.append("\t" + "Raw Counts " + conditionName);
		    sb.append("\t" + "Normalized Counts " + conditionName);
		} else {
		    for (int j=0; j<files.length; j++) sb.append("\t" + "Raw Counts " + conditionName + " Replicate " + (j+1));
		    for (int j=0; j<files.length; j++) sb.append("\t" + "Normalized Counts " + conditionName + " Replicate " + (j+1));
		}
		sb.append("\t" + "RPKM " + conditionName);
	    }
	    sb.append("\t" + "Expression " + conditionName);
	}
	for (int x=0; x<Assembler.conditionFiles.size()-1; x++) {  // First in pair
	    String conditionName1 = "" + (x+1);
	    if ((labels != null) && (labels.length == Assembler.conditionFiles.size())) conditionName1 = labels[x];
	    for (int y=x+1; y<Assembler.conditionFiles.size(); y++) {  // Second in pair
		String conditionName2 = "" + (y+1);
		if ((labels != null) && (labels.length == Assembler.conditionFiles.size())) conditionName2 = labels[y];

		if (Assembler.verbose)  // Verbose output
		    sb.append("\t" + "pValue " + conditionName1 + " vs " + conditionName2);
		sb.append("\t" + "qValue " + conditionName1 + " vs " + conditionName2);
	    }
	}
	sb.append("\n");

	// Output each transcript
	for (int z=0; z<transcripts.size(); z++) sb.append(transcripts.get(z) + "\n");
	return sb.toString();
    }



    /**************************************************
     **********   Private Instance Methods   **********
     **************************************************/

    /**
     * Determine transcripts based on DeNovoIndex.
     */
    private ArrayList<DeNovoTranscript> determineTranscripts(String allSequences, AtomicIntegerArray[] readCounts) {
	ArrayList<DeNovoTranscript> transcripts = new ArrayList<DeNovoTranscript>();
	int start = -1;      // -1 is sentinel indicating no transcript
	int end = -1;        // inclusive
	long[] readCounts_nts = new long[readCounts.length];
	for (int i=0; i<allSequences.length()-1; i++) {  // We go to length-1 to ignore final '$'
	    if (allSequences.charAt(i) == '^') {  // End of transcript.
		if ((start >= 0) && (end-start+1 >= Assembler.minTranscriptLength)) transcripts.add(new DeNovoTranscript(allSequences.substring(start, end+1), readCounts_nts));
		start = -1;
		end = -1;
		for (int j=0; j<readCounts.length; j++) readCounts_nts[j] = 0;
	    } else {  // Within transcript
		if ((start < 0) && (getReadCountAtIndex(readCounts, i) >= Assembler.MIN_READS_MAPPING)) {  // Start of new expressed transcript
		    start = i;
		    for (int j=0; j<readCounts.length; j++) readCounts_nts[j] += readCounts[j].get(i);
		}
		if (getReadCountAtIndex(readCounts, i) >= Assembler.MIN_READS_MAPPING) {  // Within expressed transcript
		    end = i;
		    for (int j=0; j<readCounts.length; j++) readCounts_nts[j] += readCounts[j].get(i);
		} else {  // It is possible expressed transcript just ended
		    if ((start >= 0) && (end-start+1 >= Assembler.minTranscriptLength)) transcripts.add(new DeNovoTranscript(allSequences.substring(start, end+1), readCounts_nts));
		    start = -1;
		    end = -1;
		    for (int j=0; j<readCounts.length; j++) readCounts_nts[j] = 0;
		}
	    }
	}
	return transcripts;
    }

    /**
     * For each transcript, compute its normalized expression, RPKM, and
     * mean expression in each condition and/or replicate.
     */
    private void computeExpression() {
	upperQuartiles = new ArrayList<ArrayList<Long>>(Assembler.conditionFiles.size());
	for (int i=0; i<Assembler.conditionFiles.size(); i++) {
	    String[] files = Assembler.conditionFiles.get(i).split(",");
	    upperQuartiles.add(new ArrayList<Long>(files.length));
	    for (int j=0; j<files.length; j++) {
		ArrayList<Long> expressions = new ArrayList<Long>(transcripts.size());
		for (int k=0; k<transcripts.size(); k++) expressions.add(transcripts.get(k).getRawCountNTs(i, j));
		long upperQuartile = Misc.select_Long(expressions, (int)(0.75*transcripts.size()));
		upperQuartiles.get(i).add(upperQuartile);
		for (int k=0; k<transcripts.size(); k++) transcripts.get(k).setNormalizedCount(i, j, 100000.0, upperQuartile);
	    }
	    for (int k=0; k<transcripts.size(); k++) transcripts.get(k).computeExpression(i);
	}
    }

    /**
     * For each transcripts, compute its expression variance and lowess
     * in each condition.
     */
    private void computeVarianceAndLowess() {
	// Compute variance
	for (int k=0; k<transcripts.size(); k++) transcripts.get(k).computeVariance(Assembler.conditionFiles.size());

	// Compute Lowess
	for (int x=0; x<Assembler.conditionFiles.size(); x++) {  // First of pair
	    String[] files = Assembler.conditionFiles.get(x).split(",");
	    for (int y=0; y<Assembler.conditionFiles.size(); y++) {  // Second of pair

		if (x == y) continue;  // No need to compute lowess for self
		/*
		if ((x > 0) && (files.length > 1)) {  // We have replicates
		    for (DeNovoTranscript transcript : transcripts)
			transcript.lowess[x][y] = transcript.lowess[0][y];
		    continue;
		}
		*/
		/* (x,y) lowess is NOT the same as (y,x) lowess
		String[] files2 = Assembler.conditionFiles.get(y).split(",");
		if ((x > y) && (files.length == 1) && (files2.length == 1)) {  // Already computed lowess[y][x], which is same as lowess[x][y]
		    for (DeNovoTranscript transcript : transcripts)
			transcript.lowess[x][y] = transcript.lowess[y][x];
		    continue;
		}
		*/

		double b = 0.0;  // Bias correction term
		for (int j=0; j<files.length; j++) {
		    b += 100000.0 / (double)upperQuartiles.get(x).get(j);
		}
		b /= files.length;

		// Create list of gene expressions and list of gene variances
		ArrayList<Long> expression = new ArrayList<Long>();
		ArrayList<Long> variance = new ArrayList<Long>();
		for (DeNovoTranscript transcript : transcripts) {
		    expression.add(transcript.getMean(x));
		    variance.add(transcript.variances[x][y]);
		}

		// Perform Lowess computation
		ArrayList<Long> lowessVariance = Lowess.lowess(expression, variance);

		// Uncomment to output data for Lowess graph to StdOut
		//System.out.println("Mean" + "\t" + "Variance" + "\t" + "Lowess");
		//for (int k=0; k<lowessVariance.size(); k++) System.out.println(expression.get(k) + "\t" + variance.get(k) + "\t" + ((long)(lowessVariance.get(k) - expression.get(k)*b)));

		// Assign each gene its lowess variances (after subtracting bias correction term)
		for (int k=0; k<transcripts.size(); k++) {
		    transcripts.get(k).lowess[x][y] = (long)(lowessVariance.get(k) - transcripts.get(k).getMean(x) * b);
		}
	    }
	}
    }

    /**
     * Computes q-values for each transcript, i.e., corrected p-values,
     * using Benjamini Hochberg correction.
     */
    private void correctPvalues() {
	if (transcripts.size() == 0) return;
	int pValue_index = 0;
	for (int x=0; x<transcripts.get(0).variances.length-1; x++) {  // First in pair
	    for (int y=x+1; y<transcripts.get(0).variances.length; y++) {  // Second in pair

		int[] indices = new int[transcripts.size()];
		double[] pvalues = new double[transcripts.size()];
		for (int j=0; j<transcripts.size(); j++) {
		    DeNovoTranscript transcript = transcripts.get(j);
		    pvalues[j] = transcript.getPvalue(pValue_index);
		    indices[j] = j;

		    // Check if there is too little expression to compute a p-value
		    double e1 = transcript.getMean(x);
		    double e2 = transcript.getMean(y);
		    e1 /= transcript.length();
		    e2 /= transcript.length();
		    if ((e1 < minDiffExpressionLevels.get(x)) && (e2 < minDiffExpressionLevels.get(y)))
			pvalues[j] = 1.0;
		}
		mergesort(pvalues, indices, 0, transcripts.size()-1);
		double previous_BH_value = 0.0;
		for (int k=0; k<pvalues.length; k++) {
		    double BH_value = pvalues[k] * transcripts.size() / (k+1);
		    BH_value = Math.min(BH_value, 1.0);  // Disallow values greater than 1

		    // To preserve monotonicity, we take maximum to ensure we don't
		    // generate a value less than the previous value.
		    BH_value = Math.max(BH_value, previous_BH_value);
		    previous_BH_value = BH_value;
		    transcripts.get(indices[k]).setQvalue(pValue_index, BH_value);
		}
		pValue_index++;
	    }
	}
    }

    /**
     * For each condition, determines the other condition that is most similar,
     * i.e., its "partner". If there are no replicate experiments then the
     * partner of each condition is used as a surrogate replicate.
     * Similarity between two conditions is measured by Pearson correlation
     * coefficient of normalized gene expression.
     * This method sets the "partner" of each condition.
     */
    private void identifySimilarConditions() {
	// Unused.
    }

    /**
     * Returns the number of reads mapping to nucleotide index i
     * summed over all sequencing read files.
     */
    private int getReadCountAtIndex(AtomicIntegerArray[] readCounts, int index) {
	int sum = 0;
	for (int i=0; i<readCounts.length; i++) sum += readCounts[i].get(index);
	return sum;
    }

    /**
     * Mergesort parallel arrays "a" and "b" based on values in "a"
     */
    private static void mergesort(double[] a, int[] b, int lo, int hi) {
	if (lo < hi) {
	    int q = (lo+hi)/2;
	    mergesort(a, b, lo, q);
	    mergesort(a, b, q+1, hi);
	    merge(a, b, lo, q, hi);
	}
    }

    /**
     * Mergesort helper method
     */
    private static void merge(double[] a, int[] b, int lo, int q, int hi) {
	double[] a1 = new double[q-lo+1];
	double[] a2 = new double[hi-q];
	int[] b1 = new int[q-lo+1];
	int[] b2 = new int[hi-q];
	for (int i=0; i<a1.length; i++) {
	    a1[i] = a[lo+i];
	    b1[i] = b[lo+i];
	}
	for (int j=0; j<a2.length; j++) {
	    a2[j] = a[q+1+j];
	    b2[j] = b[q+1+j];
	}
	int i=0;  // Index for first half arrays
	int j=0;  // Index for second half arrays
	for (int k=lo; k<=hi; k++) {
	    if (i >= a1.length) {
		a[k] = a2[j];
		b[k] = b2[j];
		j++;
	    } else if (j >= a2.length) {
		a[k] = a1[i];
		b[k] = b1[i];
		i++;
	    } else if (a1[i] <= a2[j]) {
		a[k] = a1[i];
		b[k] = b1[i];
		i++;
	    } else {
		a[k] = a2[j];
		b[k] = b2[j];
		j++;
	    }
	}
    }

}

