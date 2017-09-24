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
import java.util.Random;
import java.util.Collections;

/**
 * A Gene object represents a gene (either protein-coding or RNA).
 * A Gene object consists of a variety of data about a gene,
 * including its coordinates, strand, name, product, and 
 * expression information in each replicate (raw counts mapping
 * to gene and normalized counts mapping to gene) or in each
 * condition (mean, variance, lowess, RPKM) or in a pair of
 * conditions (p-value of differential expression).
 */
public class Gene {

    /*****************************************
     **********   CLASS VARIABLES   **********
     *****************************************/

    private static Random rand = new Random();  // Random number generator



    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private int start;  // Translation start
    private int stop;   // Translation stop
    private char strand;
    private String ID;
    private String name;
    private String synonym;
    private String product;
    private String type; // ORF or RNA
    private int tStart;  // Transcription start
    private int tStop;   // Transcription stop

    private ArrayList<ArrayList<Long>> rawCounts;  // Reads mapping to each NT of gene
    private ArrayList<ArrayList<Long>> rawCounts_reads;  // Reads mapping to gene
    private ArrayList<ArrayList<Long>> normalizedCounts;
    private ArrayList<Long> RPKMs;
    private ArrayList<Long> means;
    private long[][] variances;  // 2D array
    private long[][] lowess;     // 2D array
    private ArrayList<Double> pValues;  // Differential expression
    private ArrayList<Double> qValues;  // Differential expression (corrected)



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    /**
     * Constructs a new Genome object based on a line from a
     * gene file (either *.ptt or *.rnt).
     */
    public Gene(String line, String type) {
	String[] parse_line = line.split("\t");
	if (parse_line.length < 9) {
	    Rockhopper.output("Error - expecting 9 columns of gene information but found less than 9:\t" + line + "\n");
	    return;
	}
	String[] parse_coords = parse_line[0].split("\\.");
	this.type = type;
	this.strand = parse_line[1].charAt(0);
	this.ID = parse_line[3];
	this.name = parse_line[4];
	this.synonym = parse_line[5];
	this.product = parse_line[8];

	// Set coordinates based on strand and type
	int x = Integer.parseInt(parse_coords[0]);  // First coord
	int y = Integer.parseInt(parse_coords[2]);  // Second coord
	if (type.equals("ORF") && (strand == '+')) {  // Plus strand ORF
	    start = x;
	    stop = y;
	    tStart = 0;
	    tStop = 0;
	} else if (type.equals("ORF") && (strand == '-')) {  // Minus strand ORF
	    start = y;
	    stop = x;
	    tStart = 0;
	    tStop = 0;
	} else if (type.equals("RNA") && (strand == '+')) {  // Plus strand RNA
	    start = 0;
	    stop = 0;
	    tStart = x;
	    tStop = y;
	} else if (type.equals("RNA") && (strand == '-')) {  // Minus strand RNA
	    start = 0;
	    stop = 0;
	    tStart = y;
	    tStop = x;
	} else if (type.equals("RNA") && (strand == '?')) {  // Ambiguous strand RNA
	    start = 0;
	    stop = 0;
	    tStart = x;
	    tStop = y;
	} else
	    Rockhopper.output("Error - this case should be unreachable!\n");
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    /**
     * Return the gene's strand.
     */
    public char getStrand() {
	return strand;
    }

    /**
     * Return the gene's start coordinate of translation.
     */
    public int getStart() {
	return start;
    }

    /**
     * Return the gene's stop coordinate of translation.
     */
    public int getStop() {
	return stop;
    }

    /**
     * Reeturn the gene's name.
     */
    public String getName() {
	if (name.length() > 1) return name;
	return synonym;
    }

    /**
     * Returns true if this gene is a protein coding gene,
     * false otherwise.
     */
    public boolean isORF() {
	return type.equals("ORF");
    }

    /**
     * Returns this gene's product.
     */
    public String getProduct() {
	return product;
    }

    /**
     * Set the transcription start coordinate.
     */
    public void setStartT(int tStart) {
	this.tStart = tStart;
    }

    /**
     * Set the transcription stop coordinate.
     */
    public void setStopT(int tStop) {
	this.tStop = tStop;
    }

    /**
     * Return the gene's start coordinate of transcription.
     */
    public int getStartT() {
	return this.tStart;
    }

    /**
     * Return the gene's stop coordinate of transcription.
     */
    public int getStopT() {
	return this.tStop;
    }

    /**
     * Returns the first (smallest) coordinate of this Gene.
     * If this Gene is an ORF it returns the smallest translation coordinate.
     * If this Gene is an RNA it returns the smallest transcription coordinate.
     */
    public int getFirst() {
	if (isORF()) return Math.min(start, stop);
	return Math.min(tStart, tStop);
    }

    /**
     * Returns the last (largest) coordinate of this Gene.
     * If this Gene is an ORF it returns the largest translation coordinate.
     * If this Gene is an RNA it returns the largest transcription coordinate.
     */
    public int getLast() {
	if (isORF()) return Math.max(start, stop);
	return Math.max(tStart, tStop);
    }

    /**
     * Return the number of reads mapping to the gene in the
     * specified condition and specified replicate.
     */
    public long getRawCount(int condition, int replicate) {
	if (condition < rawCounts.size()) {
	    if (replicate < rawCounts.get(condition).size())
		return rawCounts.get(condition).get(replicate);
	}
	return 0;
    }

    /**
     * Return the number of reads mapping to the gene in the
     * specified condition and specified replicate.
     */
    public long getRawCount_reads(int condition, int replicate) {
	if (condition < rawCounts_reads.size()) {
	    if (replicate < rawCounts_reads.get(condition).size())
		return rawCounts_reads.get(condition).get(replicate);
	}
	return 0;
    }

    /**
     * Return the normalized number of reads mapping to the gene in the
     * specified condition and specified replicate.
     */
    public long getNormalizedCount(int condition, int replicate) {
	if (condition < normalizedCounts.size()) {
	    if (replicate < normalizedCounts.get(condition).size())
		return normalizedCounts.get(condition).get(replicate);
	}
	return 0;
    }

    /**
     * Return the mean reads mapping to the gene in the
     * specified condition.
     */
    public long getMean(int condition) {
	if (condition < means.size())
	    return means.get(condition);
	return 0;
    }

    /**
     * Returns the average expression of the gene (averaged
     * over the length of the gene) in the specified
     * condition.
     */
    public long getAvg(int condition) {
	long avg = 0;
	if (type.compareTo("ORF") == 0)
	    avg = means.get(condition)/(Math.max(this.start,this.stop)-Math.min(this.start,this.stop)+1);
	if (type.compareTo("RNA") == 0)
	    avg = means.get(condition)/(Math.max(this.tStart,this.tStop)-Math.min(this.tStart,this.tStop)+1);
	return avg;
    }

    /**
     * Return the number of replicates in the specified condition.
     */
    public int getNumReplicates(int condition) {
	if (condition >= rawCounts.size()) return 0;
	else return rawCounts.get(condition).size();
    }

    /**
     * Returns the minimum q-value for this Gene.
     */
    public double getMinQvalue() {
	double min = 1.0;
	for (int i=0; i<qValues.size(); i++) min = Math.min(min, qValues.get(i));
	return min;
    }

    public boolean hasQvalue(int c) {
	return ((qValues != null) && (qValues.size() > c));
    }

    /**
     * Returns a String representation of a a Gene's expression 
     * in each condition and p-values of differential expression.
     */
    public String expressionToString() {
	StringBuilder sb = new StringBuilder();
	for (int i=0; i<means.size(); i++) {
	    if (Rockhopper.verbose) {  // Verbose output
		for (int j=0; j<rawCounts.get(i).size(); j++) 
		    sb.append("\t" + rawCounts_reads.get(i).get(j));
		for (int j=0; j<normalizedCounts.get(i).size(); j++) 
		    sb.append("\t" + normalizedCounts.get(i).get(j));
		sb.append("\t" + RPKMs.get(i));
	    }
	    sb.append("\t" + getAvg(i));
	}
	for (int i=0; i<qValues.size(); i++) {
	    if (Rockhopper.verbose) sb.append("\t" + pValues.get(i).toString());
	    sb.append("\t" + qValues.get(i).toString());
	}
	return sb.toString();
    }

    /**
     * Returns a String representation of a Gene.
     */
    public String toString() {
	StringBuilder sb = new StringBuilder();
	String tStart = "";
	String start = "";
	String stop = "";
	String tStop = "";
	if (this.tStart > 0) tStart += this.tStart;
	if (this.start > 0) start += this.start;
	if (this.stop > 0) stop += this.stop;
	if (this.tStop > 0) tStop += this.tStop;
	sb.append(tStart + "\t" + start + "\t" + stop + "\t" + tStop + "\t" + strand + "\t" + name + "\t" + synonym + "\t" + product);
	return sb.toString();
    }

    /**
     * Returns the coordinate (among transcription start/stop coordinates
     * and translation start/stop coordinates) with minimum value.
     */
    public int getMinCoordinate() {
	int minCoord = Integer.MAX_VALUE;
	if (tStart > 0) minCoord = Math.min(minCoord, tStart);
	if (start > 0) minCoord = Math.min(minCoord, start);
	if (stop > 0) minCoord = Math.min(minCoord, stop);
	if (tStop > 0) minCoord = Math.min(minCoord, tStop);
	return minCoord;
    }

    /**
     * Returns the coordinate (among transcription start/stop coordinates
     * and translation start/stop coordinates) with maximum value.
     */
    public int getMaxCoordinate() {
	int maxCoord = -1;
	if (tStart > 0) maxCoord = Math.max(maxCoord, tStart);
	if (start > 0) maxCoord = Math.max(maxCoord, start);
	if (stop > 0) maxCoord = Math.max(maxCoord, stop);
	if (tStop > 0) maxCoord = Math.max(maxCoord, tStop);
	return maxCoord;
    }

    /**
     * Set the number of reads mapping to this Gene in the specified
     * replicate in the specified condition.
     */
    public void setRawCount(int condition, int replicate, long readsForGene) {
	// Initialize instance variables
	if (rawCounts == null) {
	    rawCounts = new ArrayList<ArrayList<Long>>();
	    rawCounts_reads = new ArrayList<ArrayList<Long>>();
	    normalizedCounts = new ArrayList<ArrayList<Long>>();
	    RPKMs = new ArrayList<Long>();
	    means = new ArrayList<Long>();
	}
	while (rawCounts.size() < condition+1) {
	    rawCounts.add(new ArrayList<Long>());
	    rawCounts_reads.add(new ArrayList<Long>());
	    normalizedCounts.add(new ArrayList<Long>());
	    RPKMs.add((long)0);
	    means.add((long)0);
	}
	while (rawCounts.get(condition).size() < replicate+1) {
	    rawCounts.get(condition).add((long)0);
	    rawCounts_reads.get(condition).add((long)0);
	    normalizedCounts.get(condition).add((long)0);
	}

	// Populate rawCounts
	rawCounts.get(condition).set(replicate, readsForGene);
    }

    /**
     * Set the number of reads mapping to this Gene in the specified
     * replicate in the specified condition.
     */
    public void setRawCount_reads(int condition, int replicate, long readsForGene) {
	if (condition < normalizedCounts.size()) {
	    if (replicate < normalizedCounts.get(condition).size()) {
		rawCounts_reads.get(condition).set(replicate, readsForGene);
	    }
	}
    }

    /**
     * Set the normalized number of reads mapping to this Gene in the specified
     * replicate in the specified condition.
     */
    public void setNormalizedCount(int condition, int replicate, double scalingFactor, long upperQuartile) {
	if (condition < normalizedCounts.size()) {
	    if (replicate < normalizedCounts.get(condition).size()) {
		double multiplier = scalingFactor / (double)upperQuartile;
		long normalizedCount = (long)(multiplier * getRawCount(condition, replicate));
		normalizedCounts.get(condition).set(replicate, normalizedCount);
	    }
	}
    }

    /**
     * For each condition, across all replicates of the conditions, compute
     * the mean and RPKM for this Gene.
     */
    public void computeExpression(ArrayList<Condition> conditions) {
	for (int i=0; i<means.size(); i++) {  // For each condition

	    // Compute RPKM and mean
	    long sumRawCounts = 0;
	    long totalCounts = 0;
	    long mean = 0;
	    for (int j=0; j<rawCounts.get(i).size(); j++) {
		sumRawCounts += rawCounts.get(i).get(j);
		totalCounts += conditions.get(i).getReplicate(j).getTotalReads();
		mean += normalizedCounts.get(i).get(j);
	    }
	    sumRawCounts /= rawCounts.get(i).size();
	    totalCounts /= rawCounts.get(i).size();
	    mean /= normalizedCounts.get(i).size();
	    if (totalCounts == 0) RPKMs.set(i, (long)0);  // Avoid divide-by-zero error
	    else RPKMs.set(i, (long)(1000000000 * sumRawCounts / (totalCounts * (getMaxCoordinate()-getMinCoordinate()+1))));
	    means.set(i, mean);
	}
    }

    /**
     * For each condition, across all replicates of the conditions, compute
     * the variance for this Gene.
     */
    public void computeVariance(ArrayList<Condition> conditions) {
	double varianceAdjustmentNoReplicates = 1.10;
	double varianceAdjustmentReplicates = 1.20;
	variances = new long[conditions.size()][conditions.size()];
	lowess = new long[conditions.size()][conditions.size()];

	for (int x=0; x<conditions.size(); x++) {  // First of pair
	    for (int y=0; y<conditions.size(); y++) {  // Second of pair

		if (x == y) continue;  // No need to compute variance with self

		if (conditions.get(x).numReplicates() == 1) {  // No replicates. Use partner surrogate.
		    int partner = y;
		    long mean = (means.get(x) + means.get(partner)) / 2;
		    long variance = ((means.get(x)-mean)*(means.get(x)-mean) + (means.get(partner)-mean)*(means.get(partner)-mean)) / 1;  // Sample standard deviation
		    variances[x][y] = (long)Math.pow(variance, varianceAdjustmentNoReplicates);
		} else {  // We have replicates. Use the replicates to compute the variance.
		    long variance = 0;
		    for (int j=0; j<conditions.get(x).numReplicates(); j++) {
			variance += (long)(normalizedCounts.get(x).get(j) - means.get(x)) * (long)(normalizedCounts.get(x).get(j) - means.get(x));
		    }
		    variance /= (conditions.get(x).numReplicates() - 1);  // Sample standard deviation
		    variances[x][y] = (long)Math.pow(variance, varianceAdjustmentReplicates);
		}
	    }
	}
    }

    /**
     * For each pair of conditions, compute the p-value of differntial
     * expression for this Gene.
     */
    public void computeDifferentialExpression() {
	pValues = new ArrayList<Double>();
	qValues = new ArrayList<Double>();
	for (int x=0; x<normalizedCounts.size()-1; x++) {  // First of two conditions
	    for (int y=x+1; y<normalizedCounts.size(); y++) {  // Second of two conditions

		double k_A = 0.0;
		double k_B = 0.0;
		for (int j=0; j<normalizedCounts.get(x).size(); j++) k_A += normalizedCounts.get(x).get(j);
		for (int j=0; j<normalizedCounts.get(y).size(); j++) k_B += normalizedCounts.get(y).get(j);

		if (normalizedCounts.get(x).size() < normalizedCounts.get(y).size())
		    k_B *= normalizedCounts.get(x).size() / (double)(normalizedCounts.get(y).size());
		else if (normalizedCounts.get(x).size() > normalizedCounts.get(y).size())
		    k_A *= normalizedCounts.get(y).size() / (double)(normalizedCounts.get(x).size());

		double q = k_A + k_B;

		double mean_A = q;
		double mean_B = q;
		double variance_A = lowess[x][y];
		double variance_B = lowess[y][x];

		double p_a = mean_A / variance_A;
		double p_b = mean_B / variance_B;
		double r_a = Math.max(mean_A*mean_A / (variance_A - mean_A), 1.0);  // r should not be < 1
		double r_b = Math.max(mean_B*mean_B / (variance_B - mean_B), 1.0);  // r should not be < 1

		if ((p_a < 0.0) || (p_b < 0.0) || (p_a > 1.0) || (p_b > 1.0) || (variance_A == 0.0) || (variance_B == 0.0)) {
		    pValues.add(1.0);
		    qValues.add(1.0);
		    continue;
		}

		// Compute p-value of differential expression in two conditions
		double p_ab = NegativeBinomial.pmf(r_a-1, k_A+r_a-1, p_a) * NegativeBinomial.pmf(r_b-1, k_B+r_b-1, p_b);
		long k_sum = (long)(k_A + k_B);

		// Fast p-value estimation
		double numerator = 0.0;
		double denominator = 0.0;
		long mode = (long)k_B;

		long a = mode;  // Begin near middle
		long increment = 1;
		long alpha = 1000;  // Number of times we increment by 1 (raising alpha raises precision but slows down computation)
		double previous_p = 0.0;
		while (a <= k_sum) {
		    long b = k_sum - a;
		    double current_p = NegativeBinomial.pmf(r_a-1, a+r_a-1, p_a) * NegativeBinomial.pmf(r_b-1, b+r_b-1, p_b);
		    denominator += current_p;
		    if (current_p <= p_ab) numerator += current_p;
		    if (increment > 1) {
			double average_p = (current_p+previous_p)/2.0;
			denominator += average_p * (increment-1);
			if (average_p <= p_ab) numerator += average_p * (increment-1);
		    }
		    previous_p = current_p;
		    if (a - mode >= alpha) {
			alpha *= 2;
			increment *= 2;
		    }
		    a += increment;
		}
		a = mode;  // Begin near middle
		long decrement = 1;
		alpha = 1000;  // Number of times we increment by 1 (raising alpha raises precision but slows down computation)
		previous_p = 0.0;
		while (a >= 0) {
		    long b = k_sum - a;
		    double current_p = NegativeBinomial.pmf(r_a-1, a+r_a-1, p_a) * NegativeBinomial.pmf(r_b-1, b+r_b-1, p_b);
		    denominator += current_p;
		    if (current_p <= p_ab) numerator += current_p;
		    if (decrement > 1) {
			double average_p = (previous_p + current_p) / 2.0;
			denominator += average_p * (decrement-1);
			if (average_p <= p_ab) numerator += average_p * (decrement-1);
		    }
		    previous_p = current_p;
		    if (mode - a >= alpha) {
			alpha *= 2;
			decrement *= 2;
		    }
		    a -= decrement;
		}
		double p_value = 1.0;
		if (denominator != 0.0) p_value = numerator / denominator;
		pValues.add(p_value);
		qValues.add(1.0);
	    }
	}
    }

    /**
     * Returns true is this Gene is an ORF and if it is differentially 
     * expressed in at least one pair of conditions at the specified
     * significance level. Returns false otherwise.
     */
    public boolean isDifferntiallyExpressedORF(double significance) {
	if (!isORF()) return false;
	for (int i=0; i<qValues.size(); i++) {
	    if (qValues.get(i) <= significance) return true;
	}
	return false;
    }



    /**********************************************
     **********   PUBLIC CLASS METHODS   **********
     **********************************************/

    /**
     * Computes the Lowess variance for each Gene in the Genome.
     */
    public static void setLowessVariances(ArrayList<Genome> genomes, ArrayList<Condition> conditions) {
	for (int x=0; x<conditions.size(); x++) {  // First of pair
	    for (int y=0; y<conditions.size(); y++) {  // Second of pair

		if (x == y) continue;  // No need to compute lowess for self
		/*
		if ((x > 0) && (conditions.get(x).numReplicates() > 1)) {  // We have replicates
		    for (int z=0; z<genomes.size(); z++) {
			Genome genome = genomes.get(z);
			for (int j=0; j<genome.numGenes(); j++)
			    genome.getGene(j).lowess[x][y] = genome.getGene(j).lowess[0][y];
		    }		    
		    continue;
		}
		*/
		/* (x,y) lowess is NOT the same as (y,x) lowess
		if ((x > y) && (conditions.get(x).numReplicates() == 1) && (conditions.get(y).numReplicates() == 1)) {  // Already computed lowess[y][x], which is same as lowess[x][y]
		    for (int z=0; z<genomes.size(); z++) {
			Genome genome = genomes.get(z);
			for (int j=0; j<genome.numGenes(); j++)
			    genome.getGene(j).lowess[x][y] = genome.getGene(j).lowess[y][x];
		    }
		    continue;
		}
		*/

		double b = 0.0;  // Bias correction term
		for (int j=0; j<conditions.get(x).numReplicates(); j++) {
		    b += 100000.0 / (double)conditions.get(x).getReplicate(j).getUpperQuartile();
		}
		b /= conditions.get(x).numReplicates();

		// Create list of gene expressions and list of gene variances
		ArrayList<Long> expression = new ArrayList<Long>();
		ArrayList<Long> variance = new ArrayList<Long>();
		for (int z=0; z<genomes.size(); z++) {
		    Genome genome = genomes.get(z);
		    for (int j=0; j<genome.numGenes(); j++) {
			expression.add(genome.getGene(j).means.get(x));
			variance.add(genome.getGene(j).variances[x][y]);
		    }
		}

		// Perform Lowess computation
		ArrayList<Long> lowessVariance = Lowess.lowess(expression, variance);

		// Uncomment to output data for Lowess graph to StdOut
		//System.out.println("Mean" + "\t" + "Variance" + "\t" + "Lowess");
		//for (int k=0; k<lowessVariance.size(); k++) System.out.println(expression.get(k) + "\t" + variance.get(k) + "\t" + ((long)(lowessVariance.get(k) - expression.get(k)*b)));

		// Assign each gene its lowess variances (after subtracting bias correction term)
		int previousGenomeSizes = 0;
		for (int z=0; z<genomes.size(); z++) {
		    Genome genome = genomes.get(z);
		    for (int j=0; j<genome.numGenes(); j++)
			genome.getGene(j).lowess[x][y] = (long)(lowessVariance.get(previousGenomeSizes + j) - (genome.getGene(j).getMean(x) * b));
		    previousGenomeSizes += genome.numGenes();
		}
	    }
	}
    }

    /**
     * Computes q-values for each gene, i.e., corrected p-values,
     * using Benjamini Hochberg correction.
     */
    public static void correctPvalues(ArrayList<Genome> genomes, ArrayList<Condition> conditions) {
	int totalGenes = 0;
	for (int z=0; z<genomes.size(); z++) totalGenes += genomes.get(z).numGenes();

	int pValue_index = 0;
	for (int x=0; x<genomes.get(0).getGene(0).means.size()-1; x++) {  // First in pair
	    for (int y=x+1; y<genomes.get(0).getGene(0).means.size(); y++) {  // Second in pair

		int[] indices = new int[totalGenes];
		int[] genomeIndices = new int[totalGenes];
		double[] pvalues = new double[totalGenes];
		int previousGenomeSizes = 0;
		for (int z=0; z<genomes.size(); z++) {
		    Genome genome = genomes.get(z);
		    for (int j=0; j<genome.numGenes(); j++) {
			pvalues[previousGenomeSizes + j] = genome.getGene(j).pValues.get(pValue_index);
			indices[previousGenomeSizes + j] = j;
			genomeIndices[previousGenomeSizes + j] = z;

			// Check if there is too little expression to compute a p-value
			Gene g = genome.getGene(j);
			double e1 = g.means.get(x);
			double e2 = g.means.get(y);
			if (g.isORF()) {  // ORF
			    e1 /= Math.max(g.start, g.stop) - Math.min(g.start, g.stop) + 1;
			    e2 /= Math.max(g.start, g.stop) - Math.min(g.start, g.stop) + 1;
			} else {  // RNA
			    e1 /= Math.max(g.tStart, g.tStop) - Math.min(g.tStart, g.tStop) + 1;
			    e2 /= Math.max(g.tStart, g.tStop) - Math.min(g.tStart, g.tStop) + 1;
			}
			if ((e1 < conditions.get(x).getMinDiffExpressionLevel()) &&
			    (e2 < conditions.get(y).getMinDiffExpressionLevel()))
			    pvalues[previousGenomeSizes + j] = 1.0;
		    }
		    previousGenomeSizes += genome.numGenes();
		}
		mergesort(pvalues, indices, genomeIndices, 0, totalGenes-1);
		double previous_BH_value = 0.0;
		for (int k=0; k<pvalues.length; k++) {
		    double BH_value = pvalues[k] * totalGenes / (k+1);
		    BH_value = Math.min(BH_value, 1.0);  // Disallow values greater than 1

		    // To preserve monotonicity, we take maximum to ensure we don't
		    // generate a value less than the previous value.
		    BH_value = Math.max(BH_value, previous_BH_value);
		    previous_BH_value = BH_value;
		    genomes.get(genomeIndices[k]).getGene(indices[k]).qValues.set(pValue_index, BH_value);
		}
		pValue_index++;
	    }
	}
    }



    /***********************************************
     **********   PRIVATE CLASS METHODS   **********
     ***********************************************/

    /**
     * Mergesort parallel arrays "a" and "b" and "c" based on values in "a"
     */
    private static void mergesort(double[] a, int[] b, int[] c, int lo, int hi) {
	if (lo < hi) {
	    int q = (lo+hi)/2;
	    mergesort(a, b, c, lo, q);
	    mergesort(a, b, c, q+1, hi);
	    merge(a, b, c, lo, q, hi);
	}
    }

    /**
     * Mergesort helper method
     */
    private static void merge(double[] a, int[] b, int[] c, int lo, int q, int hi) {
	double[] a1 = new double[q-lo+1];
	double[] a2 = new double[hi-q];
	int[] b1 = new int[q-lo+1];
	int[] b2 = new int[hi-q];
	int[] c1 = new int[q-lo+1];
	int[] c2 = new int[hi-q];
	for (int i=0; i<a1.length; i++) {
	    a1[i] = a[lo+i];
	    b1[i] = b[lo+i];
	    c1[i] = c[lo+i];
	}
	for (int j=0; j<a2.length; j++) {
	    a2[j] = a[q+1+j];
	    b2[j] = b[q+1+j];
	    c2[j] = c[q+1+j];
	}
	int i=0;  // Index for first half arrays
	int j=0;  // Index for second half arrays
	for (int k=lo; k<=hi; k++) {
	    if (i >= a1.length) {
		a[k] = a2[j];
		b[k] = b2[j];
		c[k] = c2[j];
		j++;
	    } else if (j >= a2.length) {
		a[k] = a1[i];
		b[k] = b1[i];
		c[k] = c1[i];
		i++;
	    } else if (a1[i] <= a2[j]) {
		a[k] = a1[i];
		b[k] = b1[i];
		c[k] = c1[i];
		i++;
	    } else {
		a[k] = a2[j];
		b[k] = b2[j];
		c[k] = c2[j];
		j++;
	    }
	}
    }

}

