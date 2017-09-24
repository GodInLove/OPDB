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

/**
 * An instance of the DeNovoTranscript class represents
 * a de novo assembled transcript.
 */
public class DeNovoTranscript {

    /********************************************
     **********   Instance Variables   **********
     ********************************************/

    private String sequence;
    private ArrayList<ArrayList<Long>> rawCounts_nts;
    private ArrayList<ArrayList<Long>> rawCounts_reads;
    private ArrayList<ArrayList<Long>> normalizedCounts;
    private ArrayList<Long> RPKMs;
    private ArrayList<Long> means;
    public long[][] variances;  // 2D array
    public long[][] lowess;     // 2D array
    private ArrayList<Double> pValues;  // Differential expression
    private ArrayList<Double> qValues;  // Differential expression (corrected)



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    public DeNovoTranscript(String sequence, long[] readCounts_nts_linear) {
	this.sequence = sequence;

	// Convert linear array to 2D ArrayList based on Conditions/Replicates.
	rawCounts_nts = new ArrayList<ArrayList<Long>>(Assembler.conditionFiles.size());
	rawCounts_reads = new ArrayList<ArrayList<Long>>(Assembler.conditionFiles.size());
	normalizedCounts = new ArrayList<ArrayList<Long>>(Assembler.conditionFiles.size());
	RPKMs = new ArrayList<Long>(Assembler.conditionFiles.size());
	means = new ArrayList<Long>(Assembler.conditionFiles.size());
	int linearIndex = 0;
	for (int i=0; i<Assembler.conditionFiles.size(); i++) {
	    String[] files = Assembler.conditionFiles.get(i).split(",");
	    rawCounts_nts.add(new ArrayList<Long>(files.length));
	    rawCounts_reads.add(new ArrayList<Long>(files.length));
	    normalizedCounts.add(new ArrayList<Long>(files.length));
	    RPKMs.add((long)0);
	    means.add((long)0);
	    for (int j=0; j<files.length; j++) {
		rawCounts_nts.get(i).add(readCounts_nts_linear[linearIndex]);
		rawCounts_reads.get(i).add(readCounts_nts_linear[linearIndex] / DeNovoTranscripts.avgLengthOfReads[linearIndex]);
		normalizedCounts.get(i).add((long)0);
		linearIndex++;
	    }
	}
    }

    /**
     * Used when reading in from compressed file.
     */
    public DeNovoTranscript(String sequence, ArrayList<ArrayList<Long>> rawCounts_nts, ArrayList<ArrayList<Long>> rawCounts_reads, ArrayList<ArrayList<Long>> normalizedCounts, ArrayList<Long> RPKMs, ArrayList<Long> means, long[][] variances, long[][] lowess, ArrayList<Double> pValues, ArrayList<Double> qValues) {
	this.sequence = sequence;
	this.rawCounts_nts = rawCounts_nts;
	this.rawCounts_reads = rawCounts_reads;
	this.normalizedCounts = normalizedCounts;
	this.RPKMs = RPKMs;
	this.means = means;
	this.variances = variances;
	this.lowess = lowess;
	this.pValues = pValues;
	this.qValues = qValues;
    }



    /*************************************************
     **********   Public Instance Methods   **********
     *************************************************/

    /**
     * Returns the length of the transcript sequence.
     */
    public int length() {
	return sequence.length();
    }

    /**
     * Returns the sequence of this transcript.
     */
    public String sequence() {
	return this.sequence;
    }

    /**
     * Returns the raw count of nts reads mapping to 
     * this transcript for the specified condition and
     * replicate.
     */
    public long getRawCountNTs(int condition, int replicate) {
	return rawCounts_nts.get(condition).get(replicate);
    }

    /**
     * Returns the raw count of full length reads mapping
     * to this transcript for the specified condition and
     * replicate.
     */
    public long getRawCountReads(int condition, int replicate) {
	return rawCounts_reads.get(condition).get(replicate);
    }

    /**
     * Returns the mean expression of this transcript in
     * the specified condition.
     */
    public long getMean(int condition) {
	return means.get(condition);
    }

    /**
     * Returns the p-value of this transcript for
     * the specified pair of conditions.
     */
    public double getPvalue(int conditionPair) {
	return pValues.get(conditionPair);
    }

    /**
     * Returns the number of p-values, i.e.,
     * the number of pairs of conditions, 
     * for this transcript.
     */
    public int getNumPvalues() {
	return pValues.size();
    }

    /**
     * Set the normalized number of reads mapping to this
     * transcript in the specified condition and replicate.
     */
    public void setNormalizedCount(int condition, int replicate, double scalingFactor, long upperQuartile) {
	double multiplier = scalingFactor / (double)upperQuartile;
	long normalizedCount = (long)(multiplier * getRawCountNTs(condition, replicate));
	normalizedCounts.get(condition).set(replicate, normalizedCount);
    }

    /**
     * Sets the q-value for the specified pair of conditions.
     */
    public void setQvalue(int conditionPair, double value) {
	qValues.set(conditionPair, value);
    }

    /**
     * Returns an array representation of rawCounts_nts
     * in the specified condition.
     */
    public long[] getRawCounts_nts_array(int condition) {
	return listToArray_long(rawCounts_nts.get(condition));
    }

    /**
     * Returns an array representation of rawCounts_reads
     * in the specified condition.
     */
    public long[] getRawCounts_reads_array(int condition) {
	return listToArray_long(rawCounts_reads.get(condition));
    }

    /**
     * Returns an array representation of normalizedCounts
     * in the specified condition.
     */
    public long[] getNormalizedCounts_array(int condition) {
	return listToArray_long(normalizedCounts.get(condition));
    }

    /**
     * Returns an array representation of RPKMs.
     */
    public long[] getRPKMs() {
	return listToArray_long(RPKMs);
    }

    /**
     * Returns an array representation of means.
     */
    public long[] getMeans() {
	return listToArray_long(means);
    }

    /**
     * Returns a 1D array representation of variances.
     */
    public long[] getVariances() {
	return array2D_to_1D_long(variances);
    }

    /**
     * Returns a 1D array representation of lowess.
     */
    public long[] getLowess() {
	return array2D_to_1D_long(lowess);
    }

    /**
     * Returns an array representation of pValues.
     */
    public double[] getPvalues() {
	return listToArray_double(pValues);
    }

    /**
     * Returns an array representation of qValues.
     */
    public double[] getQvalues() {
	return listToArray_double(qValues);
    }

    /**
     * Calculates the RPKM and mean expression for
     * this transcript in the specified condition.
     */
    public void computeExpression(int condition) {
	long sumRawCounts = 0;
	long totalCounts = 0;
	long mean = 0;
	for (int j=0; j<rawCounts_nts.get(condition).size(); j++) {
	    sumRawCounts += rawCounts_nts.get(condition).get(j);
	    totalCounts += DeNovoTranscripts.totalReads.get(condition).get(j);
	    mean += normalizedCounts.get(condition).get(j);
	}
	sumRawCounts /= rawCounts_nts.get(condition).size();
	totalCounts /= rawCounts_nts.get(condition).size();
	mean /= normalizedCounts.get(condition).size();
	if (totalCounts == 0) RPKMs.set(condition, (long)0);  // Avoid divide-by-zero error
	else RPKMs.set(condition, (long)(1000000000 * sumRawCounts / (totalCounts * sequence.length())));
	means.set(condition, mean);
    }

    /**
     * Calculates the variance for
     * this transcript in the specified condition.
     */
    public void computeVariance(int numConditions) {
	double varianceAdjustmentNoReplicates = 1.10;
	double varianceAdjustmentReplicates = 1.20;
	variances = new long[numConditions][numConditions];
	lowess = new long[numConditions][numConditions];

	for (int x=0; x<numConditions; x++) {  // First of pair
	    for (int y=0; y<numConditions; y++) {  // Second of pair

		if (x == y) continue;  // No need to compute variance with self

		if (normalizedCounts.get(x).size() == 1) {  // No replicates. Use partner surrogate.
		    int partner = y;
		    long mean = (means.get(x) + means.get(partner)) / 2;
		    long variance = ((means.get(x)-mean)*(means.get(x)-mean) + (means.get(partner)-mean)*(means.get(partner)-mean)) / 1;  // Sample standard deviation
		    //variances[x][y] = variance;
		    variances[x][y] = (long)Math.pow(variance, varianceAdjustmentNoReplicates);
		} else {  // We have replicates. Use the replicates to compute the variance.
		    long variance = 0;
		    for (int j=0; j<normalizedCounts.get(x).size(); j++) {
			variance += (long)(normalizedCounts.get(x).get(j) - means.get(x)) * (long)(normalizedCounts.get(x).get(j) - means.get(x));
		    }
		    variance /= (normalizedCounts.get(x).size() - 1);  // Sample standard deviation
		    //variances[x][y] = variance;
		    variances[x][y] = (long)Math.pow(variance, varianceAdjustmentReplicates);
		}
	    }
	}
    }

    /**
     * For each pair of conditions, compute the p-value of differntial
     * expression for this transcript.
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
     * Returns a String representation of this transcript.
     */
    public String toString() {
	StringBuilder sb = new StringBuilder();
	sb.append(sequence + "\t" + sequence.length());
	for (int i=0; i<normalizedCounts.size(); i++) {
	    if (Assembler.verbose) {  // Verbose output
		for (int j=0; j<normalizedCounts.get(i).size(); j++) sb.append("\t" + rawCounts_reads.get(i).get(j));
		for (int j=0; j<normalizedCounts.get(i).size(); j++) sb.append("\t" + normalizedCounts.get(i).get(j));
		sb.append("\t" + RPKMs.get(i));
	    }
	    sb.append("\t" + (means.get(i)/sequence.length()));
	}
	for (int i=0; i<qValues.size(); i++) {
	    if (Assembler.verbose) sb.append("\t" + pValues.get(i).toString());
	    sb.append("\t" + qValues.get(i).toString());
	}
	return sb.toString(); 
    }



    /**************************************************
     **********   Private Instance Methods   **********
     **************************************************/

    /**
     * Returns an array representation of a list of longs.
     */
    private long[] listToArray_long(ArrayList<Long> x) {
	long[] a = new long[x.size()];
	for (int i=0; i<x.size(); i++) a[i] = x.get(i);
	return a;
    }

    /**
     * Returns an array representation of a list of doubles.
     */
    private double[] listToArray_double(ArrayList<Double> x) {
	double[] a = new double[x.size()];
	for (int i=0; i<x.size(); i++) a[i] = x.get(i);
	return a;
    }

    /**
     * Returns a 1D array representation of a 2D array of longs.
     */
    private long[] array2D_to_1D_long(long[][] x) {
	long[] a = new long[x.length*x.length];
	int index = 0;
	for (int i=0; i<x.length; i++) {
	    for (int j=0; j<x.length; j++) {
		a[index] = x[i][j];
		index++;
	    }
	}
	return a;
    }

}

