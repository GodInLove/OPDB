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
import java.util.Scanner;
import java.util.Collections;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

public class Operons {

    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private SmoothDistribution lengthSame;
    private SmoothDistribution lengthOpp;
    private SmoothDistribution corrSame;
    private SmoothDistribution corrOpp;
    private double operonPrior;
    private double nonOperonPrior;
    private double[] operonIGexpression;     // Distribution of IG expression for operons
    private double[] nonOperonIGexpression;  // Distribution of IG expression for non-operons
    private ArrayList<Double> p_values;



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    public Operons(String geneFileName) {
	this(readInGenes(geneFileName), readInGenes(geneFileName));
    }

    public Operons(ArrayList<Gene> codingGenes, ArrayList<Gene> genes) {
	setOperonPrior(codingGenes);
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    /**
     * Returns the number of gene pairs predicted to be co-transcribed
     * as part of an operon.
     */
    public int getNumOperonGenePairs(ArrayList<Gene> genes) {
	int numOperonGenePairs = 0;
	for (int i=1; i<genes.size(); i++) {
	    if (isGenePairAnOperon(genes, i-1, i)) numOperonGenePairs++;
	}
	return numOperonGenePairs;
    }

    /**
     * Returns true if the two genes at the specified indices, which must be consecutive, 
     * are predicted to be co-transcribed. Returns false otherwise.
     */
    public boolean isGenePairAnOperon(ArrayList<Gene> genes, int x, int y) {
	if (x != y-1) {
	    Rockhopper.output("Error - two indices must be consecutive.\n");
	    return false;
	}
	if (genes.get(x).getStrand() != genes.get(y).getStrand()) return false;  // Genes on different strands
	int length = getIGlength(genes.get(x), genes.get(y));
	if ((length < 40) || ((length < 100) && (getCorrelation(genes.get(x), genes.get(y)) >= 0.5))) return true;
	return false;
    }

    /**
     * Output to file gene-pairs predicted to be part of the same operon.
     */
    public void outputGenePairOperons(String operonOutputFile, ArrayList<Gene> genes) {
	try {
	    PrintWriter writer = new PrintWriter(new File(operonOutputFile));
	    writer.println("Transcription Start" + "\t" + "Translation Start" + "\t" + "Translation Stop" + "\t" + "Transcription Stop" + "\t" + "Strand" + "\t" + "Name" + "\t" + "Synonym" + "\t" + "Product" + "\t" + "Transcription Start" + "\t" + "Translation Start" + "\t" + "Translation Stop" + "\t" + "Transcription Stop" + "\t" + "Strand" + "\t" + "Name" + "\t" + "Synonym" + "\t" + "Product" + "\t" + "Predicted Polycistronic?");
	    for (int i=1; i<genes.size(); i++) {
		int length = getIGlength(genes.get(i-1), genes.get(i));
		writer.print(genes.get(i-1).toString() + "\t" + genes.get(i).toString());
		if (genes.get(i-1).getStrand() != genes.get(i).getStrand()) writer.print("\n");
		else if (length < 40) writer.print("\t" + "YES" + "\n");
		else if (length >= 100) writer.print("\n");
		else if (getCorrelation(genes.get(i-1), genes.get(i)) >= 0.5) writer.print("\t" + "YES" + "\n");
		else writer.print("\n");
	    }
            writer.close();
	} catch (FileNotFoundException e) {
	    Rockhopper.output("\nError - could not open file " + operonOutputFile + "\n\n");
	    System.exit(0);
	}
    }

    /**
     * Output to file merged operons.
     * Return the number of predicted merged operons.
     */
    public int outputMergedOperons(String operonOutputFile, String browserOutputFile, String genomeName, ArrayList<Gene> genes, int size) {
	int numMergedOperons = 0;
	ArrayList<Integer> operonCoordinates = new ArrayList<Integer>(size);
	try {
	    PrintWriter writer = new PrintWriter(new File(operonOutputFile));
	    writer.println("Start" + "\t" + "Stop" + "\t" + "Strand" + "\t" + "Number of Genes" + "\t" + "Genes");
	    ArrayList<String> genesInOperon = new ArrayList<String>();
	    int start = -1;
	    int stop = -1;
	    char strand = '?';
	    for (int i=1; i<genes.size(); i++) {
		if (isGenePairAnOperon(genes, i-1, i)) {  // Operon
		    if (start == -1) {  // Start of new operon
			start = Math.min(genes.get(i-1).getFirst(), genes.get(i).getFirst());
			stop = Math.max(genes.get(i-1).getLast(), genes.get(i).getLast());
			strand = genes.get(i).getStrand();
			genesInOperon.clear();
			genesInOperon.add(genes.get(i-1).getName());
			genesInOperon.add(genes.get(i).getName());
		    } else if (strand == genes.get(i).getStrand()) {  // Within operon
			start = Math.min(start, genes.get(i).getFirst());
			stop = Math.max(stop, genes.get(i).getLast());
			genesInOperon.add(genes.get(i).getName());
		    } else {  // New operon
			writer.print(start + "\t" + stop + "\t" + strand + "\t" + genesInOperon.size() + "\t");
			for (int z = 0; z<genesInOperon.size(); z++) {
			    if (z == 0) writer.print(genesInOperon.get(z));
			    else writer.print(", " + genesInOperon.get(z));
			}
			writer.println();
			numMergedOperons++;
			for (int j=start; j<=stop; j++) {  // Keep track of multi-gene operon coords
			    while (operonCoordinates.size() <= stop) operonCoordinates.add(0);
			    if (strand == '+') operonCoordinates.set(j, 1);
			    if (strand == '-') operonCoordinates.set(j, -1);
			}
			start = -1;
			stop = -1;
			strand = '?';
			genesInOperon.clear();
			genesInOperon.add(genes.get(i-1).getName());
			genesInOperon.add(genes.get(i).getName());
		    }
		} else {  // Non-operon
		    if (start >= 0) {  // End of operon
			writer.print(start + "\t" + stop + "\t" + strand + "\t" + genesInOperon.size() + "\t");
			for (int z = 0; z<genesInOperon.size(); z++) {
			    if (z == 0) writer.print(genesInOperon.get(z));
			    else writer.print(", " + genesInOperon.get(z));
			}
			writer.println();
			numMergedOperons++;
			for (int j=start; j<=stop; j++) {  // Keep track of multi-gene operon coords
			    while (operonCoordinates.size() <= stop) operonCoordinates.add(0);
			    if (strand == '+') operonCoordinates.set(j, 1);
			    if (strand == '-') operonCoordinates.set(j, -1);
			}
			start = -1;
			stop = -1;
			strand = '?';
			genesInOperon.clear();
		    } else {  // Within non-operon
			// Do nothing.
		    }
		}
	    }
	    if (start >= 0) {  // Include final operon
		writer.print(start + "\t" + stop + "\t" + strand + "\t" + genesInOperon.size() + "\t");
		for (int z = 0; z<genesInOperon.size(); z++) {
		    if (z == 0) writer.print(genesInOperon.get(z));
		    else writer.print(", " + genesInOperon.get(z));
		}
		writer.println();
		numMergedOperons++;
		for (int j=start; j<=stop; j++) {  // Keep track of multi-gene operon coords
		    while (operonCoordinates.size() <= stop) operonCoordinates.add(0);
		    if (strand == '+') operonCoordinates.set(j, 1);
		    if (strand == '-') operonCoordinates.set(j, -1);
		}
		start = -1;
		stop = -1;
		strand = '?';
		genesInOperon.clear();
	    }
            writer.close();
	} catch (FileNotFoundException e) {
	    Rockhopper.output("\nError - could not open file " + operonOutputFile + "\n\n");
	    System.exit(0);
	}

	// Output genome browser file with multi-gene operon information
	try {
	    PrintWriter browserWriter = new PrintWriter(new File(browserOutputFile));
	    browserWriter.println("track name=" + "\"" + "Multi-gene Operons" + "\"" + " color=255,0,255 altColor=255,0,255 graphType=bar viewLimits=-1:1");
	    browserWriter.println("fixedStep chrom=" + genomeName + " start=1 step=1");
	    for (int j=1; j<operonCoordinates.size(); j++) browserWriter.println(operonCoordinates.get(j));
	    browserWriter.close();
	} catch (FileNotFoundException e) {
	    Rockhopper.output("\nError - could not open file " + browserOutputFile + "\n\n");
	    System.exit(0);
	}

	return numMergedOperons;
    }



    /**************************************************
     **********   PRIVATE INSTANCE METHODS   **********
     **************************************************/

    /**
     * Determine the distribution of lengths between gene pairs.
     * One distribution for consecutive genes on the same strand and
     * one distribution for consecutive genes on the opposite strand.
     */
    private void determineOperonLengthDistributions(ArrayList<Gene> genes) {
	ArrayList<Integer> length_same = new ArrayList<Integer>();
	ArrayList<Integer> length_opp = new ArrayList<Integer>();
	for (int i=1; i<genes.size(); i++) {
	    int length = Math.min(genes.get(i).getStart(),genes.get(i).getStop()) - Math.max(genes.get(i-1).getStart(),genes.get(i-1).getStop()) - 1;
	    if (genes.get(i).getStrand() == genes.get(i-1).getStrand()) length_same.add(length);
	    else length_opp.add(length);
	}

	// Uncomment the below lines to output operon length distributions
	/*
	lengthSame = new SmoothDistribution(length_same, 30.0, 10, -100, 300);
	lengthOpp = new SmoothDistribution(length_opp, 30.0, 10, -100, 300);
	Rockhopper.output(lengthSame.toString() + "\n" + lengthOpp.toString() + "\n");
	*/

	lengthSame = new SmoothDistribution(length_same, 30.0, 1, -50, 200);
	lengthOpp = new SmoothDistribution(length_opp, 30.0, 1);
	lengthSame.setPseudocount(0.0);
    }

    /**
     * Determine the distribution of correlations between gene pairs.
     * One distribution for consecutive genes on the same strand and
     * one distribution for consecutive genes on the opposite strand.
     */
    private void determineOperonCorrelationDistributions(ArrayList<Gene> genes) {
	ArrayList<Integer> corr_same = new ArrayList<Integer>();
	ArrayList<Integer> corr_opp = new ArrayList<Integer>();
	for (int i=1; i<genes.size(); i++) {
	    double corr = getCorrelation(genes.get(i-1), genes.get(i));
	    if (genes.get(i).getStrand() == genes.get(i-1).getStrand()) corr_same.add((int)(corr*20));
	    else corr_opp.add((int)(corr*20));
	}

	// Uncomment the below lines to output operon correlation distributions
	/*
	corrSame = new SmoothDistribution(corr_same, 5.0, 1, -20, 20);
	corrOpp = new SmoothDistribution(corr_opp, 5.0, 1, -20, 20);
	Rockhopper.output(corrSame.toString() + "\n" + corrOpp.toString() + "\n");
	*/

	corrSame = new SmoothDistribution(corr_same, 5.0, 1, -20, 20);
	corrOpp = new SmoothDistribution(corr_opp, 5.0, 1, -20, 20);
	corrSame.setPseudocount(0.0);
    }

    private double getCorrelation(Gene g1, Gene g2) {
	ArrayList<Long> e1 = new ArrayList<Long>();
	ArrayList<Long> e2 = new ArrayList<Long>();
	for (int i=0; i<Rockhopper.numConditions; i++) {
	    e1.add(g1.getAvg(i));
	    e2.add(g2.getAvg(i));
	}
	return Misc.correlation(e1, e2);
    }

    /**
     * Set prior probability of same-strand gene pair being an operon
     * (based on section 2.2 of Westover 2005).
     */
    private void setOperonPrior(ArrayList<Gene> genes) {
	int numberDirectons = 0;  // Two or more consecutive genes on the same strand
	int numberPairs = 0;  // Pair of consecutive genes on the same strand
	char previousStrand = '?';
	int currentDirectonLength = 1;  // Length of current directon under consideration
	for (int i=0; i<genes.size(); i++) {
	    Gene g = genes.get(i);
	    if ((g.getStrand() != previousStrand) && (currentDirectonLength > 1)) {
		numberDirectons += 1;
		previousStrand = g.getStrand();
		currentDirectonLength = 1;
	    } else if ((g.getStrand() != previousStrand) && (currentDirectonLength == 1)) {
		previousStrand = g.getStrand();
		currentDirectonLength = 1;
	    } else if ((g.getStrand() == previousStrand) && (currentDirectonLength > 1)) {
		numberPairs += 1;
		previousStrand = g.getStrand();
		currentDirectonLength += 1;
	    } else if ((g.getStrand() == previousStrand) && (currentDirectonLength == 1)) {
		numberPairs += 1;
		previousStrand = g.getStrand();
		currentDirectonLength += 1;
	    } else {
		// Do nothing
	    }
	}
	if (currentDirectonLength > 1) numberDirectons += 1;
	operonPrior = 1.0 - numberDirectons/(double)numberPairs;
	nonOperonPrior = 1.0 - operonPrior;
	//Rockhopper.output("Prior probability of operon:\t" + operonPrior + "\n");  // Output operon prior
    }

    /**
     * Compute probability that two consecutive genes are co-transcribed.
     * Return a list of p-values for all gene pairs.
     */
    private ArrayList<Double> operonExpression(ArrayList<Gene> genes) {
	ArrayList<Double> p_values = new ArrayList<Double>(genes.size());
	for (int i=0; i<genes.size(); i++) p_values.add(Double.MAX_VALUE);
	p_values.set(0, 0.0);
	ArrayList<ArrayList<Double>> all_p_values = new ArrayList<ArrayList<Double>>();
	for (int i=0; i<genes.size(); i++) all_p_values.add(new ArrayList<Double>());
	for (int i=0; i<Rockhopper.numConditions; i++) {
	    ArrayList<Double> means = new ArrayList<Double>(genes.size());
	    ArrayList<Double> variances = new ArrayList<Double>(genes.size());
	    for (int z=0; z<genes.size(); z++) {
		means.add(0.0);
		variances.add(0.0);
	    }
	    for (int z=1; z<genes.size(); z++) {
		int numReplicates = genes.get(z).getNumReplicates(i);

		// Compute mean expression in each condition.
		for (int j=0; j<numReplicates; j++) {
		    means.set(z, means.get(z) + 1000.0*genes.get(z).getNormalizedCount(i,j)/(Math.abs(genes.get(z).getStop()-genes.get(z).getStart())+1));
		}
		means.set(z, means.get(z) / numReplicates);

		// Compute variance.
		double varianceAdjustment = 1.15;
		// If we have NO replicates, then we use neighboring genes.
		if (numReplicates == 1) {  // We have no replicates. Use neighboring gene
		    double mean = (1000.0*genes.get(z-1).getNormalizedCount(i,0)/(Math.abs(genes.get(z-1).getStop()-genes.get(z-1).getStart())-1) + 1000.0*genes.get(z).getNormalizedCount(i,0)/(Math.abs(genes.get(z).getStop()-genes.get(z).getStart())+1)) / 2.0;
		    variances.set(z, Math.pow(1000.0*genes.get(z-1).getNormalizedCount(i,0)/(Math.abs(genes.get(z-1).getStop()-genes.get(z-1).getStart())+1) - mean, 2.0) + Math.pow(1000.0*genes.get(z).getNormalizedCount(i,0)/(Math.abs(genes.get(z).getStop()-genes.get(z).getStart())+1) - mean, 2.0));
		    variances.set(z, variances.get(z) / (2-1));
		    variances.set(z, Math.pow(variances.get(z), varianceAdjustment));
		}

		// Compute variance.
		// If we DO have replicates, then we use the replicates.
		if (numReplicates > 1) {  // We have replicates. Use them
		    for (int j=0; j<numReplicates; j++) {
			variances.set(z, Math.pow(1000.0*genes.get(z).getNormalizedCount(i,j)/(Math.abs(genes.get(z).getStop()-genes.get(z).getStart())+1) - means.get(z), 2.0));
		    }
		    variances.set(z, variances.get(z) / (numReplicates-1));
		    variances.set(z, Math.pow(variances.get(z), varianceAdjustment));
		}
	    }

	    // Generate Lowess variances
	    means.set(0, means.get(1));
	    variances.set(0, variances.get(1));
	    ArrayList<Long> means_Long = new ArrayList<Long>(means.size());
	    ArrayList<Long> variances_Long = new ArrayList<Long>(variances.size());
	    for (int w=0; w<means.size(); w++) {
		means_Long.add((long)(double)means.get(w));
		variances_Long.add((long)(double)variances.get(w));
	    }
	    ArrayList<Long> lowessVariance = Lowess.lowess(means_Long, variances_Long);

	    // Determine operon probabilities
	    for (int z=1; z<genes.size(); z++) {
		double p_value = computeProbabilityOfSameTranscript(genes.get(z-1), genes.get(z), (double)lowessVariance.get(z-1), (double)lowessVariance.get(z), i);
		p_values.set(z, Math.min(p_values.get(z), p_value));
		all_p_values.get(z).add(p_value);
	    }
	}
	for (int z=1; z<genes.size(); z++) {
	    double avg = 0.0;
	    for (int y=0; y<all_p_values.get(z).size(); y++) avg += all_p_values.get(z).get(y);
	    avg /= all_p_values.get(z).size();
	    p_values.set(z, avg);
	}
	return p_values;
    }



    /***********************************************
     **********   PRIVATE CLASS METHODS   **********
     ***********************************************/

    /**
     * Reads in a file of genes (either *.ptt or *.rnt) and returns
     * an ArrayList of gene objects.
     */
    private static ArrayList<Gene> readInGenes(String fileName) {
        ArrayList<Gene> listOfGenes = new ArrayList<Gene>();
        try {
            Scanner reader = new Scanner(new File(fileName));
            for (int i=0; i<3; i++) reader.nextLine();  // Ignore 3 header lines
            while (reader.hasNext()) {  // Continue until end of file
                listOfGenes.add(new Gene(reader.nextLine(), "ORF"));  // Create new gene
            }
            reader.close();
        } catch (FileNotFoundException e) {
            Rockhopper.output("Error - the file " + fileName + " could not be found and opened.\n");
        }
        return listOfGenes;
    }

    /**
     * Returns the probability, based on expression in the given condition "i", that "g1"
     * is NOT differentially expressed from "g2".
     * I.e., return the probability that "g1" is part of the same polycistronic transcript
     * as "g2" in the specified condition "i".
     */
    private static double computeProbabilityOfSameTranscript(Gene g1, Gene g2, double lowessVar1, double lowessVar2, int i) {
	int numReplicates = g1.getNumReplicates(i);

	double k_A = 0.0;
	double k_B = 0.0;
	for (int j=0; j<numReplicates; j++) {
	    k_A += 1000.0*g1.getNormalizedCount(i,j)/(Math.abs(g1.getStop()-g1.getStart())+1);
	    k_B += 1000.0*g2.getNormalizedCount(i,j)/(Math.abs(g2.getStop()-g2.getStart())+1);
	}
	double q = k_A + k_B;

	double mean_A = q;
	double mean_B = q;
	double variance_A = lowessVar1;
	double variance_B = lowessVar2;
	double p_a = mean_A / variance_A;
	double p_b = mean_B / variance_B;
	double r_a = Math.max(mean_A*mean_A / (variance_A - mean_A), 1.0);  // r should never be below 1
	double r_b = Math.max(mean_B*mean_B / (variance_B - mean_B), 1.0);  // r should never be below 1

	if ((p_a < 0.0) || (p_b < 0.0) || (p_a > 1.0) || (p_b > 1.0) || (variance_A == 0.0) || (variance_B == 0.0)) return 0.0;

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
		double average_p = (current_p + previous_p) / 2.0;
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
	alpha = 1000;  // Number of times we decrement by 1 (raising alpha raises precision but slows down computation)
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
	return p_value;
    }

    /**
     * Returns the number of nucleotides in the IG region between the
     * two genes. g1 precede g2. May return a negative number if the
     * genes overlap.
     */
    private static int getIGlength(Gene g1, Gene g2) {
	int IG_start = -1;
	if (g1.isORF())  IG_start = Math.max(g1.getStart(), g1.getStop()) + 1;
	else IG_start = Math.max(g1.getStartT(), g1.getStopT()) + 1;
	int IG_stop = -1;
	if (g2.isORF())  IG_stop = Math.min(g2.getStart(), g2.getStop()) - 1;
	else IG_stop = Math.min(g2.getStartT(), g2.getStopT()) - 1;
	return IG_stop - IG_start + 1;
    }

    /**
     * Returns the percentage of nucleotides in the IG region
     * that are expressed, i.e., that correspond to predicted
     * UTRs.
     */
    private static double getPercentIGexpressed(Gene g1, Gene g2) {
	if (!g1.isORF() || !g2.isORF()) return 0.0;  // Only consider two ORFs
	int IG_length = getIGlength(g1, g2);
	if (IG_length <= 0) return 0.0;  // Only consider IGs with non-zero length

	int UTR1 = 0;
	if ((g1.getStrand() == '+') && (g1.getStopT() > 0))
	    UTR1 = g1.getStopT() - g1.getStop();
	if ((g1.getStrand() == '-') && (g1.getStartT() > 0))
	    UTR1 = g1.getStartT() - g1.getStart();
	int UTR2 = 0;
	if ((g2.getStrand() == '+') && (g2.getStartT() > 0))
	    UTR2 = g2.getStart() - g2.getStartT();
	if ((g2.getStrand() == '-') && (g2.getStopT() > 0))
	    UTR2 = g2.getStop() - g2.getStopT();
	int UTR_length = UTR1 + UTR2;
	return Math.min(UTR_length / (double)IG_length, 1.0);
    }



    /*************************************
     **********   MAIN METHOD   **********
     *************************************/

    /**
     * The Main method, when invoked with the name of a gene file (*.ptt)
     * as a command line argument, computes the distribution of lengths between
     * consecutive genes on the same strand and the distribution of lengths
     * between consecutive genes on the opposite strand. Output is to same.dist
     * and opp.dist.
     */
    public static void main(String[] args) {
	if (args.length < 1) {
	    System.err.println("\nThe Operons application requires one command line argument, the name of a gene file (*.ptt). The application computes the distribution of lengths between consecutive genes on the same strand and the distribution of lengths between consecutive genes on the opposite strand. Output is to same.dist and opp.dist.");
	    System.err.println("\n\tjava Operons NC_******.ptt\n");
	    System.exit(0);
	}

	Operons ops = new Operons(args[0]);

    }

}
