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
import java.io.File;
import java.io.FileNotFoundException;

/**
 * A Replicate object represents information about a 
 * single RNA-seq experiment, including information
 * about all reads from the experiment.
 */
public class Replicate {

    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private String fileName;  // Path to sequencing reads file
    private ArrayList<String> compressedFileNames;  // List of compressed (WIG) file names
    private String name;      // Name of experiment
    private int[][] plusReads;   // Reads on the plus strand for each genome
    private int[][] minusReads;  // Reads on the minus strand for each genome
    private long totalReads;
    private long upperQuartile;  // Number of reads for 75% most expressed gene
    private ArrayList<Integer> background;  // Number of nucleotides with few or no reads
    private double backgroundParameter;  // Parameter of geometric distribution of background reads
    private double avgReads;  // Average number of reads mapping to nucleotides throughout genome
    private double minExpressionUTR;  // Minimum expression for a UTR region to be deemed expressed
    private double minExpressionRNA;  // Minimum expression for a ncRNA to be deemed expressed
    private long avgLengthReads;  // Avg length of mapping reads in this replicate file



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    /**
     * Constructs a new Replicate object based on compressed sequencing reads files.
     */
    public Replicate(ArrayList<String> fileNames, String readFile, boolean unstranded) {
	compressedFileNames = new ArrayList<String>();
	plusReads = new int[Peregrine.numSequences][];
	minusReads = new int[Peregrine.numSequences][];
	this.totalReads = 0;
	int backgroundLength = 2;
	this.background = new ArrayList<Integer>(backgroundLength);
	for (int i=0; i<backgroundLength; i++) background.add(0);
	for (int z=0; z<Peregrine.numSequences; z++) {
	    String compressedFileName = fileNames.get(z);
	    compressedFileNames.add(compressedFileName);
	    readInAlignmentFile(z, compressedFileName, unstranded);
	}
	for (int i=0; i<backgroundLength-1; i++) backgroundParameter += background.get(i);
	backgroundParameter /= (3.0*2.0*Rockhopper.genomeSize);
	if (unstranded) backgroundParameter /= (3.0*1.0*Rockhopper.genomeSize);
	avgReads = totalReads / (2.0 * Rockhopper.genomeSize);
	if (unstranded) avgReads = totalReads / (1.0 * Rockhopper.genomeSize);
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    /**
     * Return the name of this Replicate.
     */
    public String getName() {
	return this.name;
    }

    /**
     * Return the name of the compressed WIG file for this Replicate
     * for the specified genome at index z.
     */
    public String getCompressedFileName(int z) {
	return this.compressedFileNames.get(z);
    }

    /**
     * Return total number of reads in this Replicate.
     */
    public long getTotalReads() {
	return this.totalReads;
    }

    /**
     * Return average number of reads in this Replicate.
     */
    public double getAvgReads() {
	return this.avgReads;
    }

    /**
     * Return the number of reads mapping to the specified coordinate
     * on the specified strand for the specified genome at index z.
     */
    public int getReads(int z, int coord, char strand) {
	if (strand == '+') return plusReads[z][coord];
	else if (strand == '-') return minusReads[z][coord];
	else return plusReads[z][coord] + minusReads[z][coord];
    }

    /**
     * Return the number of reads on the given strand mapping
     * to the given range of genomic coordinates for the specified
     * genome at index z.
     */
    public long getReadsInRange(int z, int start, int stop, char strand) {
	if (stop < start) {  // Swap
	    int temp = start;
	    start = stop;
	    stop = temp;
	}

	long sum = 0;
	for (int i=Math.max(start, 1); i<Math.min(stop+1, plusReads[z].length); i++) {
	    if (strand == '+')  // Plus strand
		sum += this.plusReads[z][i];
	    else if (strand == '-')  // Minus strand
		sum += this.minusReads[z][i];
	    else  // Ambiguous strand
		sum += this.plusReads[z][i] + this.minusReads[z][i];
	}
	return sum;
    }

    /**
     * Return the upper quartile for this Replicate.
     */
    public long getUpperQuartile() {
	return upperQuartile;
    }

    /**
     * Set the upper quartile, i.e., the number of reads mapping
     * to the 75% most expressed gene (not counting genes with
     * zero reads mapping to them).
     */
    public void setUpperQuartile(long upperQuartile) {
	this.upperQuartile = upperQuartile;
    }

    /**
     * Sets the minimum level of expression (for a UTR region and
     * ncRNA to be considered expressed) in this Replicate based on 
     * the average number of reads per nucleotide in this Replicate
     * and the specified transcript sensitivity between 0.0 and 1.0, 
     * inclusive.
     */
    public void setMinExpression(double transcriptSensitivity) {
	minExpressionUTR = transformation(Math.pow(transcriptSensitivity, 3.3));
	minExpressionRNA = transformation(Math.pow(transcriptSensitivity, 0.4));
    }

    /**
     * Return the minimum level of expression for a UTR region to be
     * considered expressed in this Replicate.
     */
    public double getMinExpressionUTR() {
	return this.minExpressionUTR;
    }

    /**
     * Return the minimum level of expression for a ncRNA to be
     * considered expressed in this Replicate.
     */
    public double getMinExpressionRNA() {
	return this.minExpressionRNA;
    }

    /**
     * Return average length of sequencing reads in this Replicate.
     */
    public long getAvgLengthReads() {
	return this.avgLengthReads;
    }

    /**
     * Returns the probability that the given number of reads
     * at some nucleotide corresponds to the background, i.e.,
     * a non-transcript.
     */
    public double getBackgroundProb(int numReads) {
	// Based on geometric distribution
	return Math.pow(1.0 - this.backgroundParameter, numReads) * this.backgroundParameter;
    }

    /**
     * Return the mean number of reads on the given strand mapping
     * to the given range of genomic coordinates in the specified
     * genome at index z.
     */
    public double getMeanOfRange(int z, int start, int stop, char strand) {
	return getReadsInRange(z, start, stop, strand) / (double)(Math.abs(stop-start)+1);
    }

    /**
     * Return the standard deviation of reads on the given strand mapping
     * to the given range of genomic coordinates in the specified genome
     * at index z.
     */
    public double getStdevOfRange(int z, int start, int stop, char strand, double mean) {
	if (stop < start) {  // Swap
	    int temp = start;
	    start = stop;
	    stop = temp;
	}

	double stdev = 0.0;
	for (int i=Math.max(start, 1); i<Math.min(stop+1, plusReads[z].length); i++) {
	    if (strand == '+')  // Plus strand
		stdev += Math.pow(plusReads[z][i] - mean, 2.0);
	    else if (strand == '-')  // Minus strand
		stdev += Math.pow(minusReads[z][i] - mean, 2.0);
	    else  // Ambiguous strand
		stdev += Math.pow(plusReads[z][i] + minusReads[z][i] - mean, 2.0);
	}
	return stdev / Math.sqrt(Math.min(stop, plusReads[z].length-1) - Math.max(start, 1) + 1);
    }

    /**
     * Returns a String representation of this object.
     */
    public String toString() {
	return name;
    }



    /**************************************************
     **********   PRIVATE INSTANCE METHODS   **********
     **************************************************/

    /**
     * Reads in a compressed alignment file for a genome with
     * the specified index z. Stores the 
     * RNA-seq data in two lists, one for the plus strand 
     * and one for the minus strand. Also computes the 
     * total reads in the file.
     */
    private void readInAlignmentFile(int z, String inFile, boolean unstranded) {
	FileOps foo = new FileOps(inFile);
	if (!foo.isValid()) {
	    Rockhopper.output("Error - could not read in file " + inFile + "\n");
	    return;
	}

	this.fileName = foo.getReadFileName();
	int slashIndex = fileName.lastIndexOf(File.separator);
	int periodIndex = fileName.lastIndexOf(".");
	if (periodIndex < 0) periodIndex = fileName.length();
	this.name = fileName.substring(slashIndex+1, periodIndex);
	this.avgLengthReads = foo.getAvgLengthReads();
	int[] coordinates_plus = foo.getCoordinates_plus();
	int[] coordinates_minus = foo.getCoordinates_minus();
	this.plusReads[z] = new int[Rockhopper.genomeSizes.get(z)];
	this.minusReads[z] = new int[Rockhopper.genomeSizes.get(z)];
	for (int i=0; i<Rockhopper.genomeSizes.get(z); i++) {
	    int reads_plus = coordinates_plus[i];
	    int reads_minus = coordinates_minus[i];
	    plusReads[z][i] = reads_plus;
	    minusReads[z][i] = reads_minus;
	    totalReads += reads_plus + reads_minus;
	    if (!unstranded) {  // Strand specific
		if (reads_plus < background.size()) background.set(reads_plus, background.get(reads_plus) + 1);
		else background.set(background.size()-1, background.get(background.size()-1) + 1);
		if (reads_minus < background.size()) background.set(reads_minus, background.get(reads_minus) + 1);
		else background.set(background.size()-1, background.get(background.size()-1) + 1);
	    } else {  // Strand ambiguous
		if (reads_plus + reads_minus < background.size()) background.set(reads_plus + reads_minus, background.get(reads_plus + reads_minus) + 1);
		else background.set(background.size()-1, background.get(background.size()-1) + 1);
	    }
	}
    }

    /**
     * Helper method for setting the minimum expression level for UTRs and ncRNAs.
     */
    private double transformation(double transcriptSensitivity) {
	if (transcriptSensitivity <= 0.5) {  // less than or equal to 0.5
	    return Math.pow(avgReads/2, Math.pow(transcriptSensitivity/0.5, 0.25));
	} else {  // greater than 0.5
	    return Math.pow(avgReads/2, 1 + Math.pow(2.0*(transcriptSensitivity - 0.5), 2));
	}
    }



    /*************************************
     **********   MAIN METHOD   **********
     *************************************/

    /**
     * The main method is used to test the methods of this class.
     */
    public static void main(String[] args) {
	/*
	if (args.length == 0) {
	    System.err.println("\nThe Replicate application must be invoked with a single command line argument corresponding to the name of a WIG file.\n");
	    System.exit(0);
	}

	Replicate r = new Replicate(args[0], false);
	System.out.println(r);
	*/
    }
}

