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

import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;

public class Replicon {

    /********************************************
     **********   Instance Variables   **********
     ********************************************/

    private String name;      // Name of replicon (first line of FASTA file)
    private String sequence;  // Original DNA sequence with '$' appended
    private int length;       // Length of sequence
    private int[] rotations;  // Burrows-Wheeler matrix of cyclic rotations
    private String BWT;       // Borrows-Wheeler transform of original sequence

    private int[] C;          // Number of chars lexicographically less than...
    private int[][] occ;      // Number of occurrences of char up to row...



    /**************************************
     **********   Constructors   **********
     **************************************/

    public Replicon(String fastaFile) {
	this.sequence = readInFastaFile(fastaFile);
	this.sequence = this.sequence.toUpperCase();
	this.sequence = this.sequence.replace('U', 'T');
	this.sequence = replaceAmbiguousCharacters();
	this.sequence = this.sequence + "$";
	this.length = sequence.length();
	this.rotations = rotations();
	this.BWT = getTransform();

	precomputeCharacterInfo();  // Populate "C" and "occ" instance variables
    }



    /*************************************************
     **********   Public Instance Methods   **********
     *************************************************/

    // Return a row of the cyclic rotations matrix.
    public String getStringInRotationsMatrix(int row) {
	return sequence.substring(rotations[row]) + sequence.substring(0, rotations[row]);
    }

    // Reverse the BWT transform
    public String unpermute() {
	StringBuilder sb = new StringBuilder();
	int r = 0;
	while (BWT.charAt(r) != '$') {
	    sb.append(BWT.charAt(r));
	    r = stepLeft(r);
	}
	return sb.reverse().toString();
    }

    /**
     * Returns the name of the genome sequence (FASTA header line).
     */
    public String getName() {
	return this.name;
    }

    /**
     * Returns the DNA sequence of this Replicon.
     */
    public String getSequence() {
	return this.sequence;
    }

    /**
     * Returns the length of the genome sequence (plus 1).
     */
    public int getLength() {
	return this.length;
    }

    /**
     * Returns number of characters lexographically less than x.
     */
    public int getC(int x) {
	return this.C[x];
    }

    /**
     * Returns number of occurrences of character x up to row y.
     */
    public int getOcc(int x, int y) {
	return this.occ[x][y];
    }

    /**
     * Returns index in Burrows-Wheeler rotations matrix at row x.
     */
    public int getRotations(int x) {
	return this.rotations[x];
    }



    /**************************************************
     **********   Private Instance Methods   **********
     **************************************************/

    // Read in sequence from from file.
    private String readInFastaFile(String fastaFile) {
	try {
	    Scanner reader = new Scanner(new File(fastaFile));
	    StringBuilder sb = new StringBuilder();

	    // Handle header line
	    if (reader.hasNextLine()) {
		String line = reader.nextLine();
		if (line.startsWith(">")) {  // Header line
		    this.name = line.substring(1).trim();
		} else {  // Not header line
		    sb.append(line);
		}
	    }

	    // Read in file
	    while (reader.hasNextLine()) {
		sb.append(reader.nextLine());
	    }

	    reader.close();
	    return sb.toString();
	} catch (FileNotFoundException e) {
	    Peregrine.output("\nError - could not open file " + fastaFile + "\n\n");
	    System.exit(0);
	}
	return "";
    }

    // Replace any non-nucleotide character with '^'.
    // We use '^' because its ASCII value is greater than A,C,G,T
    private String replaceAmbiguousCharacters() {
	StringBuilder sb = new StringBuilder(this.sequence);
	for (int i=0; i<sb.length(); i++) {
	    if ((sb.charAt(i) != 'A') && (sb.charAt(i) != 'C') &&
		(sb.charAt(i) != 'G') && (sb.charAt(i) != 'T'))
		sb.setCharAt(i, '^');
	}
	return sb.toString();
    }

    // Compute the Burrows-Wheeler matrix of cyclic rotations.
    // Rather than return a matrix, we return an array of integers
    // where each value in the array is an index into the sequence.
    // So a "row" of the matrix corresponds to the cyclic sequence
    // beginnning at the specified index. All rows of the matrix
    // have length equal to the length of the sequence.
    private int[] rotations() {
	int[] indices= new int[this.length];
	for (int i=0; i<length; i++) indices[i] = i;
	quicksort(indices, 0, length-1);
	return indices;
    }

    // Compute the Burrows-Wheeler transform from the "rotations"
    // matrix, i.e., the last column of the matrix.
    private String getTransform() {
	StringBuilder sb = new StringBuilder();
	for (int i=0; i<length; i++) {
	    sb.append(sequence.charAt((rotations[i]-1+length) % length));;
	}
	return sb.toString();
    }

    // Precompute information about how often each character occurs in BWT
    private void precomputeCharacterInfo() {

	// Populate instance variable "C"
	this.C = new int[6];  // All characters in DNA sequence have
	                        // ASCII value less than 123.
	for (int i=0; i<length; i++) C[Index.charToInt(BWT.charAt(i))] += 1;

	// Cumulative
	int[] temp = new int[6];
	temp[0] = C[0] + C[5];  // Number of characters <= 'A'
	temp[5] = C[5];         // Number of characters <= '$'
	for (int i=1; i<5; i++) temp[i] = C[i] + temp[i-1];

	// Strictly less than
	for (int i=0; i<6; i++) C[i] = temp[i] - C[i];

	// Populate instance variable "occ"
	this.occ = new int[5][length+1];
	for (int i=0; i<length; i++) {
	    for (int j=0; j<5; j++) {
		occ[j][i+1] = occ[j][i];
	    }
	    int index = Index.charToInt(BWT.charAt(i));
	    if ((index >= 0) && (index < 5)) occ[index][i+1]++;
	}
    }

    // Helper method for Burrows-Wheeler transform.
    // Taken from appendix of Bowtie article.
    private int stepLeft(int r) {
	int index = Index.charToInt(BWT.charAt(r));
	return C[index] + occ[index][r];
    }

    // We quicksort the array. But we don't sort the values at
    // each array index, since each value is an index 
    // corresponding to a String. We sort the Strings specified
    // by each array index.
    private void quicksort(int[] a, int p, int r) {
	if (p < r) {
	    int q = partition(a, p, r);
	    quicksort(a, p, q-1);
	    quicksort(a, q+1, r);
	}
    }

    // Quicksort helper method
    private int partition(int[] a, int p, int r) {
	int x = Index.rand.nextInt(r - p + 1) + p;
	swap(a, x, r);
	x = a[r];
	int i = p-1;
	for (int j=p; j<r; j++) {
	    if (lessThanOrEqualTo(a[j], x)) {  // a[j] <= x
		i++;
		swap(a, i, j);
	    }
	}
	swap(a, i+1, r);
	return i+1;
    }

    // Quicksort partition helper method
    private void swap(int a[], int i, int j) {
	int temp = a[i];
	a[i] = a[j];
	a[j] = temp;
    }

    // Given two indices in "sequence", determine whether the
    // cyclic String corresponding to the first index is less 
    // than or equal to the cyclic String corresponding to
    // the second index.
    private boolean lessThanOrEqualTo(int i, int j) {
	for (int z=0; z<length; z++) {
	    if (sequence.charAt((i+z)%length) < sequence.charAt((j+z)%length))
		return true;
	    if (sequence.charAt((i+z)%length) > sequence.charAt((j+z)%length))
		return false;

	}
	return true;  // We have equality of the two Strings
    }



    /*************************************
     **********   Main Method   **********
     *************************************/

    public static void main(String[] args) {

        if (args.length < 1) {
            System.err.println("\nUSAGE: java Replicon <index.fna>" + "\n");
            System.err.println("Replicon takes a sequence file (such as a FASTA genome) to be indexed and it computes the Burrows-Wheeler transformed index. If run from the command line, it simply outputs the size of the index to STDOUT. The Replicon class is meant to be instantiated from another application. \n");
            System.exit(0);
        }

	Replicon index = new Replicon(args[0]);
	System.out.println("Size of index is:\t" + index.BWT.length());
    }

}

