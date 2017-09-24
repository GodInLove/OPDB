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
import java.util.concurrent.atomic.AtomicLongArray;

public class DeNovoIndex {

    /********************************************
     **********   Instance Variables   **********
     ********************************************/

    public String sequence;  // Concatenated transcripts (separated by '^' with '$' at end)
    private int length;       // Length of sequence
    private int[] rotations;  // Burrows-Wheeler matrix of cyclic rotations
    private String BWT;       // Borrows-Wheeler transform of original sequence
    private boolean stopAfterOneHit;  // Report only one mapping as opposed to all mappings
    public AtomicIntegerArray[] readCounts;  // Number of reads mapping to each nt of each transcript from each sequencing reads file
    private AtomicIntegerArray numReads;         // Num reads from each file
    private AtomicIntegerArray numMappingReads;  // Num reads from each file mapping to transcripts
    private AtomicLongArray avgLengthReads;      // Avg length of mapping reads from each file

    private int[] C;          // Number of chars lexicographically less than...
    private int[][] occ;      // Number of occurrences of char up to row...



    /*****************************************
     **********   Class Variables   **********
     *****************************************/

    private static Random rand = new Random();



    /**************************************
     **********   Constructors   **********
     **************************************/

    public DeNovoIndex(ArrayList<StringBuilder> transcripts) {
	StringBuilder sb = new StringBuilder();
	for (StringBuilder transcript : transcripts) {
	    sb.append(transcript);
	    sb.append('^');
	}
	sb.append('$');
	this.sequence = sb.toString();
	this.length = sequence.length();
	this.rotations = rotations();
	this.BWT = getTransform();
	precomputeCharacterInfo();  // Populate "C" and "occ" instance variables
    }



    /*************************************************
     **********   Public Instance Methods   **********
     *************************************************/

    /**
     * Returns true if s is in index, false otherwise.
     */
    public boolean exactMatch(String s) {
	return exactMatch(s, 0, s.length());
    }

    /**
     * Returns true if s[start:end] is in index, false otherwise.
     * Note "start" is inclusive and "end" is exclusive.
     */
    public boolean exactMatch(String s, int start, int end) {
	int i = end - 1;
	char c = s.charAt(i);
	int sp = getC(Dictionary.charToInt(c));
	int ep = -1;
	if (Dictionary.charToInt(c) < 4) ep = getC(Dictionary.charToInt(c)+1);
	else ep = length;
	i--;
	while ((sp < ep) && (i >= start)) {
	    int index = Dictionary.charToInt(s.charAt(i));
	    sp = getC(index) + getOcc(index, sp);
	    ep = getC(index) + getOcc(index, ep);
	    i--;
	}
	if (ep > sp) return true;
	return false;
    }

    /**
     * Tallies the number of full length reads mapping to 
     * each nucleotide of each transcript.
     * If full length read does not map EXACTLY, then the
     * tally is not incremented.
     */
    public void exactMatch_fullRead(String s) {
	numReads.getAndIncrement(Assembler.readsFileIndex);
	int start = 0;
	int end = s.length();
	int i = end - 1;
	char c = s.charAt(i);
	int sp = getC(Dictionary.charToInt(c));
	int ep = -1;
	if (Dictionary.charToInt(c) < 4) ep = getC(Dictionary.charToInt(c)+1);
	else ep = length;
	i--;
	while ((sp < ep) && (i >= start)) {
	    int index = Dictionary.charToInt(s.charAt(i));
	    sp = getC(index) + getOcc(index, sp);
	    ep = getC(index) + getOcc(index, ep);
	    i--;
	}
	if (!stopAfterOneHit) {  // Report all mappings of read
	    for (int j=sp; j<ep; j++) {
		for (int k=getRotations(j); k<getRotations(j)+s.length(); k++) readCounts[Assembler.readsFileIndex].incrementAndGet(k);
	    }
	    if (ep > sp) {
		avgLengthReads.getAndAdd(Assembler.readsFileIndex, s.length());
		numMappingReads.getAndIncrement(Assembler.readsFileIndex);
	    }
	} else {  // Report only one mapping of read
	    if (ep > sp) {
		int randIndex = sp + rand.nextInt(ep-sp);
		for (int k=getRotations(randIndex); k<getRotations(randIndex)+s.length(); k++) readCounts[Assembler.readsFileIndex].incrementAndGet(k);
		avgLengthReads.getAndAdd(Assembler.readsFileIndex, s.length());
		numMappingReads.getAndIncrement(Assembler.readsFileIndex);
	    }
	}
    }

    /**
     * Tallies the number of full length *paired-end* reads mapping to 
     * each nucleotide of each transcript.
     * If full length *paired-end* read does not map EXACTLY, then the
     * tally is not incremented.
     */
    public void exactMatch_fullRead(String s, String s2) {
	numReads.getAndIncrement(Assembler.readsFileIndex);

	// Read 1
	int start = 0;
	int end = s.length();
	int i = end - 1;
	char c = s.charAt(i);
	int sp = getC(Dictionary.charToInt(c));
	int ep = -1;
	if (Dictionary.charToInt(c) < 4) ep = getC(Dictionary.charToInt(c)+1);
	else ep = length;
	i--;
	while ((sp < ep) && (i >= start)) {
	    int index = Dictionary.charToInt(s.charAt(i));
	    sp = getC(index) + getOcc(index, sp);
	    ep = getC(index) + getOcc(index, ep);
	    i--;
	}
	if (sp >= ep) return;

	// Read 2
	start = 0;
	end = s2.length();
	i = end - 1;
	c = s2.charAt(i);
	int sp2 = getC(Dictionary.charToInt(c));
	int ep2 = -1;
	if (Dictionary.charToInt(c) < 4) ep2 = getC(Dictionary.charToInt(c)+1);
	else ep2 = length;
	i--;
	while ((sp2 < ep2) && (i >= start)) {
	    int index = Dictionary.charToInt(s2.charAt(i));
	    sp2 = getC(index) + getOcc(index, sp2);
	    ep2 = getC(index) + getOcc(index, ep2);
	    i--;
	}
	if (sp2 >= ep2) return;

	// Combine paired-end hits
	ArrayList<Integer> hits1 = new ArrayList<Integer>(ep-sp);
	for (int j=sp; j<ep; j++) hits1.add(getRotations(j));
	ArrayList<Integer> hits2 = new ArrayList<Integer>(ep2-sp2);
	for (int j=sp2; j<ep2; j++) hits2.add(getRotations(j));
	ArrayList<Integer> starts = new ArrayList<Integer>();
	ArrayList<Integer> ends = new ArrayList<Integer>();
	combinePairedEndHits(s.length(), s2.length(), hits1, hits2, starts, ends);

	if (!stopAfterOneHit) {  // Report all mappings of read
	    int sum = 0;
	    for (int j=0; j<starts.size(); j++) {
		for (int k=starts.get(j); k<ends.get(j); k++) readCounts[Assembler.readsFileIndex].incrementAndGet(k);
		sum += ends.get(j) - starts.get(j);
	    }
	    if (starts.size() > 0) {
		avgLengthReads.getAndAdd(Assembler.readsFileIndex, sum/starts.size());
		numMappingReads.getAndIncrement(Assembler.readsFileIndex);
	    }
	} else {  // Report only one mapping of read
	    if (starts.size() > 0) {
		int randIndex = rand.nextInt(starts.size());
		for (int k=starts.get(randIndex); k<ends.get(randIndex); k++) readCounts[Assembler.readsFileIndex].incrementAndGet(k);
		avgLengthReads.getAndAdd(Assembler.readsFileIndex, ends.get(randIndex)-starts.get(randIndex));
		numMappingReads.getAndIncrement(Assembler.readsFileIndex);
	    }
	}
    }

    /**
     * Creates infrastructure to store information about
     * full length reads mapping to the length (rather
     * than just k-mers).
     */
    public void initializeReadMapping(boolean stopAfterOneHit) {
	int numReadsFiles = Assembler.readsFileIndex;
	readCounts = new AtomicIntegerArray[numReadsFiles];
	numReads = new AtomicIntegerArray(numReadsFiles);
	avgLengthReads = new AtomicLongArray(numReadsFiles);
	numMappingReads = new AtomicIntegerArray(numReadsFiles);
	for (int i=0; i<numReadsFiles; i++) readCounts[i] = new AtomicIntegerArray(length);
	this.stopAfterOneHit = stopAfterOneHit;
    }

    public int[] getAvgLengthOfReads() {
	int[] avgLengths = new int[avgLengthReads.length()];
	for (int i=0; i<avgLengths.length; i++) {
	    if (numMappingReads.get(i) == 0) avgLengths[i] = 0;
	    else avgLengths[i] = (int)(avgLengthReads.get(i)/numMappingReads.get(i));
	}
	return avgLengths;
    }

    /**
     * In the case of unstranded (strand ambiguous) reads,
     * we invoke exactMatch_fullRead twice, once per strand,
     * thereby double counting each read. Here we fix the
     * double dip.
     */
    public void halveReads(int fileIndex) {
	numReads.set(fileIndex, numReads.get(fileIndex)/2);
    }

    /**
     * Returns an array containing the total number of full-length reads
     * in each file.
     */
    public int[] getNumReads() {
	int[] reads = new int[numReads.length()];
	for (int i=0; i<reads.length; i++)
	    reads[i] = numReads.get(i);
	return reads;
    }

    /**
     * Returns an array containing the number of full-length reads
     * in each file mapping to assembled transcripts.
     */
    public int[] getNumMappingReads() {
	int[] mappingReads = new int[numMappingReads.length()];
	for (int i=0; i<mappingReads.length; i++)
	    mappingReads[i] = numMappingReads.get(i);
	return mappingReads;
    }

    /**
     * Returns the total number of full-length reads
     * in the specified file.
     */
    public int getNumReads(int fileIndex) {
	return numReads.get(fileIndex);
    }

    /**
     * Returns the number of full-length reads
     * in the specified file mapping to assembled transcripts.
     */
    public int getNumMappingReads(int fileIndex) {
	return numMappingReads.get(fileIndex);
    }

    /**
     * Returns a 2D array representating the total number of reads mapping to every nt
     * in each replicate.
     */
    public ArrayList<ArrayList<Long>> getTotalReads() {
	ArrayList<ArrayList<Long>> totalReads = new ArrayList<ArrayList<Long>>(Assembler.conditionFiles.size());
	int linearIndex = 0;
	for (int i=0; i<Assembler.conditionFiles.size(); i++) {
	    String[] files = Assembler.conditionFiles.get(i).split(",");
	    totalReads.add(new ArrayList<Long>(files.length));
	    for (int j=0; j<files.length; j++) {
		totalReads.get(i).add((long)0);
		for (int k=0; k<readCounts[linearIndex].length(); k++) {
		    totalReads.get(i).set(j, totalReads.get(i).get(j) + readCounts[linearIndex].get(k));
		}
		linearIndex++;
	    }
	}
	return totalReads;
    }



    /**************************************************
     **********   Private Instance Methods   **********
     **************************************************/

    // Return a row of the cyclic rotations matrix.
    private String getStringInRotationsMatrix(int row) {
	return sequence.substring(rotations[row]) + sequence.substring(0, rotations[row]);
    }

    // Reverse the BWT transform
    private String unpermute() {
	StringBuilder sb = new StringBuilder();
	int r = 0;
	while (BWT.charAt(r) != '$') {
	    sb.append(BWT.charAt(r));
	    r = stepLeft(r);
	}
	return sb.reverse().toString();
    }

    /**
     * Returns number of characters lexographically less than x.
     */
    private int getC(int x) {
	return this.C[x];
    }

    /**
     * Returns number of occurrences of character x up to row y.
     */
    private int getOcc(int x, int y) {
	return this.occ[x][y];
    }

    /**
     * Returns index in Burrows-Wheeler rotations matrix at row x.
     */
    private int getRotations(int x) {
	return this.rotations[x];
    }

    /**
     * Compute the Burrows-Wheeler matrix of cyclic rotations.
     * Rather than return a matrix, we return an array of integers
     * where each value in the array is an index into the sequence.
     * So a "row" of the matrix corresponds to the cyclic sequence
     * beginnning at the specified index. All rows of the matrix
     * have length equal to the length of the sequence.
     */
    private int[] rotations() {
	int[] indices = new int[this.length];
	for (int i=0; i<length; i++) indices[i] = i;
	quicksort(indices, 0, length-1);
	return indices;
    }

    /**
     * Compute the Burrows-Wheeler transform from the "rotations"
     * matrix, i.e., the last column of the matrix.
     */
    private String getTransform() {
	StringBuilder sb = new StringBuilder(length);
	for (int i=0; i<length; i++) {
	    sb.append(sequence.charAt((rotations[i]-1+length) % length));;
	}
	return sb.toString();
    }

    /**
     * Precompute information about how often each character occurs in BWT
     */
    private void precomputeCharacterInfo() {

	// Populate instance variable "C"
	this.C = new int[6];  // All characters in DNA sequence have
	                      // ASCII value less than 123.
	for (int i=0; i<length; i++) C[Dictionary.charToInt(BWT.charAt(i))] += 1;

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
	    int index = Dictionary.charToInt(BWT.charAt(i));
	    if ((index >= 0) && (index < 5)) occ[index][i+1]++;
	}
    }

    /**
     * Helper method for Burrows-Wheeler transform.
     * Taken from appendix of Bowtie article.
     */
    private int stepLeft(int r) {
	int index = Dictionary.charToInt(BWT.charAt(r));
	return C[index] + occ[index][r];
    }

    /**
     * We quicksort the array. But we don't sort the values at
     * each array index, since each value is an index 
     * corresponding to a String. We sort the Strings specified
     * by each array index.
     */
    private void quicksort(int[] a, int p, int r) {
	if (p < r) {
	    int q = partition(a, p, r);
	    quicksort(a, p, q-1);
	    quicksort(a, q+1, r);
	}
    }

    /**
     * Quicksort helper method
     */
    private int partition(int[] a, int p, int r) {
	int x = rand.nextInt(r - p + 1) + p;
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

    /**
     * Quicksort partition helper method
     */
    private void swap(int a[], int i, int j) {
	int temp = a[i];
	a[i] = a[j];
	a[j] = temp;
    }

    /**
     * Given two indices in "sequence", determine whether the
     * cyclic String corresponding to the first index is less 
     * than or equal to the cyclic String corresponding to
     * the second index.
     */
    private boolean lessThanOrEqualTo(int i, int j) {
	for (int z=0; z<length; z++) {
	    if (sequence.charAt((i+z)%length) < sequence.charAt((j+z)%length))
		return true;
	    if (sequence.charAt((i+z)%length) > sequence.charAt((j+z)%length))
		return false;

	}
	return true;  // We have equality of the two Strings
    }



    /**********************************************
     **********   Public Class Methods   **********
     **********************************************/



    /***********************************************
     **********   Private Class Methods   **********
     ***********************************************/

    /**
     * Combines two lists of exact matches, one for each read in a paired-end read,
     * into parallel lists of coordinates of paired-end matches. A combined match occurs when
     * there are matches from each of the two input lists that are
     * within the specified distance from each other.
     */
    private static void combinePairedEndHits(int length1, int length2, ArrayList<Integer> hits1, ArrayList<Integer> hits2, ArrayList<Integer> starts, ArrayList<Integer> ends) {
	for (Integer start1 : hits1) {
	    for (int z=0; z<hits2.size(); z++) {
		Integer start2 = hits2.get(z);
		if ((start1.compareTo(start2) <= 0) && (start2+length2-start1 <= Assembler.maxPairedEndLength)) {
		    starts.add(start1);
		    ends.add(start2+length2);
		    hits2.remove(z);
		    break;
		}
		if ((start2.compareTo(start1) <= 0) && (start1+length1-start2 <= Assembler.maxPairedEndLength)) {
		    starts.add(start2);
		    ends.add(start1+length1);
		    hits2.remove(z);
		    break;
		}
	    }
	}
    }



    /*************************************
     **********   Main Method   **********
     *************************************/

    public static void main(String[] args) {

	// Test DeNovoIndex class
	StringBuilder sb1 = new StringBuilder("AAAAACCCCCGGGGGTTTTT");
	StringBuilder sb2 = new StringBuilder("AAAAAAAAAAAAAAAAAAAA");
	StringBuilder sb3 = new StringBuilder("CCCCCCCCCCCCCCCCCCCC");
	StringBuilder sb4 = new StringBuilder("GTGTGTGTGTGTGTGTGTGT");
	ArrayList<StringBuilder> v = new ArrayList<StringBuilder>();
	v.add(sb1);
	v.add(sb2);
	v.add(sb3);
	v.add(sb4);
	Dictionary d = new Dictionary(25);  // Needed to populate alphabet_2
	DeNovoIndex index = new DeNovoIndex(v);
	System.out.println(index.sequence);
	System.out.println("Length:\t" + index.length);
	System.out.println();
	System.out.println(index.BWT);
	System.out.println();
	System.out.println(index.unpermute());
	System.out.println();

	String s1 = "AAAAA";
	System.out.println("Is exact match:\t" + s1 + "\t" + index.exactMatch(s1));
	String s2 = "TGTGT";
	System.out.println("Is exact match:\t" + s2 + "\t" + index.exactMatch(s2));
	String s3 = "CCGTG";
	System.out.println("Is exact match:\t" + s3 + "\t" + index.exactMatch(s3));
	String s4 = "ACAC";
	System.out.println("Is exact match:\t" + s4 + "\t" + index.exactMatch(s4));
	String s5 = "CCCCCCCCCCCCCCCCCCCC";
	System.out.println("Is exact match:\t" + s5 + "\t" + index.exactMatch(s5));
	String s6 = "CCCCCCCCCCCCCCCCCCCCC";
	System.out.println("Is exact match:\t" + s6 + "\t" + index.exactMatch(s6));
	String s7 = "AAAAAAAAAACCCCCCCCCC";
	System.out.println("Is exact match:\t" + s7 + "\t" + index.exactMatch(s7));
	String s8 = "GGGGGAAAAA";
	System.out.println("Is exact match:\t" + s8 + "\t" + index.exactMatch(s8));
    }

}

