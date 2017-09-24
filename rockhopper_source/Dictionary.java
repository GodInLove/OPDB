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

/**
 * A Dictionary instance represents a dictionary of k-mers 
 * found in a set of sequencing reads.
 */

import java.util.HashMap;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Arrays;

public class Dictionary {

    /********************************************
     **********   Instance Variables   **********
     ********************************************/

    public Table dict;
    private int size;
    private ArrayList<StringBuilder> transcripts = new ArrayList<StringBuilder>();
    private ArrayList<KMer> seeds;
    private int seedIndex;
    public DeNovoIndex bwtIndex;



    /*****************************************
     **********   Class Variables   **********
     *****************************************/

    private static int k;  // Size of k-mer
    private static char[] alphabet_1 = {'A', 'C', 'G', 'T'};
    private static int[] alphabet_2 = new int[128];  // ASCII values



    /**************************************
     **********   Constructors   **********
     **************************************/

    public Dictionary(int k) {
	Dictionary.k = k;
	dict = new Table();
	alphabet_2['A'] = 0;
	alphabet_2['C'] = 1;
	alphabet_2['G'] = 2;
	alphabet_2['T'] = 3;
	alphabet_2['^'] = 4;
	alphabet_2['$'] = 5;
	bwtIndex = new DeNovoIndex(new ArrayList<StringBuilder>());  // Empty transcript index initially
    }



    /*************************************************
     **********   Public Instance Methods   **********
     *************************************************/

    public void add(String read) {
	String read_RC = "";
	if (Assembler.unstranded) read_RC = Assembler.reverseComplement(read);
	for (int i=0; i<read.length()-k+1; i++) {
	    if (!bwtIndex.exactMatch(read, i, i+k)) {  // k-mer not in BWT transcript index
		if (Assembler.unstranded) {  // Strand ambiguous
		    if (!bwtIndex.exactMatch(read_RC, read.length()-i-k, read.length()-i)) {  // revComp of k-mer not in index
			Long k_mer = stringToLong(read, i, i+k);
			Long k_mer_RC = stringToLong(read_RC, read.length()-i-k, read.length()-i);
			dict.add(k_mer, k_mer_RC);
		    }
		} else {  // Strand specific
		    Long k_mer = stringToLong(read, i, i+k);
		    dict.add(k_mer);
		}
	    }
	}
	if (dict.exceedsLoadFactor()) buildTranscriptsAndClearTable();  // Table is full
    }

    /**
     * Returns the number of unique elements in the dictionary.
     */
    public int getSize() {
	return size;
    }

    /**
     * Returns the number of transcripts.
     */
    public int getNumTranscripts() {
	int size = 0;
	for (StringBuilder sb : transcripts) {
	    if (sb.length() >= Assembler.minTranscriptLength) size++;
	}
	return size;
    }

    /**
     * Returns the average transcript length.
     */
    public int getAverageTranscriptLength() {
	long length = 0;
	int count = 0;
	if (transcripts.size() == 0) return 0;
	for (StringBuilder sb : transcripts) {
	    if (sb.length() >= Assembler.minTranscriptLength) {
		length += sb.length();
		count++;
	    }
	}
	if (count == 0) return 0;
	else return (int)(length/count);
    }

    public void initializeReadMapping(boolean stopAfterOneHit) {
	bwtIndex.initializeReadMapping(stopAfterOneHit);
    }

    /**
     * Map the full length read to the index of assembled transcripts.
     */
    public void mapFullLengthRead(String read) {
	bwtIndex.exactMatch_fullRead(read);
    }

    /**
     * Map the full length *paired-end* read to the index of assembled transcripts.
     */
    public void mapFullLengthRead(String read, String read2) {
	bwtIndex.exactMatch_fullRead(read, read2);
    }

    public void assembleTranscripts() {
	size = dict.size();  // Compute size

	// Extend existing transcripts
	for (int i=0; i<transcripts.size(); i++) {
	    extendSeedForward(transcripts.get(i));
	    extendSeedBackward(transcripts.get(i));
	}

	// Create new transcripts
	prepareSeeds();
	Long seed = getSeed();
	while (seed != null) {
	    StringBuilder transcript = new StringBuilder(longToString(seed));
	    dict.remove(seed);  // Remove seed from dict
	    extendSeedForward(transcript);
	    extendSeedBackward(transcript);
	    seed = getSeed();
	    transcripts.add(transcript);
	}
	dict = new Table();
	bwtIndex = new DeNovoIndex(transcripts);
	System.gc();
    }



    /**************************************************
     **********   Private Instance Methods   **********
     **************************************************/

    /**
     * Check if the table is full.
     * If so, assemble transcripts via deBruijn graph and
     * clear the dictionary.
     */
    private synchronized void buildTranscriptsAndClearTable() {
	if (dict.exceedsLoadFactor()) {  // Table is full
	    assembleTranscripts();
	}
    }

    /**
     * Returns the most frequently occurring k-mer subject
     * to the following restrictions:
     *  - The seed is sufficiently expressed, i.e., sufficiently many reads supporting the k-mer
     *  - The information content is at least 1.5
     *  - The k-mer is not palindromic (deleted).
     */
    private Long getSeed() {
	while (seedIndex < seeds.size()) {
	    Long key = seeds.get(seedIndex).key;
	    seedIndex++;
	    if (dict.containsKey(key) && (informationContent(key) >= 1.5)) return key;
	}
	return null;
    }

    /**
     * Since we need to repeatedly determine the best seed in the
     * dictionary (an expensive operation for a hashtable), we 
     * create a sorted vector so that finding the best seed is
     * efficient.
     */
    private void prepareSeeds() {
	seeds = new ArrayList<KMer>();
	for (int i=0; i<dict.capacity; i++) {
	    int value = dict.getValueAtIndex(i);
	    if (value == 0) continue;
	    Long key = dict.getKeyAtIndex(i);
	    if (value >= Assembler.minSeedExpression) seeds.add(new KMer(key, value));
	}
	Collections.sort(seeds, Collections.reverseOrder());
	seedIndex = 0;
    }

    /**
     * Extends transcript 5' to 3'.
     */
    private void extendSeedForward(StringBuilder transcript) {
	boolean done = false;
	while (!done) {
	    int max = Assembler.minExpression-1;
	    int maxIndex = -1;
	    transcript.append('?');
	    for (int i=0; i<alphabet_1.length; i++) {
		transcript.setCharAt(transcript.length()-1, alphabet_1[i]);
		Long k_mer = stringToLong(transcript, transcript.length()-k, transcript.length());
		if (dict.containsKey(k_mer) && (dict.get(k_mer) > max)) {
		    max = dict.get(k_mer);
		    maxIndex = i;
		}
		if (Assembler.unstranded) {  // Strand ambiguous
		    Long k_mer_RC = stringToLong(Assembler.reverseComplement(transcript.substring(transcript.length()-k)), 0, k);
		    if (dict.containsKey(k_mer_RC) && (dict.get(k_mer_RC) > max)) {
			max = dict.get(k_mer_RC);
			maxIndex = i;
		    }
		}
	    }
	    if (maxIndex >= 0) {  // Found an extension
		transcript.setCharAt(transcript.length()-1, alphabet_1[maxIndex]);
		Long key = stringToLong(transcript, transcript.length()-k, transcript.length());
		dict.remove(key);  // Remove k-mer from dict
		if (Assembler.unstranded) {  // Strand ambiguous
		    Long key_RC = stringToLong(Assembler.reverseComplement(transcript.substring(transcript.length()-k)), 0, k);
		    dict.remove(key_RC);  // Remove k-mer from dict
		}
	    } else {  // No extension
		transcript.deleteCharAt(transcript.length()-1);
		done = true;
	    }
	}
    }

    /**
     * Extends transcript 3' to 5'.
     */
    private void extendSeedBackward(StringBuilder transcript) {
	transcript.reverse();
	boolean done = false;
	while (!done) {
	    int max = Assembler.minExpression-1;
	    int maxIndex = -1;
	    transcript.append('?');
	    for (int i=0; i<alphabet_1.length; i++) {
		transcript.setCharAt(transcript.length()-1, alphabet_1[i]);
		Long k_mer = stringToLong_Reverse(transcript, transcript.length()-k, transcript.length());
		if (dict.containsKey(k_mer) && (dict.get(k_mer) > max)) {
		    max = dict.get(k_mer);
		    maxIndex = i;
		}
		if (Assembler.unstranded) {  // Strand ambiguous
		    Long k_mer_RC = stringToLong_Reverse(Assembler.reverseComplement(transcript.substring(transcript.length()-k)), 0, k);
		    if (dict.containsKey(k_mer_RC) && (dict.get(k_mer_RC) > max)) {
			max = dict.get(k_mer_RC);
			maxIndex = i;
		    }
		}
	    }
	    if (maxIndex >= 0) {  // Found an extension
		transcript.setCharAt(transcript.length()-1, alphabet_1[maxIndex]);
		Long key = stringToLong_Reverse(transcript, transcript.length()-k, transcript.length());
		dict.remove(key);  // Remove k-mer from dict
		if (Assembler.unstranded) {  // Strand ambiguous
		    Long key_RC = stringToLong_Reverse(Assembler.reverseComplement(transcript.substring(transcript.length()-k)), 0, k);
		    dict.remove(key_RC);  // Remove k-mer from dict
		}
	    } else {  // No extension
		transcript.deleteCharAt(transcript.length()-1);
		done = true;
	    }
	}
	transcript.reverse();
    }



    /***********************************************
     **********   Private Class Methods   **********
     ***********************************************/

    private static Long stringToLong(String s, int index1, int index2) {
	long num = 0;
	for (int i=index1; i<index2; i++) {
	    num = num << 2;
	    num += charToInt(s.charAt(i));
	}
	return new Long(num);
    }

    private static Long stringToLong(StringBuilder sb, int index1, int index2) {
	long num = 0;
	for (int i=index1; i<index2; i++) {
	    num = num << 2;
	    num += charToInt(sb.charAt(i));
	}
	return new Long(num);
    }

    private static Long stringToLong_Reverse(String s, int index1, int index2) {
	long num = 0;
	for (int i=index2-1; i>=index1; i--) {
	    num = num << 2;
	    num += charToInt(s.charAt(i));
	}
	return new Long(num);
    }

    private static Long stringToLong_Reverse(StringBuilder sb, int index1, int index2) {
	long num = 0;
	for (int i=index2-1; i>=index1; i--) {
	    num = num << 2;
	    num += charToInt(sb.charAt(i));
	}
	return new Long(num);
    }

    private static String longToString(Long key) {
	long num = key.longValue();
	StringBuilder sb = new StringBuilder(k+1);
	for (int i=0; i<k; i++) {
	    sb.append(intToChar((int)num & 3));
	    num = num >>> 2;
	}
	return sb.reverse().toString();
    }

    public static int charToInt(char c) {
	return alphabet_2[c];
    }

    private static char intToChar(int num) {
	return alphabet_1[num];
    }

    private static double informationContent(Long k_mer) {
	double[] counts = new double[4];
	long num = k_mer.longValue();
	for (int i=0; i<k; i++) {
	    counts[((int)num & 3)] += 1.0;
	    num = num >>> 2;
	}

	double sum = 0.0;
	for (int i=0; i<counts.length; i++) sum += counts[i];

	double informationContent = 0.0;
	for (int i=0; i<counts.length; i++) {
	    if (counts[i] > 0.0) {
		double frequency = counts[i] / sum;
		informationContent += frequency * Math.log(frequency) / Math.log(2.0);
	    }
	}
	return 0.0 - informationContent;
    }



    /*************************************
     **********   Main Method   **********
     *************************************/

    public static void main(String[] args) {

	int k = 25;
	Dictionary d = new Dictionary(k);

	String s = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG";
	System.out.println(s);
	for (int i=0; i<d.dict.capacity; i++) {
	    int value = d.dict.getValueAtIndex(i);
	    if (value == 0) continue;
	    Long key = d.dict.getKeyAtIndex(i);
	    System.out.println(longToString(key) + "\t" + value);
	}
    }

}



/**************************************
 **********   Helper Class   **********
 **************************************/

class KMer implements Comparable<KMer> {

    public Long key;
    public Integer count;

    public KMer(Long key, Integer count) {
	this.key = key;
	this.count = count;
    }

    public int compareTo(KMer k) {
	return count.compareTo(k.count);
    }
}
