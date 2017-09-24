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

import java.util.Collections;
import java.util.ArrayList;
import java.util.Random;

public class Index {

    /********************************************
     **********   Instance Variables   **********
     ********************************************/

    private int numReplicons;  // Number of indexed replicons
    private ArrayList<Replicon> replicons;  // Index for each replicon
    private boolean stopAfterOneHit;  // Terminate the search when the first hit is found



    /*****************************************
     **********   Class Variables   **********
     *****************************************/

    public static Random rand = new Random();  // Random number generator



    /**************************************
     **********   Constructors   **********
     **************************************/

    public Index(String[] sequenceFiles) {
	replicons = new ArrayList<Replicon>(sequenceFiles.length);
	for (int i=0; i<sequenceFiles.length; i++)
	    replicons.add(new Replicon(sequenceFiles[i]));
	this.numReplicons = replicons.size();
	this.stopAfterOneHit = true;
    }



    /*************************************************
     **********   Public Instance Methods   **********
     *************************************************/

    /**
     * Set whether the search should terminate after finding the first hit,
     * even if there are other hits with the same score.
     */
    public void setStopAfterOneHit(boolean stopAfterOneHit) {
	this.stopAfterOneHit = stopAfterOneHit;
    }

    /**
     * Return the number of indexed replicons.
     */
    public int getNumReplicons() {
	return this.numReplicons;
    }

    /**
     * Returns the ith Replicon index.
     */
    public Replicon getReplicon(int i) {
	return this.replicons.get(i);
    }

    /**
     * Return all hits of s to the genome sequence.
     */
    public void exactMatch(String s, ArrayList<Hit> hits) {
	hits.clear();
	if (!stopAfterOneHit) {  // Find all optimal hits
	    for (int i=0; i<numReplicons; i++) {
		exactMatch(s, '+', i, 0, hits);                     // Align to plus strand
		exactMatch(reverseComplement(s), '-', i, 0, hits);  // Align to minus strand
	    }
	    Collections.sort(hits);
	    return;
	} else {  // Find only one hit
	    for (int i=0; i<numReplicons; i++) {
		if (rand.nextBoolean()) {  // Randomly align first to plus strand
		    exactMatch(s, '+', i, 0, hits);
		    if (hits.isEmpty()) exactMatch(reverseComplement(s), '-', i, 0, hits);
		    return;
		} else {  // Randomly align first to minus strand
		    exactMatch(reverseComplement(s), '-', i, 0, hits);
		    if (hits.isEmpty()) exactMatch(s, '+', i, 0, hits);
		    return;
		}
	    }
	}
	return;  // This line should never be executed.
    }

    /**
     * Return all hits of s to the genome sequence.
     */
    public void inexactMatch(String s, int[] qualities, Alignment a, ArrayList<Hit> hits, ArrayList<Hit> seedHits) {
	hits.clear();
	seedHits.clear();
	inexactMatch(s, qualities, 0, a, hits, seedHits);
	Collections.sort(hits);
	return;
    }



    /**************************************************
     **********   Private Instance Methods   **********
     **************************************************/

    // Find indices of all occurrences of s in the sequence.
    private void exactMatch(String s, char strand, int repliconIndex, int errors, ArrayList<Hit> hits) {
	exactMatch(s, strand, repliconIndex, errors, hits, this.stopAfterOneHit);
    }

    // Find indices of all occurrences of s in the sequence.
    private void exactMatch(String s, char strand, int repliconIndex, int errors, ArrayList<Hit> hits, boolean stopAfterOneHit_local) {
	Replicon r = replicons.get(repliconIndex);
	int i = s.length() - 1;
	char c = s.charAt(i);
	int sp = r.getC(charToInt(c));
	int ep = -1;
	if (charToInt(c) < 4) ep = r.getC(charToInt(c)+1);
	else ep = r.getLength();
	i--;
	while ((sp < ep) && (i >= 0)) {
	    int index = charToInt(s.charAt(i));
	    sp = r.getC(index) + r.getOcc(index, sp);
	    ep = r.getC(index) + r.getOcc(index, ep);
	    i--;
	}
	if (!stopAfterOneHit_local) {  // Find all optimal hits
	    for (int z=sp; z<ep; z++) hits.add(new Hit(r.getRotations(z)+1, r.getRotations(z)+1+s.length()-1, strand, errors, repliconIndex));  // 1-indexed
	} else {  // Find only one hit
	    if (ep > sp) {
		int randIndex = sp + rand.nextInt(ep-sp);
		hits.add(new Hit(r.getRotations(randIndex)+1, r.getRotations(randIndex)+1+s.length()-1, strand, errors, repliconIndex));  // 1-indexed
	    }
	}
    }

    /**
     * Find all indices in the sequence (with the same 
     * best score) that match the specified query String 
     * with errors that sum to a value at or below the threshold.
     */
    private void inexactMatch(String s, int[] qualities, int threshold, Alignment a, ArrayList<Hit> hits, ArrayList<Hit> seedHits) {

	// Exact match search
	if (!stopAfterOneHit) {  // Find all optimal hits
	    for (int i=0; i<numReplicons; i++) {
		exactMatch(s, '+', i, 0, hits);                     // Align to plus strand
		exactMatch(reverseComplement(s), '-', i, 0, hits);  // Align to minus strand
	    }
	} else {  // Find only one hit
	    for (int i=0; i<numReplicons; i++) {
		if (rand.nextBoolean()) {  // Randomly align first to plus strand
		    exactMatch(s, '+', i, 0, hits);
		    if (hits.isEmpty()) exactMatch(reverseComplement(s), '-', i, 0, hits);
		} else {  // Randomly align first to minus strand
		    exactMatch(reverseComplement(s), '-', i, 0, hits);
		    if (hits.isEmpty()) exactMatch(s, '+', i, 0, hits);
		}
		if (!hits.isEmpty()) return;
	    }
	}
	if (!hits.isEmpty()) return;

	// Inexact search via alignment
	int qualityThreshold = (int)(40 * s.length() * Peregrine.percentMismatches);
	int mismatches = (int)(s.length() * Peregrine.percentMismatches);
	int minSeedLength = Math.max((int)(s.length() * Peregrine.percentSeedLength), Peregrine.minReadLength);
	a.setNumErrors(mismatches);
	a.setThreshold(qualityThreshold);
	if ((mismatches == 0) || (s.length() < minSeedLength)) return;
	int minScore = Integer.MAX_VALUE;
	int numSeeds = Math.min(mismatches+1, (int)Math.ceil((s.length()+1)/(double)minSeedLength));
	int increment = s.length()/numSeeds;
	int seedLength = Math.max(increment, minSeedLength);
	int[] reverseQualities = new int[qualities.length];
	for (int q=0; q<qualities.length; q++) reverseQualities[qualities.length-1-q] = qualities[q];
	for (int i=0; i<numReplicons; i++) {  // Find seed in each replicon

	    // Determine all seeds from read
	    int x = 0;
	    int count = 0;
	    while (count < numSeeds) {  // Try each seed from read
		int start = x;
		int stop = x+seedLength-1;
		if (count == numSeeds-1) {  // Final read
		    start = s.length()-seedLength;
		    stop = s.length()-1;
		}
		x += increment;
		count++;

		seedHits.clear();
		exactMatch(s.substring(start, stop+1), '+', i, 0, seedHits, false);
		exactMatch(reverseComplement(s.substring(start, stop+1)), '-', i, 0, seedHits, false);  // Align to minus strand

		// Align read to replicon sequence based on seedHits
		for (int j=0; j<seedHits.size(); j++) {
		    if (seedHits.get(j).getStrand() == '+')  // Plus strand
			a.align(replicons.get(i).getSequence(), seedHits.get(j).getStart()-2, seedHits.get(j).getStop(), s, start-1, stop+1, qualities);
		    else  // Minus strand
			a.align(replicons.get(i).getSequence(), seedHits.get(j).getStart()-2, seedHits.get(j).getStop(), reverseComplement(s), s.length()-1-stop-1, s.length()-1-start+1, reverseQualities);
		    // Add alignment (if found and if among *best* scoring hits found so far) to "hits"
		    if ((a.getScore() >= 0) && (a.getScore() <= qualityThreshold) && (a.getErrors() <= mismatches)) {  // Valid alignment
			if (a.getScore() < minScore) {  // This is the new best scoring alignment
			    hits.clear();
			    hits.add(new Hit(a.getStart(), a.getStop(), seedHits.get(j).getStrand(), a.getErrors(), i));
			    minScore = a.getScore();
			} else if (a.getScore() == minScore) {  // This is another best scoring alignment
			    hits.add(new Hit(a.getStart(), a.getStop(), seedHits.get(j).getStrand(), a.getErrors(), i));
			} else {  // We have previously found a better scoring alignment
			    // Do nothing.
			}
		    }
		}
	    }

	    // Remove neighboring Hits
	    Collections.sort(hits);
	    int j=1;
	    while (j < hits.size()) {
		if (hits.get(j).getStart() - hits.get(j-1).getStop() < s.length()) hits.remove(j);
		else j++;
	    }
	}

	// Randomly choose one best hit to return (as opposed to all best hits)
	if ((hits.size() > 1) && stopAfterOneHit) {  // Find only one hit (not all optimal hits)
	    Hit randomBest = hits.get(rand.nextInt(hits.size()));
	    hits.clear();
	    hits.add(randomBest);
	}
    }



    /**********************************************
     **********   Public Class Methods   **********
     **********************************************/

    /**
     * Set seed for pseudorandom number generator.
     */
    public static void setRandomSeed(long seed) {
	rand = new Random(seed);
    }

    // Convert a nucleotide character to an integer
    public static int charToInt(char c) {
	if (c == 'A') return 0;
	if (c == 'C') return 1;
	if (c == 'G') return 2;
	if (c == 'T') return 3;
	if (c == '^') return 4;
	if (c == '$') return 5;
	return -1;
    }

    /**
     * Return the reverse complement of the input sequence.
     */
    public static String reverseComplement(String s) {
	StringBuilder sb = new StringBuilder(s.length());
	for (int i=s.length()-1; i>=0; i--) {
	    if (s.charAt(i) == 'A') sb.append('T');
	    else if (s.charAt(i) == 'C') sb.append('G');
	    else if (s.charAt(i) == 'G') sb.append('C');
	    else if (s.charAt(i) == 'T') sb.append('A');
	    else if (s.charAt(i) == '^') sb.append('^');
	    else Peregrine.output("Invalid DNA nucleotide character: " + s.charAt(i) + "\n");
	}
	return sb.toString();
    }



    /***********************************************
     **********   Private Class Methods   **********
     ***********************************************/



    /*************************************
     **********   Main Method   **********
     *************************************/

    public static void main(String[] args) {

        if (args.length < 1) {
            System.err.println("\nUSAGE: java Index <index1.fna,index2.fna,index3.fna>" + "\n");
            System.err.println("Index takes a comma-separated set of sequence files (such as FASTA genome files) to be indexed and it computes the Burrows-Wheeler transformed index for each sequence. If run from the command line, it simply outputs the size of each index to STDOUT. The Index class is meant to be instantiated from another application. \n");
            System.exit(0);
        }

	Index index = new Index(args[0].split(","));
	for (int i=0; i<index.getNumReplicons(); i++)
	    System.out.println(index.getReplicon(i).getName() + " has size " + index.getReplicon(i).getLength());
    }

}

