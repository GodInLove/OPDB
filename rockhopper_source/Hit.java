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

/**
 * A Hit represents an alignment between a read and a genome.
 */
public class Hit implements Comparable<Hit> {

    private int start;    // Start coordinate of alignment in genome
    private int stop;     // Stop coordinate of alignment in genome
    private char strand;  // Strand of alignment in genome
    private int errors;   // Number of errors/mismatches in alignment
    private int repliconIndex;  // The index of the replicon to which the read was mapped

    public Hit(int start, int stop, char strand, int errors, int repliconIndex) {
	this.start = start;
	this.stop = stop;
	this.strand = strand;
	this.errors = errors;
	this.repliconIndex = repliconIndex;
    }

    public String getFormattedHit(String sequence, String read) {
	StringBuilder sb = new StringBuilder(read.length()*3 + 40);
	sb.append("\t" + (start+1) + "\t");
	String subSequence = sequence.substring(start-1, stop);
	if (strand == '-') subSequence = Index.reverseComplement(subSequence);
	sb.append(subSequence);
	sb.append("\t" + (stop+1) + "\t" + "(" + strand + ")" + "\n" + "\t" + "\t");
	for (int i=0; i<read.length(); i++) {
	    if (subSequence.charAt(i) == read.charAt(i)) sb.append('|');
	    else sb.append(' ');
	}
	sb.append("\t" + errors + "\n" + "\t" + "\t");
	sb.append(read);
	sb.append("\n");
	return sb.toString();
    }

    public String toString() {
	return start + "\t" + stop + "\t" + strand + "\t" + errors;
    }

    public int getStart() {
	return start;
    }

    public int getStop() {
	return stop;
    }

    public char getStrand() {
	return strand;
    }

    public int getErrors() {
	return errors;
    }

    public int getRepliconIndex() {
	return this.repliconIndex;
    }

    public int getLength() {
	return stop - start + 1;
    }

    public int compareTo(Hit h) {
	if (this.repliconIndex < h.repliconIndex) return -1;
	else if (this.repliconIndex > h.repliconIndex) return 1;
	else {  // Same replicon
	    if (this.start < h.start) return -1;
	    else if (this.start > h.start) return 1;
	    else {  // start coordinates are equal
		if (this.stop < h.stop) return -1;
		else if (this.stop > h.stop) return -1;
		else {  // stop coordinates are equal
		    if ((this.strand == '+') && (h.strand == '-')) return -1;
		    else if ((this.strand == '-') && (h.strand == '+')) return 1;
		    else return 0;
		}
	    }
	}
    }



    /**********************************************
     **********   Public Class Methods   **********
     **********************************************/    

    /**
     * Combines two lists of Hits, one for each read in a paired-end read,
     * into a list of paired-end Hits. A combined Hit occurs when
     * there are Hits from each of the two input lists that are
     * on the same strand and within the specified distance from each other.
     */
    public static void combinePairedEndHits(ArrayList<Hit> combinedHits, ArrayList<Hit> hits1, ArrayList<Hit> hits2, int maxPairedEndLength) {
	combinedHits.clear();
	for (Hit h1 : hits1) {
	    for (Hit h2 : hits2) {
		if ((h1.strand == h2.strand) && (h1.repliconIndex == h2.repliconIndex)) {
		    if (Math.max(h1.stop, h2.stop) - Math.min(h1.start, h2.start) <= maxPairedEndLength) {
			combinedHits.add(new Hit(Math.min(h1.start, h2.start), Math.max(h1.stop, h2.stop), h1.strand, Math.max(h1.errors, h2.errors), h1.repliconIndex));
			break;
		    }
		}
	    }
	}
    }

}

