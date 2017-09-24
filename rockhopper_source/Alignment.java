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

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/*************************************************************
 * An instance of the Alignment class represents an alignment
 * of a sequence read to a reference (genome) sequence.
 *************************************************************/
public class Alignment {

    /*********************************************************
     ****************** INSTANCE VARIABLES *******************
     *********************************************************/

    // Parameters
    private int numErrors;
    private int threshold;

    private String seq1;
    private int start1;
    private int stop1;
    private String seq2;
    private int start2;
    private int stop2;
    private int[] qualities;
    private int[][] table;
    private int[][] errors;
    private int[][] backtrack;  // 1=LEFT, 2=DIAGONAL, 3=ABOVE, -1=DONE
    private int size;
    private int optimalScore_post;
    private int optimalRow_post;
    private int optimalCol_post;
    private int optimalError_post;
    private int optimalScore_pre;
    private int optimalRow_pre;
    private int optimalCol_pre;
    private int optimalError_pre;
    private StringBuilder formattedAlignment;



    /*********************************************************
     ****************** CONSTRUCTORS *************************
     *********************************************************/

    public Alignment() {
	this(0, 0);
    }

    public Alignment(int numErrors, int threshold) {
	this.size = 120;
	this.table = new int[size][size];
	this.errors = new int[size][size];
	this.numErrors = numErrors;
	this.threshold = threshold;
    }



    /*********************************************************
     ****************** PUBLIC INSTANCE METHODS **************
     *********************************************************/

    public void setNumErrors(int numErrors) {
	this.numErrors = numErrors;
    }

    public void setThreshold(int threshold) {
	this.threshold = threshold;
    }

    public void align(File file1, File file2) {
	String genome = readSequenceFromFile(file1);
	String read = readSequenceFromFile(file2);
	int[] qualityScores = new int[read.length()];
	for (int i=0; i<qualityScores.length; i++) qualityScores[i] = 40;
	this.align(genome, -1, 0, read, -1, 0, qualityScores);
    }

    /**
     * Align seq1 (genome) and seq2 (read). The exclusive (start1:stop1) and (start2:stop2)
     * coordinates represent a perfect match seed between the two sequences.
     */
    public void align(String seq1, int start1, int stop1, String seq2, int start2, int stop2, int[] qualities) {
	optimalError_post = 0;
	optimalError_pre = 0;

	if (seq1.length()-stop1 < seq2.length()-stop2) {  // genome is not long enough for read
	    optimalScore_post = -2;
	    optimalScore_pre = -2;
	    return;
	}
	if (start1 < start2) {  // genome is not long enough for read
	    optimalScore_post = -2;
	    optimalScore_pre = -2;
	    return;
	}

	this.seq1 = seq1;
	this.start1 = start1;
	this.stop1 = stop1;
	this.seq2 = seq2;
	this.start2 = start2;
	this.stop2 = stop2;
	this.qualities = qualities;
	if (this.seq2.length() > this.size - this.numErrors) {
	    this.size = this.seq2.length() + this.numErrors;
	    this.table = new int[this.size][this.size];
	    this.errors = new int[this.size][this.size];
	}

	// Align post-seed
	this.optimalScore_post = Integer.MAX_VALUE;
	this.optimalRow_post = -2;
	this.optimalCol_post = -2;
	this.optimalError_post = 0;
	postAlign();

	// Align pre-seed
	this.optimalScore_pre = Integer.MAX_VALUE;
	this.optimalRow_pre = -2;
	this.optimalCol_pre = -2;
	this.optimalError_pre = 0;
	preAlign();
    }

    public int getScore() {
	if ((optimalScore_pre == -2) || (optimalScore_post == -2)) return -2;
	return optimalScore_pre + optimalScore_post;
    }

    public int getErrors() {
	return optimalError_pre + optimalError_post;
    }

    /**
     * Return the starting coordinate in the genome (seq1) of the alignment.
     */
    public int getStart() {
	if (optimalRow_pre >= 0)  // There is a region preceding the seed
	    return start1 - optimalRow_pre + 1;
	else  // No region in front of the seed
	    return start1 + 2;
    }

    /**
     * Return the stopping coordinate in the genome (seq1) of the alignment.
     */
    public int getStop() {
	if (optimalRow_post >= 0)  // There is a region following the seed
	    return stop1 + optimalRow_post + 1;
	else  // No region following the seed
	    return stop1;
    }

    /**
     * Returns a String representation of this Alignment.
     */
    public String toString() {
	if ((optimalScore_pre == -2) || (optimalScore_post == -2)) return "No Alignment.";
	return "\n" + "Errors:\t" + (optimalError_post+optimalError_pre) + "\n" + "Score:\t" + (optimalScore_post+optimalScore_pre);
    }



    /*********************************************************
     ****************** PUBLIC CLASS METHODS *****************
     *********************************************************/



    /*********************************************************
     ****************** PRIVATE INSTANCE METHODS *************
     *********************************************************/

    /**
     * Align region after seed.
     */
    private void postAlign() {

	// Initialize first row and column in table
	if ((stop1 >= seq1.length()) || (stop2 >= seq2.length())) {
	    optimalScore_post = 0;
	    optimalError_post = 0;
	    optimalRow_post = -2;
	    return;
	}
	table[0][0] = mismatch(stop1, stop2)*qualities[stop2];
	errors[0][0] = mismatch(stop1, stop2);
	int j = 0;
	for (j=1; j<Math.min(numErrors+1, seq2.length()-stop2); j++) {  // First row
	    int gapLeft = table[0][j-1] + qualities[stop2+j];
	    int match = j*qualities[stop2+j] + mismatch(stop1, stop2+j)*qualities[stop2+j];
	    if (optimalError_pre + errors[0][j-1]+1 > numErrors) gapLeft = Integer.MAX_VALUE;
	    if (optimalError_pre + j+mismatch(stop1, stop2+j) > numErrors) match = Integer.MAX_VALUE;
	    int minimum = Math.min(gapLeft, match);
	    table[0][j] = minimum;
	    if (match == minimum) {
		errors[0][j] = j + mismatch(stop1, stop2+j);
	    } else {  // Gap left
		errors[0][j] = errors[0][j-1] + 1;
	    }
	}
	for (int i=1; i<Math.min(numErrors+1, seq1.length()-stop1); i++) {  // First column
	    int gapAbove = table[i-1][0] + qualities[stop2];
	    int match = i*qualities[stop2] + mismatch(stop1+i, stop2)*qualities[stop2];
	    if (optimalError_pre + errors[i-1][0]+1 > numErrors) gapAbove = Integer.MAX_VALUE;
	    if (optimalError_pre + i+mismatch(stop1+i, stop2) > numErrors) match = Integer.MAX_VALUE;
	    int minimum = Math.min(gapAbove, match);
	    table[i][0] = minimum;
	    if (match == minimum) {
		errors[i][0] = i + mismatch(stop1+i, stop2);
	    } else {  // Gap above
		errors[i][0] = errors[i-1][0] + 1;
	    }
	}

	// Check if this is optimal alignment (i.e., final column in table)
	if (stop2 == seq2.length()-1) {  // Only one column in table
	    optimalScore_post = table[0][0];
	    optimalRow_post = 0;
	    optimalCol_post = 0;
	    optimalError_post = errors[0][0];
	    for (int i=1; i<Math.min(numErrors+1, seq1.length()-stop1); i++) {
		if (table[i][0] < optimalScore_post) {
		    optimalScore_post = table[i][0];
		    optimalRow_post = i;
		    optimalCol_post = 0;
		    optimalError_post = errors[i][0];
		}
	    }
	} else if (stop2+j-1 == seq2.length()-1) {
	    optimalScore_post = table[0][j-1];
	    optimalRow_post = 0;
	    optimalCol_post = j-1;
	    optimalError_post = errors[0][j-1];
	}

	// Populate table
	for (int i=1; i<Math.min(seq2.length()-stop2+numErrors, seq1.length()-stop1); i++) {
	    for (j=Math.max(i-numErrors, 1); j<Math.min(i+numErrors+1, seq2.length()-stop2); j++) {

		int gapLeft = Integer.MAX_VALUE;
		int gapAbove = Integer.MAX_VALUE;
		int gapLeft_errors = Integer.MAX_VALUE;
		int gapAbove_errors = Integer.MAX_VALUE;
		if (j == i-numErrors) {  // Left table entry is unavailable
		    gapAbove = table[i-1][j] + qualities[stop2+j];
		    gapAbove_errors = errors[i-1][j] + 1;
		} else if (j == i+numErrors) {  // Above table entry is unavailable
		    gapLeft = table[i][j-1] + qualities[stop2+j];
		    gapLeft_errors = errors[i][j-1] + 1;
		} else {  // Both left and above table entries are available
		    gapAbove = table[i-1][j] + qualities[stop2+j];
		    gapLeft = table[i][j-1] + qualities[stop2+j];
		    gapAbove_errors = errors[i-1][j] + 1;
		    gapLeft_errors = errors[i][j-1] + 1;
		}
		int match = table[i-1][j-1] + mismatch(stop1+i, stop2+j)*qualities[stop2+j];
		int match_errors = errors[i-1][j-1] + mismatch(stop1+i, stop2+j);
		if ((match_errors < gapAbove_errors) || (gapLeft_errors < gapAbove_errors) || (gapAbove_errors + optimalError_pre > numErrors)) gapAbove = Integer.MAX_VALUE;
		if ((match_errors < gapLeft_errors) || (gapAbove_errors < gapLeft_errors) || (gapLeft_errors + optimalError_pre > numErrors)) gapLeft = Integer.MAX_VALUE;
		if ((gapAbove_errors < match_errors) || (gapLeft_errors < match_errors) || (match_errors + optimalError_pre > numErrors)) match = Integer.MAX_VALUE;
		int minimum = min(match, gapAbove, gapLeft);

		// Assign value to table
		table[i][j] = minimum;
		if (match == minimum) {
		    errors[i][j] = errors[i-1][j-1] + mismatch(stop1+i, stop2+j);
		} else if (gapAbove == minimum) {
		    errors[i][j] = errors[i-1][j] + 1;
		} else if (gapLeft == minimum) {
		    errors[i][j] = errors[i][j-1] + 1;
		} else {
		    // Impossible case. Do nothing.
		}

		// Check if this is optimal alignment (i.e., final column in table)
		if ((stop2 + j == seq2.length()-1) && (table[i][j] < optimalScore_post)) {
		    optimalScore_post = table[i][j];
		    optimalRow_post = i;
		    optimalCol_post = j;
		    optimalError_post = errors[i][j];
		}
	    }
	}
    }

    /**
     * Align region before seed.
     */
    private void preAlign() {

	// Initialize first row and column in table
	if ((start1 < 0) || (start2 < 0)) {
	    optimalScore_pre = 0;
	    optimalError_pre = 0;
	    optimalRow_pre = -2;
	    return;
	}
	table[0][0] = mismatch(start1, start2)*qualities[start2];
	errors[0][0] = mismatch(start1, start2);
	int j = 0;
	for (j=1; j<Math.min(numErrors+1, start2+1); j++) {  // First row
	    int gapLeft = table[0][j-1] + qualities[start2-j];
	    int match = j*qualities[start2-j] + mismatch(start1, start2-j)*qualities[start2-j];
	    if (optimalError_post + errors[0][j-1]+1 > numErrors) gapLeft = Integer.MAX_VALUE;
	    if (optimalError_post + j+mismatch(start1, start2-j) > numErrors) match = Integer.MAX_VALUE;
	    int minimum = Math.min(gapLeft, match);
	    table[0][j] = minimum;
	    if (match == minimum) {
		errors[0][j] = j + mismatch(start1, start2-j);
	    } else {  // Gap left
		errors[0][j] = errors[0][j-1] + 1;
	    }
	}
	for (int i=1; i<Math.min(numErrors+1, start1+1); i++) {  // First column
	    int gapAbove = table[i-1][0] + qualities[start2];
	    int match = i*qualities[start2] + mismatch(start1-i, start2)*qualities[start2];
	    if (optimalError_post + errors[i-1][0]+1 > numErrors) gapAbove = Integer.MAX_VALUE;
	    if (optimalError_post + i+mismatch(start1-i, start2) > numErrors) match = Integer.MAX_VALUE;
	    int minimum = Math.min(gapAbove, match);
	    table[i][0] = minimum;
	    if (match == minimum) {
		errors[i][0] = i + mismatch(start1-i, start2);
	    } else {  // Gap above
		errors[i][0] = errors[i-1][0] + 1;
	    }
	}

	// Check if this is optimal alignment (i.e., final column in table)
	if (start2 == 0) {  // Only one column in table
	    optimalScore_pre = table[0][0];
	    optimalRow_pre = 0;
	    optimalCol_pre = 0;
	    optimalError_pre = errors[0][0];
	    for (int i=1; i<Math.min(numErrors+1, start1+1); i++) {
		if (table[i][0] < optimalScore_pre) {
		    optimalScore_pre = table[i][0];
		    optimalRow_pre = i;
		    optimalCol_pre = 0;
		    optimalError_pre = errors[i][0];
		}
	    }
	} else if (start2-j+1 == 0) {
	    optimalScore_pre = table[0][j-1];
	    optimalRow_pre = 0;
	    optimalCol_pre = j-1;
	    optimalError_pre = errors[0][j-1];
	}

	// Populate table
	for (int i=1; i<Math.min(start2+1+numErrors, start1+1); i++) {
	    for (j=Math.max(i-numErrors, 1); j<Math.min(i+numErrors+1, start2+1); j++) {

		int gapLeft = Integer.MAX_VALUE;
		int gapAbove = Integer.MAX_VALUE;
		int gapLeft_errors = Integer.MAX_VALUE;
		int gapAbove_errors = Integer.MAX_VALUE;
		if (j == i-numErrors) {  // Left table entry is unavailable
		    gapAbove = table[i-1][j] + qualities[start2-j];
		    gapAbove_errors = errors[i-1][j] + 1;
		} else if (j == i+numErrors) {  // Above table entry is unavailable
		    gapLeft = table[i][j-1] + qualities[start2-j];
		    gapLeft_errors = errors[i][j-1] + 1;
		} else {  // Both left and above table entries are available
		    gapAbove = table[i-1][j] + qualities[start2-j];
		    gapLeft = table[i][j-1] + qualities[start2-j];
		    gapAbove_errors = errors[i-1][j] + 1;
		    gapLeft_errors = errors[i][j-1] + 1;
		}
		int match = table[i-1][j-1] + mismatch(start1-i, start2-j)*qualities[start2-j];
		int match_errors = errors[i-1][j-1] + mismatch(start1-i, start2-j);
		if ((match_errors < gapAbove_errors) || (gapLeft_errors < gapAbove_errors) || (gapAbove_errors + optimalError_post > numErrors)) gapAbove = Integer.MAX_VALUE;
		if ((match_errors < gapLeft_errors) || (gapAbove_errors < gapLeft_errors) || (gapLeft_errors + optimalError_post > numErrors)) gapLeft = Integer.MAX_VALUE;
		if ((gapAbove_errors < match_errors) || (gapLeft_errors < match_errors) || (match_errors + optimalError_post > numErrors)) match = Integer.MAX_VALUE;
		int minimum = min(match, gapAbove, gapLeft);

		// Assign value to table
		table[i][j] = minimum;
		if (match == minimum) {
		    errors[i][j] = errors[i-1][j-1] + mismatch(start1-i, start2-j);
		} else if (gapAbove == minimum) {
		    errors[i][j] = errors[i-1][j] + 1;
		} else if (gapLeft == minimum) {
		    errors[i][j] = errors[i][j-1] + 1;
		} else {
		    // Impossible case. Do nothing.
		}

		// Check if this is optimal alignment (i.e., final column in table)
		if ((start2 - j == 0) && (table[i][j] < optimalScore_pre)) {
		    optimalScore_pre = table[i][j];
		    optimalRow_pre = i;
		    optimalCol_pre = j;
		    optimalError_pre = errors[i][j];
		}
	    }
	}
    }

    /**
     * Once an alignment score has been computed, the backtracking table is searched to
     * determine a StringBuilder representation of the alignment.
     */
    private void backtrack_post(StringBuilder s1_alignment, StringBuilder alignment, StringBuilder s2_alignment) {

        // Begin backtracking search at optimal alignment score (found in last column in table)
        int i = optimalRow_post;
        int j = optimalCol_post;  // Last column in table

        // Backtrack through table to beginning of alignment
        while ((i >= 0) || (j >= 0)) {
	    if (backtrack[i][j] == 2) {  // Diagonal
		s1_alignment.append(seq1.charAt(stop1+i));
		alignment.append(getMatchMismatchChar(seq1.charAt(stop1+i), seq2.charAt(stop2+j)));
		s2_alignment.append(seq2.charAt(stop2+j));
		i--;
		j--;
	    } else if (backtrack[i][j] == 3) {  // Gap above
		s1_alignment.append(seq1.charAt(stop1+i));
		alignment.append(' ');
		s2_alignment.append('-');
		i--;
	    } else if (backtrack[i][j] == 1) {  // Gap left
		s1_alignment.append('-');
		alignment.append(' ');
		s2_alignment.append(seq2.charAt(stop2+j));
		j--;
	    } else if (backtrack[i][j] == -1) {  // DONE
		s1_alignment.append(seq1.charAt(stop1+i));
		alignment.append(getMatchMismatchChar(seq1.charAt(stop1+i), seq2.charAt(stop2+j)));
		s2_alignment.append(seq2.charAt(stop2+j));
		while (i != 0) {  // We're in first column. Alignment begins with gaps.
		    s1_alignment.append(seq1.charAt(stop1+i-1));
		    alignment.append(' ');
		    s2_alignment.append('-');
		    i--;
		}
		while (j != 0) {  // We're in first row. Alignment begins with gaps.
		    s1_alignment.append('-');
		    alignment.append(' ');
		    s2_alignment.append(seq2.charAt(stop2+j-1));
		    j--;
		}
		i--;
		j--;
	    } else {  // Impossible case.
		Peregrine.output("\nThere was an error when backtracking.\n\n");
		i = -1;
		j = -1;
	    }
        }

	// Add seed region to alignment
	i = stop1 - 1;
	j = stop2 - 1;
	while ((i > start1) && (j > start2)) {
	    s1_alignment.append(seq1.charAt(i));
	    alignment.append(getMatchMismatchChar(seq1.charAt(i), seq2.charAt(j)));
	    s2_alignment.append(seq2.charAt(j));
	    i--;
	    j--;
	}

	// Reverse the alignment
	s1_alignment.reverse();
	alignment.reverse();
	s2_alignment.reverse();
    }

    /**
     * Once an alignment score has been computed, the backtracking table is searched to
     * determine a StringBuilder representation of the alignment.
     */
    private void backtrack_pre(StringBuilder s1_alignment, StringBuilder alignment, StringBuilder s2_alignment) {

        // Begin backtracking search at optimal alignment score (found in last column in table)
        int i = optimalRow_pre;
        int j = optimalCol_pre;  // Last column in table

        // Backtrack through table to beginning of alignment
        while ((i >= 0) || (j >= 0)) {
	    if (backtrack[i][j] == 2) {  // Diagonal
		s1_alignment.append(seq1.charAt(start1-i));
		alignment.append(getMatchMismatchChar(seq1.charAt(start1-i), seq2.charAt(start2-j)));
		s2_alignment.append(seq2.charAt(start2-j));
		i--;
		j--;
	    } else if (backtrack[i][j] == 3) {  // Gap above
		s1_alignment.append(seq1.charAt(start1-i));
		alignment.append(' ');
		s2_alignment.append('-');
		i--;
	    } else if (backtrack[i][j] == 1) {  // Gap left
		s1_alignment.append('-');
		alignment.append(' ');
		s2_alignment.append(seq2.charAt(start2-j));
		j--;
	    } else if (backtrack[i][j] == -1) {  // DONE
		s1_alignment.append(seq1.charAt(start1-i));
		alignment.append(getMatchMismatchChar(seq1.charAt(start1-i), seq2.charAt(start2-j)));
		s2_alignment.append(seq2.charAt(start2-j));
		while (i != 0) {  // We're in first column. Alignment begins with gaps.
		    s1_alignment.append(seq1.charAt(start1-i+1));
		    alignment.append(' ');
		    s2_alignment.append('-');
		    i--;
		}
		while (j != 0) {  // We're in first row. Alignment begins with gaps.
		    s1_alignment.append('-');
		    alignment.append(' ');
		    s2_alignment.append(seq2.charAt(start2-j+1));
		    j--;
		}
		i--;
		j--;
	    } else {  // Impossible case.
		Peregrine.output("\nThere was an error when backtracking.\n\n");
		i = -1;
		j = -1;
	    }
        }
    }

    /**
     * Returns 1 if the characters at the specified indices in
     * the two sequences mismatch.
     * Returns 0 if the characters at the specified indices in
     * the two sequence are the same.
     */
    private int mismatch(int x, int y) {
	if (seq1.charAt(x) == seq2.charAt(y)) return 0;
	return 1;
    }

    /**
     * Returns an alignment character if the two specified characters are the same.
     * Returns a space character otherwise.
     */
    private char getMatchMismatchChar(char a, char b) {
	if (a == b) return '|';
	else return ' ';
    }

    /**
     * Returns a String representation of a 2D integer array.
     */
    private String tableToString_post(int[][] a) {
	StringBuilder sb = new StringBuilder();
	for (int j=0; j<seq2.length()-stop2; j++)
	    sb.append("\t" + seq2.charAt(stop2 + j));
	sb.append("\n");
	for (int i=0; i<seq1.length()-stop1; i++) {
	    sb.append(seq1.charAt(stop1 + i));
	    for (int j=0; j<seq2.length()-stop2; j++) {
		sb.append("\t" + a[i][j]);
	    }
	    sb.append("\n");
	}
	return sb.toString();
    }

    /**
     * Returns a String representation of a 2D integer array.
     */
    private String tableToString_pre(int[][] a) {
	StringBuilder sb = new StringBuilder();
	for (int j=0; j<start2+1; j++)
	    sb.append("\t" + seq2.charAt(start2 - j));
	sb.append("\n");
	for (int i=0; i<start1+1; i++) {
	    sb.append(seq1.charAt(start1 - i));
	    for (int j=0; j<start2+1; j++) {
		sb.append("\t" + a[i][j]);
	    }
	    sb.append("\n");
	}
	return sb.toString();
    }

    /**
     * Returns the minimum of three integers.
     */
    private int min(int a, int b, int c) {
	return Math.min(a, Math.min(b, c));
    }



    /*********************************************************
     ****************** PRIVATE CLASS METHODS ****************
     *********************************************************/

    /**
     * Reads in and returns a genomic sequence from the specified FASTA file.
     *
     * @param   f   a <code>File</code> object referring to a FASTA file containing a genomic sequence
     * @return      the genomic sequence read in from a FASTA file
     */
    private static String readSequenceFromFile(File f) {
	StringBuilder sequence = new StringBuilder();
	try {
	    Scanner reader = new Scanner(f);
	    String header = reader.nextLine();  // Header line of FASTA file
	    if ((header.length() == 0) || (header.charAt(0) != '>')) {
		Peregrine.output("Error - first line of file " + f + " is not in FASTA format.\n");
		return sequence.toString();
	    }
	    while (reader.hasNext()) {  // continue until we reach end of file
		sequence.append(reader.nextLine());
	    }
	    reader.close();
	} catch (FileNotFoundException e) {
	    Peregrine.output("Error - the file " + f + " could not be found and opened.\n");
	    return sequence.toString();
	}
	return sequence.toString();
    }



    /*********************************************************
     ****************** MAIN METHOD **************************
     *********************************************************/

    public static void main(String [] args) {

	if (args.length < 2) {
	    System.err.println("\nWhen executing this program, please enter the name of two files,");
	    System.err.println("each containing a sequence. The program will align the two sequences.\n");
	    System.err.println("\tjava Alignment file1.txt file2.txt\n");
	    System.exit(0);
	}

	Alignment a = new Alignment(2, 100);
	a.align(new File(args[0]), new File(args[1]));
	System.out.println(a);
    }

}

