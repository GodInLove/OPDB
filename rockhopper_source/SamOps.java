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
import java.util.HashMap;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.FileInputStream;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.zip.GZIPInputStream;
import java.lang.reflect.Field;

/**
 * An instance of the SamOps class supports operations
 * on SAM/BAM files. It is assumed that a SAM/BAM file contains
 * information about a single RNA-seq experiment.
 * Java 7 is required for properly handling BAM files.
 */
public class SamOps {

    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private String fileName;
    private boolean isValidSam;
    private boolean isValidBam;
    private boolean isPairedEnd;
    private boolean isDenovo;
    private HashMap<String, Integer> replicons;
    private ArrayList<Integer> repliconLengths;
    public AtomicIntegerArray alignedReads;
    public AtomicInteger unalignedReads;
    public AtomicInteger badReads;
    public AtomicLong avgLengthReads;
    public AtomicIntegerArray[] coordinates_plus;
    public AtomicIntegerArray[] coordinates_minus;
    private int maxQuality;
    private long pairedEndCount;  // Temp counter used when reading in first lines
    private String qualities;  // Temp String used when reading in unaligned reads
    private int[] cigarMap;  // Indexed by ASCII char values {M,I,D,N,S,H,P,=,X}
    private char[] cigarChars = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};
    private static byte[] intBuffer = new byte[4];  // Temp array used for reading bam files
    private static byte[] buffer1 = new byte[1];    // Temp array used for reading bam files
    private static char[] NT_code = {'=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};
    public static final String[] emptyStringArray = new String[0];
    private GZIPInputStream in;  // Used for reading BAM files
    public static int maxPairedEndLength = 500;



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    /**
     * Constructor for reading in a subset of lines in
     * SAM/BAM file to check if file is valid.
     */
    public SamOps(String fileName, long numLines, boolean isDenovo) {
	this.fileName = fileName;
	this.replicons = new HashMap<String, Integer>();
	this.repliconLengths = new ArrayList<Integer>();
	this.alignedReads = new AtomicIntegerArray(1);
	this.unalignedReads = new AtomicInteger(0);
	this.badReads = new AtomicInteger(0);
	this.avgLengthReads = new AtomicLong(0);
	this.isPairedEnd = false;
	this.isDenovo = isDenovo;
	this.pairedEndCount = 0;
	this.qualities = "";
	if (isTextFile(fileName)) this.isValidSam = prepSAM(numLines);
	else this.isValidBam = prepBAM(numLines);
	if (isValidSam || isValidBam) {
	    if ((replicons.size() > 0) && (!isDenovo)) {
		alignedReads = new AtomicIntegerArray(replicons.size());
		coordinates_plus = new AtomicIntegerArray[replicons.size()];
		coordinates_minus = new AtomicIntegerArray[replicons.size()];
		for (int i=0; i<repliconLengths.size(); i++) {
		    coordinates_plus[i] = new AtomicIntegerArray(repliconLengths.get(i)+1);
		    coordinates_minus[i] = new AtomicIntegerArray(repliconLengths.get(i)+1);
		}
	    } else {
		coordinates_plus = new AtomicIntegerArray[1];
		coordinates_minus = new AtomicIntegerArray[1];
		coordinates_plus[0] = new AtomicIntegerArray(0);
		coordinates_minus[0] = new AtomicIntegerArray(0);
	    }
	}
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    /**
     * Returns true if the first lines of the input SAM
     * file are valid (correspond to what we are expecting
     * from a SAM file) and false otherwise.
     */
    public boolean isValidSam() {
	return this.isValidSam;
    }

    /**
     * Returns true if the first lines of the input BAM
     * file are valid (correspond to what we are expecting
     * from a BAM file) and false otherwise.
     */
    public boolean isValidBam() {
	return this.isValidBam;
    }

    /**
     * Returns true if file contains paired end reads.
     * False otherwise.
     * Only defined if file is valid.
     */
    public boolean isPairedEnd() {
	return isPairedEnd;
    }

    /**
     * Returns the maximum quality score found in the first lines of the file.
     */
    public int getMaxQuality() {
	return maxQuality;
    }

    /**
     * Returns the length of the replicon at the specified index.
     */
    public int getRepliconLength(int i) {
	return repliconLengths.get(i);
    }

    /**
     * Returns an array of replicon names. The array is ordered
     * so that the indices in the array are the same as those
     * for the repliconLengths list.
     */
    public String[] getReplicons() {
	String[] names = new String[replicons.size()];
	int[] indices = new int[replicons.size()];
	int i = 0;
	for (String key : replicons.keySet()) {
	    names[i] = key;
	    indices[i] = replicons.get(key);
	    i++;
	}
	mergesort(indices, names, 0, replicons.size()-1);
	return names;
    }

    /**
     * Return a String containing stats/information about this SamOps object.
     */
    public String toString() {
	if (!isValidSam() && !isValidBam()) return "Not a valid SAM/BAM file.\n";
	StringBuilder sb = new StringBuilder();
	sb.append("\n");
	sb.append("Sam/Bam file contains reads mapping to " + replicons.size() + " replicons." + "\n");
	if (replicons.size() == 0) {  // No replicons in SAM/BAM header
	    int final_NT = 0;
	    int repliconIndex = 0;
	    int count_plus = 0;
	    for (int i=0; i<coordinates_plus[repliconIndex].length(); i++) {
		if (coordinates_plus[repliconIndex].get(i) > 0) {
		    count_plus++;
		    if (i > final_NT) final_NT = i;
		}
	    }
	    int count_minus = 0;
	    for (int i=0; i<coordinates_minus[repliconIndex].length(); i++) {
		if (coordinates_minus[repliconIndex].get(i) > 0) {
		    count_minus++;
		    if (i > final_NT) final_NT = i;
		}
	    }
	    sb.append("\t\tPercent of plus strand covered: \t" + (int)(100.0*count_plus/final_NT) + "%" + "\n");
	    sb.append("\t\tPercent of minus strand covered:\t" + (int)(100.0*count_minus/final_NT) + "%" + "\n");
	} else {  // We have replicons in SAM/BAM header
	    for (String name : replicons.keySet()) {
		int repliconIndex = replicons.get(name);
		sb.append("\tReplicon " + name + " with length " + repliconLengths.get(repliconIndex) + "\n");
		int count = 0;
		for (int i=0; i<coordinates_plus[repliconIndex].length(); i++) {
		    if (coordinates_plus[repliconIndex].get(i) > 0) count++;
		}
		sb.append("\t\tPercent of plus strand covered: \t" + (int)(100.0*count/coordinates_plus[repliconIndex].length()) + "%" + "\n");
		count = 0;
		for (int i=0; i<coordinates_minus[repliconIndex].length(); i++) {
		    if (coordinates_minus[repliconIndex].get(i) > 0) count++;
		}
		sb.append("\t\tPercent of minus strand covered:\t" + (int)(100.0*count/coordinates_minus[repliconIndex].length()) + "%" + "\n");
	    }
	}
	sb.append("\n");
	sb.append("Bad reads:      \t\t" + badReads.get() + "\n");
	sb.append("Unaligned reads:\t\t" + unalignedReads.get() + "\n");
	String[] orderedReplicons = getReplicons();
	for (int i=0; i<orderedReplicons.length; i++)
	    sb.append("Aligned reads in " + orderedReplicons[i] + ":  \t" + alignedReads.get(i) + "\n");
	int count = 0;
	return sb.toString();
    }

    /**
     * Parse one line of alignment information from the SAM file.
     * Returns empty String array if line should be skipped.
     * Returns parsed String array with "U" prepended at index 0 if
     * read is unaligned and alignment should be attempted.
     * Returns parsed String array with "A" prepended at index 0 if
     * read has already been aligned.
     */
    public String[] parseAlignmentLine_SAM(String line) {
	try {
	    if (isHeaderLine(line)) return emptyStringArray;  // We have a header line
	    String[] parse_line = line.split("\t");
	    //String qname = parse_line[0];  // QNAME
	    int flag = Integer.parseInt(parse_line[1]);  // FLAG
	    String rname = parse_line[2];  // RNAME
	    int pos = Integer.parseInt(parse_line[3]);  // POS
	    //int mapq = Integer.parseInt(parse_line[4]);  // MAPQ
	    //String rnext = parse_line[6];  // RNEXT
	    //int pnext = Integer.parseInt(parse_line[7]);  // PNEXT
	    int tlen = Integer.parseInt(parse_line[8]);  // TLEN
	    String seq = parse_line[9];  // SEQ
	    String qual = parse_line[10];  // QUAL

	    if (flag >= 256) {  // Bad or chimeric read
		badReads.incrementAndGet();
		return emptyStringArray;
	    } else if ((((flag >>> 2) & 1) == 1) || (rname.equals("*")) || (pos == 0) || (((flag & 1) == 1) && (tlen == 0)) || isDenovo) {  // Unaligned read
		unalignedReads.incrementAndGet();
		if (seq.equals("*")) return emptyStringArray;
		if (isDenovo && (((flag >>> 2) & 1) == 0) && (((flag >>> 4) & 1) == 0)) parse_line[9] = reverseComplement(seq);
		parse_line[0] = "U" + parse_line[0];
		if (qual.length() == seq.length()) return parse_line;
		else if (qualities.length() != seq.length()) qualities = new String(new char[seq.length()]).replace('\0', '(');
		parse_line[10] = qualities;
		return parse_line;
	    } else if ((flag & 1) == 0) {  // Single-end aligned read
		int length = seq.length();
		boolean isRevComp = (((flag >>> 4) & 1) == 1);
		int repliconIndex = 0;
		if (replicons.containsKey(rname)) repliconIndex = replicons.get(rname);
		mapReadToCoords(repliconIndex, pos, length, isRevComp);
		alignedReads.incrementAndGet(repliconIndex);
		avgLengthReads.addAndGet(length);
		parse_line[0] = "A" + parse_line[0];
		return parse_line;
	    } else if (((flag & 1) == 1) && (tlen != 0) && (((flag >>> 6) & 1) == 1) && (((flag >>> 7) & 1) == 0) && (Math.abs(tlen) <= maxPairedEndLength)) {  // First of paired-end aligned read
		boolean isRevComp = (((flag >>> 4) & 1) == 1);
		int repliconIndex = 0;
		if (replicons.containsKey(rname)) repliconIndex = replicons.get(rname);
		mapReadToCoords(repliconIndex, pos, tlen, isRevComp);
		alignedReads.incrementAndGet(repliconIndex);
		avgLengthReads.addAndGet(Math.abs(tlen));
		parse_line[0] = "A" + parse_line[0];
		return parse_line;
	    } else if (((flag & 1) == 1) && (tlen != 0) && (((flag >>> 6) & 1) == 0) && (((flag >>> 7) & 1) == 1) && (Math.abs(tlen) <= maxPairedEndLength)) {  // Second of paired-end aligned read
		parse_line[0] = "A" + parse_line[0];
		return parse_line;
	    }
	    return emptyStringArray;
	} catch (Exception e) {
	    return emptyStringArray;
	}
    }

    /**
     * Parse one line of alignment information from the BAM file.
     * Return null if there is an error or if we reach the end
     * of the file. 
     * Returns empty String array if line should be skipped.
     * Returns parsed String array with "U" prepended at index 0 if
     * read is unaligned and alignment should be attempted.
     * Returns parsed String array with "A" prepended at index 0 if
     * read has already been aligned.
     */
    public String[] parseAlignmentLine_BAM() {
	try {
	    if (this.in == null) {
		this.in = new GZIPInputStream(new FileInputStream(fileName));
		disregardHeader_BAM(this.in);
	    }
	    int lengthOfAlignmentSection = readInt(in);
	    int repliconIndex = readInt(in);
	    int pos = readInt(in) + 1;
	    int bin_mq_nl = readInt(in);
	    //int bin = bin_mq_nl >>> 16;
	    int mapq = (bin_mq_nl >>> 8) & 0xFF;
	    int l_read_name = bin_mq_nl & 0xFF;
	    int flag_nc = readInt(in);
	    int flag = flag_nc >>> 16;
	    int numCigarOps = flag_nc & 0xFFFF;
	    int seqLength = readInt(in);
	    int nextRefID = readInt(in);
	    int pnext = readInt(in) + 1;
	    int tlen = readInt(in);
	    String readName = readString(in, l_read_name);
	    String qname = readName.substring(0, readName.length()-1);
	    String cigar = readCigarString(in, numCigarOps);
	    String seq = readSequence(in, seqLength);
	    String qual = readString(in, seqLength);
	    int partialLengthOfAlignmentSection = 8*4 + l_read_name + numCigarOps*4 + ((seqLength+1)/2) + seqLength;
	    String remainder = readString(in, lengthOfAlignmentSection - partialLengthOfAlignmentSection);

	    if (flag >= 256) {  // Bad or chimeric read
		badReads.incrementAndGet();
		return emptyStringArray;
	    } else if ((((flag >>> 2) & 1) == 1) || (repliconIndex < 0) || (pos == 0) || (((flag & 1) == 1) && (tlen == 0)) || isDenovo) {  // Unaligned read
		unalignedReads.incrementAndGet();
		if (seq.length() <= 1) return emptyStringArray;
		if (isDenovo && (((flag >>> 2) & 1) == 0) && (((flag >>> 4) & 1) == 0)) seq = reverseComplement(seq);
		String[] parse_line = new String[11];
		parse_line[0] = "U" + qname;  // QNAME; Sentinel: "unaligned"
		parse_line[1] = flag + "";  // FLAG
		parse_line[2] = repliconIndex + "";  // RNAME
		parse_line[3] = pos + "";  // POS
		parse_line[4] = mapq + "";  // MAPQ
		parse_line[5] = cigar;  // CIGAR
		parse_line[6] = nextRefID +  "";  // RNEXT
		parse_line[7] = pnext + "";  // PNEXT
		parse_line[8] = tlen + "";  // TLEN
		parse_line[9] = seq;  // SEQ
		parse_line[10] = qual;  // QUAL
		if ((qual.length() == seq.length()) && (qual.indexOf('\t') == -1)) return parse_line;
		else if (qualities.length() != seq.length()) qualities = new String(new char[seq.length()]).replace('\0', '(');
		parse_line[10] = qualities;
		return parse_line;
	    } else if ((flag & 1) == 0) {  // Single-end aligned read
		int length = seq.length();
		boolean isRevComp = (((flag >>> 4) & 1) == 1);
		mapReadToCoords(repliconIndex, pos, length, isRevComp);
		alignedReads.incrementAndGet(repliconIndex);
		avgLengthReads.addAndGet(length);
		String[] parse_line = new String[11];
		parse_line[0] = "A" + qname;  // QNAME; Sentinel: "aligned"
		parse_line[1] = flag + "";  // FLAG
		parse_line[2] = repliconIndex + "";  // RNAME
		parse_line[3] = pos + "";  // POS
		parse_line[4] = mapq + "";  // MAPQ
		parse_line[5] = cigar;  // CIGAR
		parse_line[6] = nextRefID +  "";  // RNEXT
		parse_line[7] = pnext + "";  // PNEXT
		parse_line[8] = tlen + "";  // TLEN
		parse_line[9] = seq;  // SEQ
		parse_line[10] = qual;  // QUAL
		return parse_line;
	    } else if (((flag & 1) == 1) && (tlen != 0) && (((flag >>> 6) & 1) == 1) && (((flag >>> 7) & 1) == 0) && (Math.abs(tlen) <= maxPairedEndLength)) {  // First of paired-end aligned read
		boolean isRevComp = (((flag >>> 4) & 1) == 1);
		mapReadToCoords(repliconIndex, pos, tlen, isRevComp);
		alignedReads.incrementAndGet(repliconIndex);
		avgLengthReads.addAndGet(Math.abs(tlen));
		String[] parse_line = new String[11];
		parse_line[0] = "A" + qname;  // QNAME; Sentinel: "aligned"
		parse_line[1] = flag + "";  // FLAG
		parse_line[2] = repliconIndex + "";  // RNAME
		parse_line[3] = pos + "";  // POS
		parse_line[4] = mapq + "";  // MAPQ
		parse_line[5] = cigar;  // CIGAR
		parse_line[6] = nextRefID +  "";  // RNEXT
		parse_line[7] = pnext + "";  // PNEXT
		parse_line[8] = tlen + "";  // TLEN
		parse_line[9] = seq;  // SEQ
		parse_line[10] = qual;  // QUAL
		return parse_line;
	    } else if (((flag & 1) == 1) && (tlen != 0) && (((flag >>> 6) & 1) == 0) && (((flag >>> 7) & 1) == 1) && (Math.abs(tlen) <= maxPairedEndLength)) {  // Second of paired-end aligned read
		String[] parse_line = new String[11];
		parse_line[0] = "A" + qname;  // QNAME; Sentinel: "aligned"
		parse_line[1] = flag + "";  // FLAG
		parse_line[2] = repliconIndex + "";  // RNAME
		parse_line[3] = pos + "";  // POS
		parse_line[4] = mapq + "";  // MAPQ
		parse_line[5] = cigar;  // CIGAR
		parse_line[6] = nextRefID +  "";  // RNEXT
		parse_line[7] = pnext + "";  // PNEXT
		parse_line[8] = tlen + "";  // TLEN
		parse_line[9] = seq;  // SEQ
		parse_line[10] = qual;  // QUAL
		return parse_line;
	    }
	    return emptyStringArray;
	} catch (Exception e) {
	    try {
		if (this.in != null) {
		    in.close();
		    in = null;
		}
	    } catch (IOException e2) {
		return null;
	    }
	    return null;
	}
    }



    /**************************************************
     **********   PRIVATE INSTANCE METHODS   **********
     **************************************************/

    /**
     * Read in the first lines of the SAM file, extract useful header info,
     * and check validity.
     */
    private boolean prepSAM(long numLines) {
	cigarMap = new int[91];  // Converts cigar characters to counts
	File f = new File(fileName);
	if (!f.exists()) {  // File does not exist
	    //output("Error - file " + fileName + " does not exist.\n");
	    return false;
	}
	try {
	    Scanner reader = new Scanner(f);
	    long lineCount = (long)0;
	    int totalAlignmentLines = 0;  // Total number of non-header lines
	    while (reader.hasNextLine() && (lineCount < numLines)) {
		String line = reader.nextLine();
		boolean isValidLine;
		if (isHeaderLine(line)) isValidLine = parseHeaderLine_SAM(line);
		else {
		    isValidLine = checkAlignmentLine_SAM(line);
		    totalAlignmentLines++;
		}
		if (!isValidLine) {
		    if (pairedEndCount == totalAlignmentLines) isPairedEnd = true;
		    return false;
		}
		lineCount++;
	    }
	    reader.close();
	    if (pairedEndCount == totalAlignmentLines) isPairedEnd = true;
	    return true;
	} catch (FileNotFoundException e) {
	    //output("Error - could not read in from file " + fileName + "\n");
	    return false;
	}
    }

    /**
     * Read in the first lines of the BAM file, extract useful header info,
     * and check validity.
     */
    private boolean prepBAM(long numLines) {
	File f = new File(fileName);
	if (!f.exists()) {  // File does not exist
	    //output("Error - file " + fileName + " does not exist.\n");
	    return false;
	}
	try {
	    this.in = new GZIPInputStream(new FileInputStream(fileName));

	    // Read in header info
	    String magicString = readString(in, 4);
	    if (!magicString.substring(0, 3).toLowerCase().equals("bam")) {  // Not a bam file
		in.close();
		in = null;
		return false;
	    }
	    int headerTextLength = readInt(in);
	    String samHeaderText = readString(in, headerTextLength);
	    int numRefSeqs = readInt(in);
	    // Read in each replicon (reference sequence) name and length
	    for (int i=0; i<numRefSeqs; i++) {
		int repliconNameLength = readInt(in);
		String repliconName = readString(in, repliconNameLength);
		int repliconLength = readInt(in);
		replicons.put(repliconName, replicons.size());
		repliconLengths.add(repliconLength);
	    }

	    long lineCount = (long)0;
	    boolean isValidLine = true;
	    while (isValidLine && (lineCount < numLines)) {
		lineCount++;
		isValidLine = checkAlignmentLine_BAM();
		if (!isValidLine) {
		    if (pairedEndCount == lineCount) isPairedEnd = true;
		    in.close();
		    in = null;
		    return false;
		}
	    }
	    in.close();
	    in = null;
	    if (pairedEndCount == lineCount) isPairedEnd = true;
	    return true;
	} catch (IOException e) {
	    //output("Error - could not read in from file " + fileName + "\n");
	    try { in.close();
	    } catch (Exception e2) { }
	    in = null;
	    return false;
	}
    }

    /**
     * Read in the entire SAM file.
     */
    private boolean readInSAM() {
	try {
	    File f = new File(fileName);
	    Scanner reader = new Scanner(f);
	    while (reader.hasNextLine()) {
		parseAlignmentLine_SAM(reader.nextLine());
	    }
	    reader.close();
	    return true;
	} catch (FileNotFoundException e) {
	    output("Error - could not read in from file " + fileName + "\n");
	    return false;
	}
    }

    /**
     * Read in the entire BAM file.
     */
    private void readInBAM() {
	while (parseAlignmentLine_BAM() != null) { }
    }

    /**
     * Parse one line of header information from the SAM file.
     */
    private boolean parseHeaderLine_SAM(String line) {
	try {
	    String[] parse_line = line.split("\t");
	    if (parse_line[0].equals("@HD")) {
		// Information in header line is currently unused
	    } else if (parse_line[0].equals("@SQ")) {
		String repliconName = null;
		Integer repliconLength = null;
		for (int i=1; i<parse_line.length; i++) {
		    String[] pair = parse_line[i].split(":");
		    String tag = pair[0];
		    String value = pair[1];
		    if (tag.equals("SN")) repliconName = value;
		    if (tag.equals("LN")) repliconLength = Integer.parseInt(value);
		}
		if ((repliconName != null) && (repliconLength != null)) {
		    replicons.put(repliconName, replicons.size());
		    repliconLengths.add(repliconLength);
		} else {
		    //output("Error - reference sequence header line in SAM file is not formatted as expected:\n");
		    //output("\t" + line + "\n\n");
		    return false;
		}
	    } else if (parse_line[0].equals("@RG")) {
		// Information in read group line is currently unused
	    } else if (parse_line[0].equals("@PG")) {
		// Information in program line is currently unused
	    } else if (parse_line[0].equals("@CO")) {
		// Information in comment line is currently unused
	    } else {  // Unrecognized header record type
		//output("Error - unrecognized SAM record type " + parse_line[0] + " in header line:\n");
		//output("\t" + line + "\n\n");
		return false;
	    }
	    return true;
	} catch (Exception e) {
	    //output("Error - line not in expected SAM format:\n");
	    //output("\t" + line + "\n\n");
	    return false;
	}
    }

    /**
     * Read in and disregard header of BAM file.
     */
    private void disregardHeader_BAM(GZIPInputStream in) throws IOException {
	String magicString = readString(in, 4);
	int headerTextLength = readInt(in);
	String samHeaderText = readString(in, headerTextLength);
	int numRefSeqs = readInt(in);

	// Read in each replicon (reference sequence) name and length
	for (int i=0; i<numRefSeqs; i++) {
	    int repliconNameLength = readInt(in);
	    String repliconName = readString(in, repliconNameLength);
	    int repliconLength = readInt(in);
	}
    }

    /**
     * Parse one line of alignment information from the SAM file and check its validity.
     */
    private boolean checkAlignmentLine_SAM(String line) {
	try {
	    String[] parse_line = line.split("\t");
	    //String qname = parse_line[0];  // QNAME
	    int flag = Integer.parseInt(parse_line[1]);  // FLAG
	    String rname = parse_line[2];  // RNAME
	    int pos = Integer.parseInt(parse_line[3]);  // POS
	    //int mapq = Integer.parseInt(parse_line[4]);  // MAPQ
	    //String rnext = parse_line[6];  // RNEXT
	    //int pnext = Integer.parseInt(parse_line[7]);  // PNEXT
	    int tlen = Integer.parseInt(parse_line[8]);  // TLEN
	    String seq = parse_line[9];  // SEQ
	    String qual = parse_line[10];  // QUAL
	    if ((flag & 1) == 1) pairedEndCount++;

	    if (flag >= 256) {  // Bad or chimeric read
		// Do nothing
	    } else if ((((flag >>> 2) & 1) == 1) || (rname.equals("*")) || (pos == 0) || (((flag & 1) == 1) && (tlen == 0))) {  // Unaligned read
		if ((!seq.equals("*")) && (!qual.equals("*")) && (seq.length() != qual.length())) {
		    //output("Error - unrecognized format in SAM file - expecting seq length to equal qual length:\n");
		    //output("\t" + line + "\n\n");
		    return false;
		}
	    } else if ((flag & 1) == 0) {  // Single-end aligned read
		if (!rname.equals("*") && (replicons.size() > 0) && (!replicons.containsKey(rname))) {
		    //output("Error - unrecognized format in SAM file - expecting RNAME to match a SQ-SN tag:\n");
		    //output("\t" + line + "\n\n");
		    return false;
		}
		parseCigar(parse_line[5]);  // CIGAR
		int cigarLength = (cigarMap['M'] + cigarMap['I'] + cigarMap['S'] + cigarMap['='] + cigarMap['X']);
		if (!seq.equals("*") && (seq.length() != cigarLength)) {
		    //output("Error - unrecognized format in SAM file - expecting sequence length to match CIGAR length:\n");
		    //output("\t" + line + "\n\n");
		    return false;
		}
	    } else if (((flag & 1) == 1) && (tlen > 0) && ((((flag >>> 6) & 1) == 0) || (((flag >>> 7) & 1) == 0))) {  // Paired-end aligned read (leftmost mate-pair)
		if (!rname.equals("*") && (replicons.size() > 0) && (!replicons.containsKey(rname))) {
		    //output("Error - unrecognized format in SAM file - expecting RNAME to match a SQ-SN tag:\n");
		    //output("\t" + line + "\n\n");
		    return false;
		}
	    }
	    return true;
	} catch (Exception e) {
	    //output("Error - line not in expected SAM format:\n");
	    //output("\t" + line + "\n\n");
	    return false;
	}
    }

    /**
     * Parse one line of alignment information from the BAM file and check its validity.
     */
    private boolean checkAlignmentLine_BAM() {
	try {
	    int lengthOfAlignmentSection = readInt(in);
	    int repliconIndex = readInt(in);
	    int pos = readInt(in) + 1;
	    int bin_mq_nl = readInt(in);
	    //int bin = bin_mq_nl >>> 16;
	    //int mapq = (bin_mq_nl >>> 8) & 0xFF;
	    int l_read_name = bin_mq_nl & 0xFF;
	    int flag_nc = readInt(in);
	    int flag = flag_nc >>> 16;
	    int numCigarOps = flag_nc & 0xFFFF;
	    int seqLength = readInt(in);
	    int nextRefID = readInt(in);
	    int pnext = readInt(in) + 1;
	    int tlen = readInt(in);
	    String readName = readString(in, l_read_name);
	    //String qname = readName.substring(0, readName.length()-1);
	    for (int i=0; i<numCigarOps; i++) readInt(in);  // Ignore CIGAR ops
	    String seq = readSequence(in, seqLength);
	    String qual = readString(in, seqLength);
	    int partialLengthOfAlignmentSection = 8*4 + l_read_name + numCigarOps*4 + ((seqLength+1)/2) + seqLength;
	    String remainder = readString(in, lengthOfAlignmentSection - partialLengthOfAlignmentSection);

	    if ((flag & 1) == 1) pairedEndCount++;
	    //this.maxReadLength = Math.max(this.maxReadLength, seq.length());
	    for (int i=0; i<qual.length(); i++)
		this.maxQuality = Math.max(this.maxQuality, (int)(qual.charAt(i)));

	    if (flag >= 256) {  // Bad or chimeric read
		// Do nothing
	    } else if ((((flag >>> 2) & 1) == 1) || (repliconIndex < 0) || (pos == 0) || (((flag & 1) == 1) && (tlen == 0))) {  // Unaligned read
		if ((!seq.equals("*")) && (!qual.equals("*")) && (seq.length() != qual.length())) {
		    //output("Error - unrecognized format in BAM file - expecting seq length to equal qual length:\n");
		    //output("\t" + line + "\n\n");
		    return false;
		}
	    } else if ((flag & 1) == 0) {  // Single-end aligned read
		if (repliconIndex >= replicons.size()) {
		    //output("Error - unrecognized format in BAM file - expecting valid replicon index:\n");
		    //output("\t" + line + "\n\n");
		    return false;
		}
	    } else if (((flag & 1) == 1) && (tlen > 0) && ((((flag >>> 6) & 1) == 0) || (((flag >>> 7) & 1) == 0))) {  // Paired-end aligned read (leftmost mate-pair)
		if (repliconIndex >= replicons.size()) {
		    //output("Error - unrecognized format in BAM file - expecting RNAME to match a SQ-SN tag:\n");
		    //output("\t" + line + "\n\n");
		    return false;
		}
	    }
	    return true;
	} catch (Exception e) {
	    return false;
	}
    }

    /**
     * Parse a CIGAR string from a SAM file and populate the cigarMap.
     */
    private void parseCigar(String s) {
	for (char ch : cigarChars) cigarMap[ch] = 0;  // Initialize cigar map
	if (s.equals("*")) return;
	int index = 0;
	while (index < s.length()) {
	    StringBuilder sb = new StringBuilder();
	    while (Character.isDigit(s.charAt(index))) {
		sb.append(s.charAt(index));
		index++;
	    }
	    cigarMap[s.charAt(index)] += Integer.parseInt(sb.toString());
	    index++;
	}
    }

    /**
     * Read in CIGAR information from a BAM file and return the CIGAR string.
     */
    private String readCigarString(GZIPInputStream in, int numCigarOps) throws IOException {
	String cigar = "";
	for (int i=0; i<numCigarOps; i++) {
	    int bin_cigar_op = readInt(in);
	    int op = bin_cigar_op & 0xF;
	    if (op > 9) return "";
	    int op_len = bin_cigar_op >>> 4;
	    cigar += op_len + "" + cigarChars[op];
	}
	return cigar;
    }

    /**
     * Map a read to the appropriate loci in the appropriate coords array.
     */
    private void mapReadToCoords(int repliconIndex, int pos, int length, boolean isRevComp) {
	if (length < 0) {
	    pos = pos - length + 1;
	    if (pos < 0) pos = 0;
	    length = -length;
	}
	if (coordinates_plus[repliconIndex].length() < (pos+length)) {  // Increase size of arrays
	    synchronized (coordinates_plus) {
		if (coordinates_plus[repliconIndex].length() < (pos+length)) {  // Double check if we want to increase size of arrays
		    int new_size = Math.max(pos+length, 2*coordinates_plus[repliconIndex].length());
		    AtomicIntegerArray temp_plus = new AtomicIntegerArray(new_size);
		    AtomicIntegerArray temp_minus = new AtomicIntegerArray(new_size);
		    for (int i=0; i<coordinates_plus[repliconIndex].length(); i++) {
			temp_plus.set(i, coordinates_plus[repliconIndex].get(i));
			temp_minus.set(i, coordinates_minus[repliconIndex].get(i));
		    }
		    coordinates_plus[repliconIndex] = temp_plus;
		    coordinates_minus[repliconIndex] = temp_minus;
		}
	    }
	}
	if (!isRevComp) {  // Plus strand
	    for (int i=pos; i<pos+length; i++) coordinates_plus[repliconIndex].incrementAndGet(i);
	} else {  // Minus strand
	    for (int i=pos; i<pos+length; i++) coordinates_minus[repliconIndex].incrementAndGet(i);
	}
    }



    /**********************************************
     **********   PUBLIC CLASS METHODS   **********
     **********************************************/

    /**
     * Returns true if line is a SAM header line, false otherwise.
     */
    public static boolean isHeaderLine(String line) {
	return line.startsWith("@");
    }

    /**
     * Returns true if the file is a text file.
     * Returns false if the file is a binary file (e.g., 
     * gzip or bam).
     */
    public static boolean isTextFile(String fileName) {
	try {
	    FileInputStream in = new FileInputStream(new File(fileName));
	    int size = in.available();
	    if (size > 1024) size = 1024;
	    byte[] data = new byte[size];
	    in.read(data);
	    in.close();

	    int ascii = 0;
	    int other = 0;
	    for (int i=0; i<data.length; i++) {
		byte b = data[i];
		if (b < 0x09) return false;
		if ((b == 0x09) || (b == 0x0A) || (b == 0x0C) || (b == 0x0D)) ascii++;
		else if ((b >= 0x20) && (b <= 0x7E)) ascii++;
		else other++;
	    }
	    return ((100 * ascii / size) > 95);
	} catch (IOException e) {
	    return false;
	}
    }



    /***********************************************
     **********   PRIVATE CLASS METHODS   **********
     ***********************************************/

    private static void output(String s) {
	System.err.print(s);
    }

    /**
     * Converts the first four bytes in the specified array into an integer.
     * Assumes little-endianness.
     */
    private static int byteArrayToInt(byte[] b) {
	if (b.length < 4) return 0;
	//return (b[0] << 24) + ((b[1] & 0xFF) << 16) + ((b[2] & 0xFF) << 8) + (b[3] & 0xFF);
	return ((b[3] & 0xFF) << 24) + ((b[2] & 0xFF) << 16) + ((b[1] & 0xFF) << 8) + (b[0] & 0xFF);
    }

    /**
     * Reads an integer (little-endian) from a gzip/bam file.
     */
    private static int readInt(GZIPInputStream in) throws IOException {
	in.read(buffer1);
	intBuffer[0] = buffer1[0];
	in.read(buffer1);
	intBuffer[1] = buffer1[0];
	in.read(buffer1);
	intBuffer[2] = buffer1[0];
	in.read(buffer1);
	intBuffer[3] = buffer1[0];
	return byteArrayToInt(intBuffer);
    }

    /**
     * Reads a String of the specified length from a gzip/bam file.
     */ 
    private static String readString(GZIPInputStream in, int length) throws IOException {
	byte[] s = new byte[length];
	for (int i=0; i<length; i++) {
	    in.read(buffer1);
	    s[i] = buffer1[0];
	}
	return (new String(s));
    }

    /**
     * Reads a nucleotide sequence of the given length from a gzip/bam file.
     */
    private static String readSequence(GZIPInputStream in, int length) throws IOException {
	char[] sequence = new char[length];
	int NT_index = 15;
	for (int i=0; i<(length+1)/2; i++) {
	    in.read(buffer1);
	    NT_index = (buffer1[0] >>> 4) & 0xF;
	    sequence[2*i] = NT_code[NT_index];
	    NT_index = buffer1[0] & 0xF;
	    if (2*i+1 < length) sequence[2*i+1] = NT_code[NT_index];
	}
	return (new String(sequence));
    }

    /**
     * Mergesort parallel arrays "a" and "b" based on values in "a"
     */
    private static void mergesort(int[] a, String[] b, int lo, int hi) {
	if (lo < hi) {
	    int q = (lo+hi)/2;
	    mergesort(a, b, lo, q);
	    mergesort(a, b, q+1, hi);
	    merge(a, b, lo, q, hi);
	}
    }

    /**
     * Mergesort helper method
     */
    private static void merge(int[] a, String[] b, int lo, int q, int hi) {
	int[] a1 = new int[q-lo+1];
	int[] a2 = new int[hi-q];
	String[] b1 = new String[q-lo+1];
	String[] b2 = new String[hi-q];
	for (int i=0; i<a1.length; i++) {
	    a1[i] = a[lo+i];
	    b1[i] = b[lo+i];
	}
	for (int j=0; j<a2.length; j++) {
	    a2[j] = a[q+1+j];
	    b2[j] = b[q+1+j];
	}
	int i=0;  // Index for first half arrays
	int j=0;  // Index for second half arrays
	for (int k=lo; k<=hi; k++) {
	    if (i >= a1.length) {
		a[k] = a2[j];
		b[k] = b2[j];
		j++;
	    } else if (j >= a2.length) {
		a[k] = a1[i];
		b[k] = b1[i];
		i++;
	    } else if (a1[i] <= a2[j]) {
		a[k] = a1[i];
		b[k] = b1[i];
		i++;
	    } else {
		a[k] = a2[j];
		b[k] = b2[j];
		j++;
	    }
	}
    }

    /**
     * Returns a reverse complemented version of String s.
     */
    private static String reverseComplement(String s) {
	StringBuilder sb = new StringBuilder(s.length());
	for (int i=s.length()-1; i>=0; i--) {
	    if (s.charAt(i) == 'A') sb.append('T');
	    else if (s.charAt(i) == 'C') sb.append('G');
	    else if (s.charAt(i) == 'G') sb.append('C');
	    else if (s.charAt(i) == 'T') sb.append('A');
	    else sb.append(s.charAt(i));
	}
	return sb.toString();
    }



    /*************************************
     **********   MAIN METHOD   **********
     *************************************/

    public static void main(String[] args) {
	if (args.length < 1) {
	    System.err.println("\n" + "SamOps requires a SAM/BAM file as a command line argument and it outputs information about the SAM/BAM file." + "\n");
	    System.err.println("\t" + "java SampOps <foo.sam or foo.bam>" + "\n");
	    System.exit(0);
	}
	long numLines = 10000;
	SamOps S = new SamOps(args[0], numLines, false);  // Read in first few lines to check validity
	if (S.isValidSam()) {
	    System.out.println("\n" + "File appears to be in valid SAM format.");
	    System.out.println("Does file contain paired-end reads?\t" + S.isPairedEnd() + "\n");
	    S.readInSAM();
	    System.out.println(S.toString());
	}
	if (S.isValidBam()) {
	    System.out.println("\n" + "File appears to be in valid BAM format.");
	    System.out.println("Does file contain paired-end reads?\t" + S.isPairedEnd() + "\n");
	    S.readInBAM();
	    System.out.println(S.toString());
	}
    }

}

