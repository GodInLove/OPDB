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
import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLong;
import javax.swing.JTextArea;

public class Peregrine {

    /********************************************
     **********   Instance Variables   **********
     ********************************************/

    private Index index;
    private AtomicIntegerArray[] coordinates_plus;   // Number of reads mapping to each coordinate in each sequence
    private AtomicIntegerArray[] coordinates_minus;  // Number of reads mapping to each coordinate in each sequence
    private int PHRED_offset;
    private BufferedReader reader;
    private BufferedReader reader2;  // Mate-pair file if using paired-end reads
    private PrintWriter samWriter;  // Output reads to SAM file
    private boolean isCompressedFile;  // True if gzip file. False otherwise.
    private SamOps sam;  // Used if reads file is in SAM/BAM format
    private boolean firstLineDone;  // Used when reading in the first line of FASTA files
    private String parametersString;
    private String[] genomeNames;  // Used when outputting SAM files
    private long runningTime;  // Time taken to execute program

    // Results
    public AtomicInteger totalReads;
    private AtomicInteger invalidQualityReads;      // Reads containing an invalid quality score
    private AtomicInteger ambiguousReads;           // Reads with ambiguous nucleotides
    private AtomicInteger lowQualityReads;          // Reads with quality scores too low for mapping
    private AtomicInteger mappedOnceReads;          // Reads mapping to one location
    private AtomicInteger mappedMoreThanOnceReads;  // Reads mapping to multiple locations
    public AtomicIntegerArray exactMappedReads;     // Reads mapping without mismatches/errors to each sequence
    public AtomicIntegerArray inexactMappedReads;   // Reads mapping with mismatches/errors to each sequence
    public AtomicLong avgLengthReads;               // Avg length of mapping reads
    private AtomicInteger numMappingReads;          // Number of reads mapping
    private AtomicInteger nameCount;                // Used when outputting SAM file and reads have no name



    /*****************************************
     **********   Class Variables   **********
     *****************************************/

    // Parameters
    public static JTextArea output;  // Output to GUI or, if null, to System.out
    public static PrintWriter summaryWriter = null;  // For outputting summary file
    public static int numSequences;
    private static int readsInputFormat;  // 0=FASTQ, 1=QSEQ, 2=FASTA, 3=SAM(single), 4=SAM(paired), 5=BAM(single), 6=BAM(paired)
    public static String[] sequenceFiles;
    public static String[][] annotationsPlus = null;  // Gene/RNA coords for each replicon
    public static String[][] annotationsMinus = null;  // Gene/RNA coords for each replicon
    public static String readsFile;
    public static double percentMismatches = 0.15;
    public static boolean stopAfterOneHit = true;
    public static double percentSeedLength = 0.33;  // Minimum allowed seed (read) length
    public static int minReadLength = 20;
    public static boolean isPairedEnd = false;  // Do we have paired-end reads?
    public static String pairedEndFile;  // Name of mate-pair file if using paired-end reads
    public static int maxPairedEndLength = 500;
    public static boolean singleEndOrientationReverseComplement = false;  // RevComp reads
    public static String pairedEndOrientation = "fr";  // ff, fr, rf, rr
    public static int numThreads = 1;  // Number of threads
    public static boolean outputSAM = false;  // Output reads to SAM file
    public static String outputDIR = "Peregrine_Results/";  // Output directory
    private static long randomSeed = -1;  // Seed for the random number generator
    private static boolean time = false;  // Output time taken to execute program

    public static String browserDIR = "genomeBrowserFiles/";  // Files for genome browser
    public static ArrayList<String> wigFiles;  // List of comma-delimited set of wig files for IGV
    public static boolean outputBrowserFile = true;  // Should we output a browser file?



    /**************************************
     **********   Constructors   **********
     **************************************/

    public Peregrine() {

	runningTime = System.currentTimeMillis();
	this.coordinates_plus = new AtomicIntegerArray[Peregrine.numSequences];
	this.coordinates_minus = new AtomicIntegerArray[Peregrine.numSequences];
	this.exactMappedReads = new AtomicIntegerArray(Peregrine.numSequences);
	this.inexactMappedReads = new AtomicIntegerArray(Peregrine.numSequences);
	if (wigFiles == null) {
	    wigFiles = new ArrayList<String>(numSequences);
	    for (int i=0; i<numSequences; i++) wigFiles.add("");
	}
	parametersString = parametersToString();
	// For paired-end reads, we must consider all optimal alignments
	if (Peregrine.isPairedEnd) Peregrine.stopAfterOneHit = false;

	// Check if reads hits have been pre-computed
	boolean isValid = true;
	for (int i=0; i<numSequences; i++) {
	    String inFileName = getCompressedFileName(readsFile, sequenceFiles[i]);
	    FileOps foo = new FileOps(inFileName, readsFile, sequenceFiles[i]);
	    if (foo.isValid()) {
		coordinates_plus[i] = new AtomicIntegerArray(foo.getCoordinates_plus());
		coordinates_minus[i] = new AtomicIntegerArray(foo.getCoordinates_minus());
		totalReads = new AtomicInteger(foo.getTotalReads());
		invalidQualityReads = new AtomicInteger(foo.getInvalidQualityReads());
		ambiguousReads = new AtomicInteger(foo.getAmbiguousReads());
		lowQualityReads = new AtomicInteger(foo.getLowQualityReads());
		mappedOnceReads = new AtomicInteger(foo.getMappedOnceReads());
		mappedMoreThanOnceReads = new AtomicInteger(foo.getMappedMoreThanOnceReads());
		avgLengthReads = new AtomicLong(foo.getAvgLengthReads());
		numMappingReads = new AtomicInteger(0);
		nameCount = new AtomicInteger(1);
		exactMappedReads.set(i, foo.getExactMappedReads());
		inexactMappedReads.set(i, foo.getInexactMappedReads());
		runningTime = System.currentTimeMillis() - runningTime;
		if (i == numSequences-1) outputResults();
	    } else {
		isValid = false;
		totalReads = new AtomicInteger(0);
		invalidQualityReads = new AtomicInteger(0);
		ambiguousReads = new AtomicInteger(0);
		lowQualityReads = new AtomicInteger(0);
		mappedOnceReads = new AtomicInteger(0);
		mappedMoreThanOnceReads = new AtomicInteger(0);
		avgLengthReads = new AtomicLong(0);
		numMappingReads = new AtomicInteger(0);
		nameCount = new AtomicInteger(1);
		exactMappedReads.set(i, 0);
		inexactMappedReads.set(i, 0);
	    }
	    if (wigFiles.get(i).length() > 0) wigFiles.set(i, wigFiles.get(i) + ",");
	    wigFiles.set(i, wigFiles.get(i) + outputDIR + browserDIR + inFileName + ".plus" + ".wig");
	    wigFiles.set(i, wigFiles.get(i) + "," + outputDIR + browserDIR + inFileName + ".minus" + ".wig");
	}

	if (isValid && outputSAM) output("Rockhopper did not output a SAM file since sequencing read alignments were previously cached. If you want a SAM file to be output, you must first clear the cache." + "\n\n");

	if (!isValid) {
	    this.index = new Index(sequenceFiles);
	    readsInputFormat = getReadsInputFileType(readsFile);
	    int numThreads_OLD = numThreads;
	    if ((readsInputFormat == 5) || (readsInputFormat == 6)) numThreads = 1;  // BAM file
	    for (int i=0; i<numSequences; i++) {
		coordinates_plus[i] = new AtomicIntegerArray(index.getReplicon(i).getLength());
		coordinates_minus[i] = new AtomicIntegerArray(index.getReplicon(i).getLength());
	    }
	    this.PHRED_offset = getPhredOffset(readsFile);
	    firstLineDone = false;
	    if (outputSAM) outputSamHeader();

	    // Set parameters
	    index.setStopAfterOneHit(Peregrine.stopAfterOneHit);
	    if (randomSeed >= 0) Index.setRandomSeed(Peregrine.randomSeed);

	    mapReads();
	    runningTime = System.currentTimeMillis() - runningTime;
	    numThreads = numThreads_OLD;  // Reset the number of threads
	    if (numMappingReads.get() > 0) avgLengthReads.set(avgLengthReads.get() / numMappingReads.get());
	    avgLengthReads.set(Math.max(avgLengthReads.get(), 72));  // Minimum of 72 nts

	    // Incorporate results from SAM file
	    if (sam != null) {
		int samAlignedReads = 0;
		for (int i=0; i<sam.alignedReads.length(); i++) samAlignedReads += sam.alignedReads.get(i);
		if (samAlignedReads > 0) avgLengthReads.set(Math.max(avgLengthReads.get(), sam.avgLengthReads.get() / samAlignedReads));
		totalReads.addAndGet(sam.badReads.get());
		invalidQualityReads.addAndGet(sam.badReads.get());
		for (int i=0; i<Math.min(exactMappedReads.length(), sam.alignedReads.length()); i++) {
		    totalReads.addAndGet(sam.alignedReads.get(i));
		    exactMappedReads.addAndGet(i, sam.alignedReads.get(i));
		    for (int j=0; j<Math.min(coordinates_plus[i].length(), sam.coordinates_plus[i].length()); j++) {
			coordinates_plus[i].addAndGet(j, sam.coordinates_plus[i].get(j));
			coordinates_minus[i].addAndGet(j, sam.coordinates_minus[i].get(j));
		    }
		}
	    }

	    if (samWriter != null) samWriter.close();
	    outputResults();
	    outputReadsToFile();
	}
    }



    /*************************************************
     **********   Public Instance Methods   **********
     *************************************************/

    /**
     * Read from reads file and map any reads from file to 
     * reference index.
     */
    public void processReads() {
	String line = "";;
	String name = "";
	String nextName = "";  // Special case (FASTA)
	String read = "";
	String quality = "";
	String quality_fasta = "";  // Special case (FASTA)
	int[] qualities = new int[0];
	String[] parse_line = null;
	ArrayList<Hit> hits = new ArrayList<Hit>();
	ArrayList<Hit> seedHits = new ArrayList<Hit>();
	Alignment a = new Alignment();
	// Using paired-end reads
	String line2 = "";;
	String name2 = "";
	String nextName2 = "";  // Special case (FASTA)
	String read2 = "";
	String quality2 = "";
	String quality_fasta2 = "";  // Special case (FASTA)
	int[] qualities2 = new int[0];
	String[] parse_line2 = null;
	ArrayList<Hit> hits2 = new ArrayList<Hit>();
	ArrayList<Hit> seedHits2 = new ArrayList<Hit>();
	ArrayList<Hit> combinedHits = new ArrayList<Hit>();
	Alignment a2 = new Alignment();

	try {
	    synchronized(reader) {
		if ((readsInputFormat == 2) && !firstLineDone) {  // Special case (FASTA)
		    line = reader.readLine();  // Special case (FASTA)
		    if (line.startsWith(">")) nextName = line.substring(1);
		    if (isPairedEnd) {
			line2 = reader2.readLine();
			if (line2.startsWith(">")) nextName2 = line2.substring(1);
		    }
		}
		firstLineDone = true;
	    }

	    while (line != null) {

		synchronized(reader) {  // Thread safe reading from file
		    if (readsInputFormat < 5) line = reader.readLine();  // Not a BAM file
		    if (line == null) continue;
		    if (readsInputFormat == 0) {  // FASTQ file
			name = line;
			read = reader.readLine();
			line = reader.readLine();
			quality = reader.readLine();
		    } else if (readsInputFormat == 1) {  // QSEQ file
			parse_line = line.split("\t");
			name = parse_line[0] + ":" + parse_line[2] + ":" + parse_line[3] + ":" + parse_line[4] + ":" + parse_line[5] + ":" + parse_line[6] + ":" + parse_line[7];
			read = parse_line[8];
			quality = parse_line[9];
		    } else if (readsInputFormat == 2) {  // FASTA file
			name = nextName;
			read = "";
			while ((line != null) && ((line.length() == 0) || (line.charAt(0) != '>'))) {
			    read += line;
			    line = reader.readLine();
			}
			if (line == null) nextName = "";
			else if (line.startsWith(">")) nextName = line.substring(1);
			else nextName = "";  // Case should never be reached
			if (quality_fasta.length() != read.length()) quality_fasta = new String(new char[read.length()]).replace('\0', '(');
			quality = quality_fasta;
		    } else if (readsInputFormat == 3) {  // SAM file (single-end)
			if (SamOps.isHeaderLine(line)) continue;
		    } else if (readsInputFormat == 4) {  // SAM file (paired-end)
			if (SamOps.isHeaderLine(line)) continue;
			line2 = reader.readLine();
		    } else if (readsInputFormat == 5) {  // BAM file (single-end)
			parse_line = sam.parseAlignmentLine_BAM();
			if (parse_line == null) {
			    line = null;
			    continue;
			}
			if (parse_line == SamOps.emptyStringArray) continue;
			if (parse_line[0].startsWith("U")) {
			    read = parse_line[9];
			    quality = mapPhred(parse_line[10]);
			} else if (samWriter != null) {
			    name = parse_line[0].substring(1);
			    while (name.startsWith("@")) name = name.substring(1);
			    mapBamNames(parse_line);
			    synchronized(samWriter) {
				samWriter.println(name + "\t" + parse_line[1] + "\t" + parse_line[2] + "\t" + parse_line[3] + "\t" + parse_line[4] + "\t" + parse_line[5] + "\t" + parse_line[6] + "\t" + parse_line[7] + "\t" + parse_line[8] + "\t" + parse_line[9] + "\t" + mapPhred(parse_line[10]));
				continue;
			    }
			} else continue;
		    } else if (readsInputFormat == 6) {  // BAM file (paired-end)
			parse_line = sam.parseAlignmentLine_BAM();
			if (parse_line == null) {
			    line = null;
			    continue;
			}
			parse_line2 = sam.parseAlignmentLine_BAM();
			if (parse_line2 == null) {
			    line = null;
			    continue;
			}
			if ((parse_line == SamOps.emptyStringArray) || (parse_line2 == SamOps.emptyStringArray)) continue;
			if ((parse_line[0].startsWith("U")) || parse_line2[0].startsWith("U")) {
			    read = parse_line[9];
			    quality = mapPhred(parse_line[10]);
			    read2 = parse_line2[9];
			    quality2 = mapPhred(parse_line2[10]);
			} else if (samWriter != null) {
			    name = parse_line[0].substring(1);
			    while (name.startsWith("@")) name = name.substring(1);
			    name2 = parse_line2[0].substring(1);
			    while (name2.startsWith("@")) name2 = name2.substring(1);
			    mapBamNames(parse_line);
			    mapBamNames(parse_line2);
			    synchronized(samWriter) {
				samWriter.println(name + "\t" + parse_line[1] + "\t" + parse_line[2] + "\t" + parse_line[3] + "\t" + parse_line[4] + "\t" + parse_line[5] + "\t" + parse_line[6] + "\t" + parse_line[7] + "\t" + parse_line[8] + "\t" + parse_line[9] + "\t" + mapPhred(parse_line[10]));
				samWriter.println(name2 + "\t" + parse_line2[1] + "\t" + parse_line2[2] + "\t" + parse_line2[3] + "\t" + parse_line2[4] + "\t" + parse_line2[5] + "\t" + parse_line2[6] + "\t" + parse_line2[7] + "\t" + parse_line2[8] + "\t" + parse_line2[9] + "\t" + mapPhred(parse_line2[10]));
				continue;
			    }
			} else continue;
		    }

		    if (isPairedEnd) {  // Using paired-end reads
			synchronized(reader2) {  // Thread safe reading from file
			    line2 = reader2.readLine();
			    if (line2 == null) {
				output("Error - the two files of mate-pair reads contain different numbers of sequencing reads.");
				continue;
			    }
			    if (readsInputFormat == 0) {  // FASTQ file
				name2 = line;
				read2 = reader2.readLine();
				line2 = reader2.readLine();
				quality2 = reader2.readLine();
			    } else if (readsInputFormat == 1) {  // QSEQ file
				parse_line2 = line2.split("\t");
				name2 = parse_line2[0] + ":" + parse_line2[2] + ":" + parse_line2[3] + ":" + parse_line2[4] + ":" + parse_line2[5] + ":" + parse_line2[6] + ":" + parse_line2[7];
				read2 = parse_line2[8];
				quality2 = parse_line2[9];
			    } else if (readsInputFormat == 2) {  // FASTA file
				name2 = nextName2;
				read2 = "";
				while ((line2 != null) && ((line2.length() == 0) || (line2.charAt(0) != '>'))) {
				    read2 += line2;
				    line2 = reader2.readLine();
				}
				if (line2 == null) nextName2 = "";
				else if (line2.startsWith(">")) nextName2 = line2.substring(1);
				else nextName2 = "";  // Case should never be reached
				if (quality_fasta2.length() != read2.length()) quality_fasta2 = new String(new char[read2.length()]).replace('\0', '(');
				quality2 = quality_fasta2;
			    }
			}
		    }
		}

		if (readsInputFormat == 3) {  // SAM file (single-end)
		    parse_line = sam.parseAlignmentLine_SAM(line);
		    if (parse_line == SamOps.emptyStringArray) continue;
		    if (parse_line[0].startsWith("U")) {
			read = parse_line[9];
			quality = parse_line[10];
		    } else if (samWriter != null) {  // Output to SAM file
			name = parse_line[0].substring(1);
			while (name.startsWith("@")) name = name.substring(1);
			synchronized(samWriter) {
			    samWriter.println(name + "\t" + parse_line[1] + "\t" + parse_line[2] + "\t" + parse_line[3] + "\t" + parse_line[4] + "\t" + parse_line[5] + "\t" + parse_line[6] + "\t" + parse_line[7] + "\t" + parse_line[8] + "\t" + parse_line[9] + "\t" + parse_line[10]);
			}
			continue;
		    } else continue;
		} else if (readsInputFormat == 4) {  // SAM File (paired-end)
		    parse_line = sam.parseAlignmentLine_SAM(line);
		    parse_line2 = sam.parseAlignmentLine_SAM(line2);
		    if ((parse_line == SamOps.emptyStringArray) || (parse_line2 == SamOps.emptyStringArray)) continue;
		    if ((parse_line[0].startsWith("U")) || parse_line2[0].startsWith("U")) {
			read = parse_line[9];
			quality = parse_line[10];
			read2 = parse_line2[9];
			quality2 = parse_line2[10];
		    } else if (samWriter != null) {
			name = parse_line[0].substring(1);
			while (name.startsWith("@")) name = name.substring(1);
			name2 = parse_line2[0].substring(1);
			while (name2.startsWith("@")) name2 = name2.substring(1);
			synchronized(samWriter) {
			    samWriter.println(name + "\t" + parse_line[1] + "\t" + parse_line[2] + "\t" + parse_line[3] + "\t" + parse_line[4] + "\t" + parse_line[5] + "\t" + parse_line[6] + "\t" + parse_line[7] + "\t" + parse_line[8] + "\t" + parse_line[9] + "\t" + parse_line[10]);
			    samWriter.println(name2 + "\t" + parse_line2[1] + "\t" + parse_line2[2] + "\t" + parse_line2[3] + "\t" + parse_line2[4] + "\t" + parse_line2[5] + "\t" + parse_line2[6] + "\t" + parse_line2[7] + "\t" + parse_line2[8] + "\t" + parse_line2[9] + "\t" + parse_line2[10]);
			}
			continue;
		    } else continue;
		}

		// Process read
		read = read.toUpperCase();
		boolean ambiguous = false;
		for (int i=0; i<read.length(); i++) {
		    if ((read.charAt(i) != 'A') && (read.charAt(i) != 'C') &&
			(read.charAt(i) != 'G') && (read.charAt(i) != 'T'))
			ambiguous = true;
		}
		if (isPairedEnd || ((sam != null) && (sam.isPairedEnd()))) {  // Using paired-end reads
		    read2 = read2.toUpperCase();
		    for (int i=0; i<read2.length(); i++) {
			if ((read2.charAt(i) != 'A') && (read2.charAt(i) != 'C') &&
			    (read2.charAt(i) != 'G') && (read2.charAt(i) != 'T'))
			    ambiguous = true;
		    }

		    // Orient paired-end reads (0=ff, 1=fr, 2=rf, 3=rr)
		    if (pairedEndOrientation.charAt(0) == 'r') {  // Reverse first read
			read = reverseComplement(read);
			quality = reverse(quality);
		    }
		    if (pairedEndOrientation.charAt(1) == 'r') {  // Reverse second read
			read2 = reverseComplement(read2);
			quality2 = reverse(quality2);
		    }
		} else {  // Using single-end reads
		    // Orient single-end reads (true=reverse_complement, false=nothing)
		    if (singleEndOrientationReverseComplement) {  // Reverse complement read
			read = reverseComplement(read);
			quality = reverse(quality);
		    }
		}

		// Process quality
		if (quality.length() != read.length()) {
		    output("Quality sequence has different length than read sequence.\n");
		    continue;
		}
		int endIndex = quality.length()-1;
		while ((endIndex >= 0) && ((int)(quality.charAt(endIndex))-PHRED_offset == 2)) endIndex--;
		quality = quality.substring(0, endIndex+1);
		read = read.substring(0, endIndex+1);
		if (qualities.length != read.length()) qualities = new int[read.length()];
		boolean invalid = false;
		for (int i=0; i<quality.length(); i++) {
		    qualities[i] = (int)(quality.charAt(i)) - PHRED_offset;
		    if ((qualities[i] < 0) || (qualities[i] > 93)) invalid = true;
		}
		if (isPairedEnd || ((sam != null) && (sam.isPairedEnd()))) {  // Using paired-end reads
		    if (quality2.length() != read2.length()) {
			output("Quality sequence has different length than read sequence.\n");
			continue;
		    }
		    int endIndex2 = quality2.length()-1;
		    while ((endIndex2 >= 0) && ((int)(quality2.charAt(endIndex2))-PHRED_offset == 2)) endIndex2--;
		    quality2 = quality2.substring(0, endIndex2+1);
		    read2 = read2.substring(0, endIndex2+1);
		    if (qualities2.length != read2.length()) qualities2 = new int[read2.length()];
		    for (int i=0; i<quality2.length(); i++) {
			qualities2[i] = (int)(quality2.charAt(i)) - PHRED_offset;
			if ((qualities2[i] < 0) || (qualities2[i] > 93)) invalid = true;
		    }
		}

		// Align read to genome
		totalReads.getAndIncrement();
		if (ambiguous) ambiguousReads.getAndIncrement();
		else if ((quality.length() < minReadLength) || (isPairedEnd && (quality2.length() < minReadLength)) || ((sam != null) && (sam.isPairedEnd()) && (quality2.length() < minReadLength))) lowQualityReads.getAndIncrement();  // Minimum read length
		else if (invalid || (read.length() < minReadLength)) invalidQualityReads.getAndIncrement();  // Invalid quality score
		else {
		    index.inexactMatch(read, qualities, a, hits, seedHits);
		    if ((isPairedEnd && (hits.size()>0)) || ((sam != null) && (sam.isPairedEnd()) && (hits.size()>0))) {  // Using paired-end reads
			index.inexactMatch(read2, qualities2, a2, hits2, seedHits2);
			Hit.combinePairedEndHits(combinedHits, hits, hits2, maxPairedEndLength);
			// Swap "hits" and "combinedHits"
			ArrayList<Hit> temp = hits;
			hits = combinedHits;
			combinedHits = temp;
		    }
		    for (int i=0; i<hits.size(); i++) {  // Update coordinates
			mapHitsToCoordinates(hits.get(i));
		    }
		    if (hits.size() > 0) {  // Update stats
			if (hits.get(0).getErrors() == 0) exactMappedReads.getAndIncrement(hits.get(0).getRepliconIndex());
			else inexactMappedReads.getAndIncrement(hits.get(0).getRepliconIndex());
			if (hits.size() == 1) mappedOnceReads.getAndIncrement();
			else mappedMoreThanOnceReads.getAndIncrement();
			avgLengthReads.addAndGet(hits.get(0).getLength());
			numMappingReads.incrementAndGet();
		    }
		}

		// Output to SAM file
		if (samWriter != null) {
		    synchronized(samWriter) {
			int flag = 0;
			int flag2 = 0;
			if ((quality.length() < minReadLength) || (isPairedEnd && (quality2.length() < minReadLength)) || ((sam != null) && (sam.isPairedEnd()) && (quality2.length() < minReadLength))) {  // READ TOO SHORT, DON'T USE
			    // Do nothing.
			}
			else if (ambiguous || invalid || (hits.size() == 0)) {  // UNALIGNED
			    if (isPairedEnd || ((sam != null) && (sam.isPairedEnd()))) {  // Paired-end
				if (readsInputFormat == 6) {  // BAM
				    mapBamNames(parse_line);
				    mapBamNames(parse_line2);
				    parse_line[10] = mapPhred(parse_line[10]);
				    parse_line2[10] = mapPhred(parse_line2[10]);
				}
				if ((readsInputFormat == 4) || (readsInputFormat == 6)) {  // SAM/BAM
				    name = parse_line[0].substring(1);
				    while (name.startsWith("@")) name = name.substring(1);
				    name2 = parse_line2[0].substring(1);
				    while (name2.startsWith("@")) name2 = name2.substring(1);
				    samWriter.println(name + "\t" + parse_line[1] + "\t" + parse_line[2] + "\t" + parse_line[3] + "\t" + parse_line[4] + "\t" + parse_line[5] + "\t" + parse_line[6] + "\t" + parse_line[7] + "\t" + parse_line[8] + "\t" + parse_line[9] + "\t" + parse_line[10]);
				    samWriter.println(name2 + "\t" + parse_line2[1] + "\t" + parse_line2[2] + "\t" + parse_line2[3] + "\t" + parse_line2[4] + "\t" + parse_line2[5] + "\t" + parse_line2[6] + "\t" + parse_line2[7] + "\t" + parse_line2[8] + "\t" + parse_line2[9] + "\t" + parse_line2[10]);
				} else {  // FASTA, QSEQ, FASTQ
				    flag = 77;  // 01001101
				    if (pairedEndOrientation.charAt(0) == 'r') {  // Reverse first read
					read = reverseComplement(read);
					quality = reverse(quality);
				    }
				    while (name.startsWith("@")) name = name.substring(1);
				    if (name.length() == 0) {
					name = "r" + nameCount.getAndIncrement();
					name2 = name;
				    }
				    samWriter.println(name + "\t" + flag + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + read + "\t" + quality);
				    flag2 = 141;  // 10001101
				    if (pairedEndOrientation.charAt(1) == 'r') {  // Reverse second read
					read2 = reverseComplement(read2);
					quality2 = reverse(quality2);
				    }
				    while (name2.startsWith("@")) name2 = name2.substring(1);
				    samWriter.println(name2 + "\t" + flag2 + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + read2 + "\t" + quality2);
				}
			    } else {  // Single-end
				if (readsInputFormat == 5) {  // BAM
				    mapBamNames(parse_line);
				    parse_line[10] = mapPhred(parse_line[10]);
				}
				if ((readsInputFormat == 3) || (readsInputFormat == 5)) {  // SAM/BAM
				    name = parse_line[0].substring(1);
				    while (name.startsWith("@")) name = name.substring(1);
				    samWriter.println(name + "\t" + parse_line[1] + "\t" + parse_line[2] + "\t" + parse_line[3] + "\t" + parse_line[4] + "\t" + parse_line[5] + "\t" + parse_line[6] + "\t" + parse_line[7] + "\t" + parse_line[8] + "\t" + parse_line[9] + "\t" + parse_line[10]);
				} else {  // FASTA, QSEQ, FASTQ
				    flag = 4;  // Unmapped
				    if (singleEndOrientationReverseComplement) {  // Reverse complement read
					read = reverseComplement(read);
					quality = reverse(quality);
				    }
				    while (name.startsWith("@")) name = name.substring(1);
				    if (name.length() == 0) name = "r" + nameCount.getAndIncrement();
				    samWriter.println(name + "\t" + flag + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + read + "\t" + quality);
				}
			    }
			} else {  // ALIGNED
			    if (isPairedEnd || ((sam != null) && (sam.isPairedEnd()))) {  // Paired-end (FASTA, QSEQ, FASTQ, SAM, BAM)
				flag = 67;   // 01000011
				flag2 = 131;  // 10000011
				if (pairedEndOrientation.charAt(0) == 'r') {  // Reverse first read
				    read = reverseComplement(read);
				    quality = reverse(quality);
				    if (hits.get(0).getStrand() == '+') {
					flag += 16;
					flag2 += 32;
				    }
				} else {
				    if (hits.get(0).getStrand() == '-') {
					flag += 16;
					flag2 += 32;
				    }
				}
				if (pairedEndOrientation.charAt(1) == 'r') {  // Reverse second read
				    read2 = reverseComplement(read2);
				    quality2 = reverse(quality2);
				    if (hits.get(0).getStrand() == '+') {
					flag2 += 16;
					flag += 32;
				    }
				} else {
				    if (hits.get(0).getStrand() == '+') {
					flag2 += 16;
					flag += 32;
				    }
				}
				if ((readsInputFormat == 4) || (readsInputFormat == 6)) {  // SAM/BAM
				    name = parse_line[0].substring(1);
				    name2 = parse_line2[0].substring(1);
				}
				while (name.startsWith("@")) name = name.substring(1);
				while (name2.startsWith("@")) name2 = name2.substring(1);
				if (name.length() == 0) {
				    name = "r" + nameCount.getAndIncrement();
				    name2 = name;
				}
				samWriter.println(name + "\t" + flag + "\t" + genomeNames[hits.get(0).getRepliconIndex()] + "\t" + hits.get(0).getStart() + "\t" + "255" + "\t" + read.length() + "M" + "\t" + "=" + "\t" + (hits.get(0).getStop()-read2.length()+1) + "\t" + hits.get(0).getLength() + "\t" + read + "\t" + quality);
				samWriter.println(name2 + "\t" + flag2 + "\t" + genomeNames[hits.get(0).getRepliconIndex()] + "\t" + (hits.get(0).getStop()-read2.length()+1) + "\t" + "255" + "\t" + read2.length() + "M" + "\t" + "=" + "\t" + hits.get(0).getStart() + "\t" + (0-hits.get(0).getLength()) + "\t" + read2 + "\t" + quality2);
			    } else {  // Single-end (FASTA, QSEQ, FASTQ, SAM, BAM)
				flag = 0;
				if (singleEndOrientationReverseComplement) {  // Reverse complement read
				    read = reverseComplement(read);
				    quality = reverse(quality);
				    if (hits.get(0).getStrand() == '+') flag += 16;
				} else {
				    if (hits.get(0).getStrand() == '-') flag += 16;
				}
				if ((readsInputFormat == 3) || (readsInputFormat == 5))  // SAM/BAM
				    name = parse_line[0].substring(1);
				while (name.startsWith("@")) name = name.substring(1);
				if (name.length() == 0) name = "r" + nameCount.getAndIncrement();
				samWriter.println(name + "\t" + flag + "\t" + genomeNames[hits.get(0).getRepliconIndex()] + "\t" + hits.get(0).getStart() + "\t" + "255" + "\t" + read.length() + "M" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + read + "\t" + quality);
			    }
			}
		    }
		}
	    }
	} catch (IOException e) {
	    output("Error - could not read in from file " + readsFile + "\n\n");
	    System.exit(0);
	}
    }

    /**
     * Return the base file name for the compressed sequencing reads file.
     */
    public String getCompressedFileName(String readsFile, String genomeFile) {
	if (!isPairedEnd) return getFileNameBase(readsFile) + "_" + getFileNameBase(genomeFile) + parametersString;
	else return getFileNameBase(readsFile) + "_" + getFileNameBase(pairedEndFile) + "_" + getFileNameBase(genomeFile) + parametersString;
    }

    /**
     * Return a list of the base file names for the compressed sequencing reads files.
     */
    public ArrayList<String> getListOfCompressedFileNames(String readsFile) {
	ArrayList<String> fileNames = new ArrayList<String>(numSequences);
	for (int i=0; i<numSequences; i++) {
	    fileNames.add(getCompressedFileName(readsFile, sequenceFiles[i]));
	}
	return fileNames;
    }



    /**********************************************
     **********   Public Class Methods   **********
     **********************************************/

    /**
     * Output mapped reads to WIG file that can be loaded
     * by a genome browser.
     */
    public static void outputBrowserFile(String readsFile, String genomeFile, String dir1, String dir2, String outFileName, int[] coordinates_plus, int[] coordinates_minus, int trackRange) {

	try {
	    // Set up directory
	    File d = new File(dir1);
	    if (!d.isDirectory()) d.mkdir();
	    d = new File(dir1 + dir2);
	    if (!d.isDirectory()) d.mkdir();
	    PrintWriter out = new PrintWriter(new File(dir1 + dir2 + outFileName + ".plus" + ".wig"));
	    out.println("track name=" + getFileNameBase(readsFile) + "_PLUS color=0,0,255 graphType=bar viewLimits=0:" + trackRange);
	    out.println("fixedStep chrom=" + getGenomeName(genomeFile) + " start=1 step=1");
	    for (int j=1; j<coordinates_plus.length; j++) out.println(coordinates_plus[j]);
	    out.close();
	    out = new PrintWriter(new File(dir1 + dir2 + outFileName + ".minus" + ".wig"));
	    out.println("track name=" + getFileNameBase(readsFile) + "_MINUS color=255,0,0 graphType=bar viewLimits=0:" + trackRange);
	    out.println("fixedStep chrom=" + getGenomeName(genomeFile) + " start=1 step=1");
	    for (int j=1; j<coordinates_minus.length; j++) out.println(coordinates_minus[j]);
	    out.close();
	} catch (IOException e) {
	    output("Could not write wiggle file in directory " + (dir1 + dir2) + "\n");
	}
    }

    public static void output(String s) {
	if (summaryWriter != null) summaryWriter.print(s);
	if (output == null) System.out.print(s);
	else output.append(s);
    }



    /**************************************************
     **********   Private Instance Methods   **********
     **************************************************/

    /**
     * Read in reads file and map all reads to reference
     * index using multiple threads.
     */
    private void mapReads() {

        try {
	    if (isCompressedFile) reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readsFile))));
	    else reader = new BufferedReader(new FileReader(readsFile));
	    if (isPairedEnd) {
		if (FileOps.isGZIP(pairedEndFile)) reader2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(pairedEndFile))));
		else reader2 = new BufferedReader(new FileReader(pairedEndFile));
	    }
	    Thread[] threads = new Thread[numThreads];
	    for (int i=0; i<numThreads; i++) {
		threads[i] = new ProcessReads_Thread(this);
		threads[i].start();
	    }
	    for (int i=0; i<numThreads; i++) threads[i].join();
	    reader.close();
	    if (isPairedEnd) reader2.close();
	} catch (IOException e) {
	    output("Error - could not read from " + readsFile + "\n\n");
	    System.exit(0);
	} catch (InterruptedException e) {
	    output("Error - thread execution was interrupted.\n\n");
	    System.exit(0);
	}
    }

    /**
     * Output the reads aligning to each replicon to a file.
     */
    private void outputReadsToFile() {
	for (int i=0; i<numSequences; i++) {
	    String outFileName = getCompressedFileName(readsFile, sequenceFiles[i]);
	    FileOps.writeCompressedFile(outFileName, readsFile, sequenceFiles[i], index.getReplicon(i).getName(), index.getReplicon(i).getLength(), totalReads.get(), invalidQualityReads.get(), ambiguousReads.get(), lowQualityReads.get(), mappedOnceReads.get(), mappedMoreThanOnceReads.get(), exactMappedReads.get(i), inexactMappedReads.get(i), avgLengthReads.get(), toArray(coordinates_plus[i]), toArray(coordinates_minus[i]));

	    // Write Wiggle files
	    if (Peregrine.outputBrowserFile) {
		int[] plusCoords = new int[coordinates_plus[i].length()];
		int[] minusCoords = new int[coordinates_minus[i].length()];
		for (int j=0; j<plusCoords.length; j++) plusCoords[j] = coordinates_plus[i].get(j);
		for (int j=0; j<minusCoords.length; j++) minusCoords[j] = coordinates_minus[i].get(j);
		int lengthOfAllReplicons = 0;
		for (int k=0; k<index.getNumReplicons(); k++) lengthOfAllReplicons += index.getReplicon(k).getLength();
		long trackRange = 2 * avgLengthReads.get() * (long)(exactMappedReads.get(i) + inexactMappedReads.get(i)) / lengthOfAllReplicons;
		outputBrowserFile(readsFile, sequenceFiles[i], outputDIR, browserDIR, outFileName, plusCoords, minusCoords, (int)trackRange);
	    }
	}
    }

    /**
     * Outputs to STDOUT stats on how many reads align to the genomic sequence.
     */
    private void outputResults() {
	output("Total reads:            \t" + totalReads + "\n");
	//output("Invalid quality reads:  \t" + invalidQualityReads + "\t" + Math.round((100.0*invalidQualityReads.get())/totalReads.get()) + "%" + "\n");
	//output("Ambiguous reads:        \t" + ambiguousReads + "\t" + Math.round((100.0*ambiguousReads.get())/totalReads.get()) + "%" + "\n");
	//output("Low quality reads:      \t" + lowQualityReads + "\t" + Math.round((100.0*lowQualityReads.get())/totalReads.get()) + "%" + "\n");
	//output("Aligned 1 time reads:   \t" + mappedOnceReads + "\t" + Math.round((100.0*mappedOnceReads.get())/totalReads.get()) + "%" + "\n");
	//output("Aligned >1 time reads:   \t" + mappedMoreThanOnceReads + "\t" + Math.round((100.0*mappedMoreThanOnceReads.get())/totalReads.get()) + "%" + "\n");
	//output("Exact alignment reads:  \t" + exactMappedReads + "\t" + Math.round((100.0*exactMappedReads.get())/totalReads.get(i)) + "%" + "\n");
	//output("Inexact alignment reads:\t" + inexactMappedReads + "\t" + Math.round((100.0*inexactMappedReads.get(i))/totalReads.get()) + "%" + "\n");
	//output("\n");
	for (int i=0; i<Peregrine.numSequences; i++) {
	    output("Successfully aligned reads:\t" + (exactMappedReads.get(i) + inexactMappedReads.get(i)) + "\t" + Math.round((100.0*(exactMappedReads.get(i)+inexactMappedReads.get(i)))/totalReads.get()) + "%" + "\t" + "(" + getInformalName(sequenceFiles[i]) + ")" + "\n");
	    outputAnnotationStats(i);  // Number of reads aligning genes and RNAs
	}
	output("\n");

	if (samWriter != null) output("Sequencing reads written to SAM file:\t" + outputDIR + getFileNameBase(readsFile) + ".sam" + "\n\n");

	if (time) {  // Output time taken to execute program
	    output("Time taken to align reads:\t" + (runningTime/60000) + " minutes " + ((runningTime%60000)/1000) + " seconds\n");
	}
    }

    /**
     * Determines the file format of the input reads file.
     * Returns 0 if the file is in FASTQ format.
     * Returns 1 if the file is in QSEQ format.
     * Returns 2 if the file is in FASTA format.
     * Returns 3 if the file is in SAM single-end format.
     * Returns 4 if the file is in SAM paired-end format.
     * Returns 5 if the file is in BAM single-end format.
     * Returns 6 if the file is in BAM paired-end format.
     */
    private int getReadsInputFileType(String readsFile) {
	int MAX_LINES = 10000;
	int lineCounter = 0;
	int fastqCount = 0;
	int qseqCount = 0;
	int fastaCount = 0;
	try {
	    // First check if we have a SAM file
	    sam = new SamOps(readsFile, MAX_LINES, false);
	    if (sam.isValidSam() || sam.isValidBam()) {
		String[] replicons = sam.getReplicons();
		output("\nUser should make sure that reference replicons input by the user correspond to those in SAM/BAM file.\n");
		output("\tReference replicons input by user:\n");
		for (int i=0; i<numSequences; i++) output("\t\t" + getGenomeName(sequenceFiles[i]) + " with length " + (index.getReplicon(i).getLength()-1) + " nts\n");
		output("\tReference replicons found in SAM file:\n");
		for (int i=0; i<replicons.length; i++) output("\t\t" + replicons[i] + " with length " + sam.getRepliconLength(i) + " nts\n");
		output("\n");
		if (sam.isValidSam() && !sam.isPairedEnd()) return 3;  // Single-end SAM file
		else if (sam.isValidSam()) {  // Paired-end SAM file
		    Peregrine.stopAfterOneHit = false;
		    return 4;
		} else if (sam.isValidBam() && !sam.isPairedEnd()) return 5;  // Single-end BAM file
		else if (sam.isValidBam()) {  // Paired-end BAM file
		    Peregrine.stopAfterOneHit = false;
		    return 6;
		} else return -1;  // This case should never be reached
	    } else sam = null;

	    Scanner reader = null;
	    isCompressedFile = FileOps.isGZIP(readsFile);
	    if (isCompressedFile) reader = new Scanner(new GZIPInputStream(new FileInputStream(readsFile)));
	    else reader = new Scanner(new File(readsFile));
	    String line;
	    while (reader.hasNextLine() && (lineCounter < MAX_LINES)) {
		line = reader.nextLine();
		if (line.length() > 0) {
		    if (line.charAt(0) == '@') fastqCount++;
		    else if (line.split("\t").length >= 11) qseqCount++;
		    else if (line.charAt(0) == '>') fastaCount++;
		}
		lineCounter++;
	    }
	    reader.close();
	    if (fastqCount > Math.max(qseqCount, fastaCount)) return 0;
	    else if (qseqCount > fastaCount) return 1;
	    else if (fastaCount > 0) return 2;
	    else {
		output("Error - could not recognize format of file " + readsFile + "\n\n");
		System.exit(0);
	    }
	} catch (IOException e) {
	    output("\nError - could not open file " + readsFile + "\n\n");
	    System.exit(0);
	}
	return -1;
    }

    /**
     * Depending on version of Solexa machine used, the quality
     * score may range from 35 (or 37) up to 73, or it may range
     * from 66 up to 104. To obtain the Phred score, we subtract 
     * either 33 or 64, depending on the version. Here, we get the
     * maximum quality score over all nucleotides in the first
     * 100,000 reads (it should be either 73 or 104). The offset 
     * is the number we must subtract from this maximum to 
     * obtain 40, i.e., either 33 or 64.
     */
    private int getPhredOffset(String readsFile) {
	int lineCounter = 0;
	try {
	    int maxQuality = Integer.MIN_VALUE;
	    if ((readsInputFormat == 5) || (readsInputFormat == 6)) {  // BAM file
		maxQuality = sam.getMaxQuality();
	    } else {
		Scanner reader = null;
		if (isCompressedFile) reader = new Scanner(new GZIPInputStream(new FileInputStream(readsFile)));
		else reader = new Scanner(new File(readsFile));
		String line;
		String read = "";
		String quality = "";

		while (reader.hasNextLine() && (lineCounter < 100000)) {

		    line = reader.nextLine();
		    if (readsInputFormat == 0) {  // FASTQ file
			read = reader.nextLine();
			line = reader.nextLine();
			quality = reader.nextLine();
		    } else if (readsInputFormat == 1) {  // QSEQ file
			String[] parse_line = line.split("\t");
			read = parse_line[8];
			quality = parse_line[9];
		    } else if (readsInputFormat == 2) {  // FASTA file
			read = reader.nextLine();
			quality = "(";  // ASCII value of 40
		    } else if ((readsInputFormat == 3) || (readsInputFormat == 4)) {  // SAM file
			if (line.startsWith("@")) {
			    read = "";
			    quality = "";
			} else {
			    String[] parse_line = line.split("\t");
			    read = parse_line[9];
			    quality = parse_line[10];
			    if (quality.equals("*")) quality = "";
			}
		    }

		    // Process quality
		    for (int i=0; i<quality.length(); i++)
			maxQuality = Math.max(maxQuality, (int)(quality.charAt(i)));

		    lineCounter++;
		}
		reader.close();
	    }

	    // Are we Phred+33 or Phred+64?
	    if (Math.abs(maxQuality - 40 - 33) < Math.abs(maxQuality - 40 - 64))
		return 33;
	    else 
		return 64;
	} catch (IOException e) {
	    output("\nError - could not open file " + readsFile + "\n\n");
	    System.exit(0);
	}
	return -1;
    }

    /**
     * Given the hit of a read mapping to the genome, we update
     * the coordinates corresponding to the read hit.
     */
    private void mapHitsToCoordinates(Hit a) {
	AtomicIntegerArray coordinates = this.coordinates_plus[a.getRepliconIndex()];  // Plus strand
	if (a.getStrand() == '-') coordinates = this.coordinates_minus[a.getRepliconIndex()];  // Minus strand
	for (int i=a.getStart(); i<=a.getStop(); i++)
	    coordinates.getAndIncrement(i);
    }

    /**
     * Output number of reads aligning sense/antisense to
     * annotated genes and RNAs in replicon i.
     */
    private void outputAnnotationStats(int i) {
	if ((annotationsPlus == null) || annotationsMinus == null) return;
	if ((i >= annotationsPlus.length) || (i >= annotationsMinus.length)) return;

	long gene_sense = 0;
	long gene_antisense = 0;
	long rRNA_sense = 0;
	long rRNA_antisense = 0;
	long tRNA_sense = 0;
	long tRNA_antisense = 0;
	long miscRNA_sense = 0;
	long miscRNA_antisense = 0;
	long IG = 0;

	String[] annotationPlus = annotationsPlus[i];
	String[] annotationMinus = annotationsMinus[i];
	AtomicIntegerArray coordinate_plus = coordinates_plus[i];
	AtomicIntegerArray coordinate_minus = coordinates_minus[i];
	for (int j=1; j<coordinate_plus.length(); j++) {  // 1-indexed
	    // There are five possibilities on each strand (rRNA, tRNA, RNA, gene, empty).
	    // We handle all twenty-five cases.
	    if (annotationPlus[j].equals("rRNA")) {  // Plus rRNA
		rRNA_sense += coordinate_plus.get(j);
		if (annotationMinus[j].length() == 0) rRNA_antisense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("rRNA")) rRNA_sense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("tRNA")) tRNA_sense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("RNA")) miscRNA_sense += coordinate_minus.get(j);
		else gene_sense += coordinate_minus.get(j);
	    } else if (annotationPlus[j].equals("tRNA")) {  // Plus tRNA
		tRNA_sense += coordinate_plus.get(j);
		if (annotationMinus[j].length() == 0) tRNA_antisense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("rRNA")) rRNA_sense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("tRNA")) tRNA_sense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("RNA")) miscRNA_sense += coordinate_minus.get(j);
		else gene_sense += coordinate_minus.get(j);
	    } else if (annotationPlus[j].equals("RNA")) {  // Plus misc RNA
		miscRNA_sense += coordinate_plus.get(j);
		if (annotationMinus[j].length() == 0) miscRNA_antisense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("rRNA")) rRNA_sense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("tRNA")) tRNA_sense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("RNA")) miscRNA_sense += coordinate_minus.get(j);
		else gene_sense += coordinate_minus.get(j);
	    } else if (annotationPlus[j].length() > 0) {  // Plus protein-coding gene
		gene_sense += coordinate_plus.get(j);
		if (annotationMinus[j].length() == 0) gene_antisense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("rRNA")) rRNA_sense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("tRNA")) tRNA_sense += coordinate_minus.get(j);
		else if (annotationMinus[j].equals("RNA")) miscRNA_sense += coordinate_minus.get(j);
		else gene_sense += coordinate_minus.get(j);
	    } else {  // Plus no annotation
		if (annotationMinus[j].equals("rRNA")) {
		    rRNA_sense += coordinate_minus.get(j);
		    rRNA_antisense += coordinate_plus.get(j);
		} else if (annotationMinus[j].equals("tRNA")) {
		    tRNA_sense += coordinate_minus.get(j);
		    tRNA_antisense += coordinate_plus.get(j);
		} else if (annotationMinus[j].equals("RNA")) {
		    miscRNA_sense += coordinate_minus.get(j);
		    miscRNA_antisense += coordinate_plus.get(j);
		} else if (annotationMinus[j].length() > 0) {
		    gene_sense += coordinate_minus.get(j);
		    gene_antisense += coordinate_plus.get(j);
		} else {
		    IG += coordinate_minus.get(j);
		    IG += coordinate_plus.get(j);
		}
	    }
	}

	// Output stats
	long total = gene_sense + gene_antisense + rRNA_sense + rRNA_antisense + tRNA_sense + tRNA_antisense + miscRNA_sense + miscRNA_antisense + IG;
	output("\t" + "Aligning (sense) to protein-coding genes:\t" + Math.round(100.0*gene_sense/total) + "%" + "\n");
	output("\t" + "Aligning (antisense) to protein-coding genes:\t" + Math.round(100.0*gene_antisense/total) + "%" + "\n");
	output("\t" + "Aligning (sense) to ribosomal RNAs:\t" + Math.round(100.0*rRNA_sense/total) + "%" + "\n");
	output("\t" + "Aligning (antisense) to ribosomal RNAs:\t" + Math.round(100.0*rRNA_antisense/total) + "%" + "\n");
	output("\t" + "Aligning (sense) to transfer RNAs:\t" + Math.round(100.0*tRNA_sense/total) + "%" + "\n");
	output("\t" + "Aligning (antisense) to transfer RNAs:\t" + Math.round(100.0*tRNA_antisense/total) + "%" + "\n");
	output("\t" + "Aligning (sense) to miscellaneous RNAs:\t" + Math.round(100.0*miscRNA_sense/total) + "%" + "\n");
	output("\t" + "Aligning (antisense) to miscellaneous RNAs:\t" + Math.round(100.0*miscRNA_antisense/total) + "%" + "\n");
	output("\t" + "Aligning to unannotated regions:        \t" + Math.round(100.0*IG/total) + "%" + "\n");
    }

    /**
     * Open SAM file for writing.
     * Output header information to SAM file.
     */
    private void outputSamHeader() {
	String samFileName = outputDIR + getFileNameBase(readsFile) + ".sam";
	if (readsFile.equals(samFileName)) return;  // Do not read and write same file at same time
	try {
	    samWriter = new PrintWriter(new File(samFileName));
	    samWriter.println("@HD" + "\t" + "VN:1.0" + "\t" + "SO:unsorted");
	    genomeNames = new String[sequenceFiles.length];
	    for (int i=0; i<sequenceFiles.length; i++) {
		samWriter.println("@SQ" + "\t" + "SN:" + getGenomeName(sequenceFiles[i]) + "\t" + "LN:" + (index.getReplicon(i).getLength()-1) + "\t" + "SP:" + getInformalName(sequenceFiles[i]));
		genomeNames[i] = getGenomeName(sequenceFiles[i]);
	    }
	    samWriter.println("@PG" + "\t" + "ID:Rockhopper" + "\t" + "PN:Rockhopper" + "\t" + "VN:" + Rockhopper.version);
	} catch (IOException f) {
	    samWriter = null;
	}
    }

    /**
     * Takes a String array representing one line from a BAM file and
     * we map the repliconIndex (at parse_line[2] and parse_line[6])
     * to the corresponding replicon name.
     */
    private void mapBamNames(String[] parse_line) {
	if (parse_line[2].equals("*")) parse_line[6] = "*";
	else {
	    int repliconIndex = 0;
	    try {
		repliconIndex = Integer.parseInt(parse_line[2]);
	    } catch (NumberFormatException e) {
		parse_line[2] = "*";
		parse_line[6] = "*";
		return;
	    }
	    if (repliconIndex < 0) {
		parse_line[2] = "*";
		parse_line[6] = "*";
	    } else parse_line[2] = genomeNames[repliconIndex];
	    if (!parse_line[6].equals("*")) {
		if (parse_line[2].equals(parse_line[6])) parse_line[6] = "=";
		else {
		    try {
			repliconIndex = Integer.parseInt(parse_line[6]);
		    } catch (NumberFormatException e) {
			parse_line[6] = "*";
			return;
		    }
		    if (repliconIndex < 0) parse_line[6] = "*";
		    else parse_line[6] = genomeNames[repliconIndex];
		}
	    }
	}
    }

    /**
     * Return a String representation of the parameter settings.
     */
    private String parametersToString() {
	return "_v" + numSequences + "_m" + (int)(100*percentMismatches) + "_a" + ("" + stopAfterOneHit).toUpperCase().charAt(0) + "_d" + maxPairedEndLength + "_l" + (int)(100*percentSeedLength) + "_" + pairedEndOrientation + "_c" + ("" + singleEndOrientationReverseComplement).toUpperCase().charAt(0);
    }

    /**
     * Returns an int[] version of an atomic integer array.
     */
    private int[] toArray(AtomicIntegerArray a) {
	int[] b = new int[a.length()];
	for (int i=0; i<a.length(); i++) b[i] = a.get(i);
	return b;
    }

    /**
     * Returns a reversed version of String s.
     */
    private String reverse(String s) {
	return (new StringBuilder(s)).reverse().toString();
    }

    /**
     * Returns a reverse complemented version of String s.
     */
    private String reverseComplement(String s) {
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



    /***********************************************
     **********   Private Class Methods   **********
     ***********************************************/

    /**
     * Free up memory by deleting class variables that
     * we are finished with.
     */
    public static void releaseMemory() {
	annotationsPlus = null;
	annotationsMinus = null;
    }



    /***********************************************
     **********   Private Class Methods   **********
     ***********************************************/

    /**
     * Returns the base of the specified file name, i.e., 
     * the file extension and path to the file are excluded.
     */
    private static String getFileNameBase(String fileName) {
	int slashIndex = fileName.lastIndexOf(File.separator);
	if (slashIndex == -1) slashIndex = fileName.lastIndexOf("/");
	int periodIndex = fileName.indexOf('.', slashIndex+1);
	if (periodIndex == -1) periodIndex = fileName.length();
	return fileName.substring(slashIndex+1, periodIndex);
    }

    /**
     * Returns the name of the replicon from the FASTA header line,
     * e.g., gi|59800473|ref|NC_002946.2|
     */ 
    private static String getGenomeName(String fileName) {
	String genomeName = "";
	try {
	    Scanner reader = new Scanner(new File(fileName));
	    String line = reader.nextLine();
	    String[] parse_line = line.split("\\s+");
	    genomeName = parse_line[0].substring(1);
	    reader.close();
	} catch (IOException e) {
	    output("Error - could not read in from file " + fileName + "\n\n");
	    System.exit(0);
	}
	return genomeName;
    }

    /**
     * Returns the informal name of the replicon from the FASTA header line,
     * e.g., Escherichia coli str. K-12 substr. MG1655 chromosome
     */
    private static String getInformalName(String fileName) {
	String genomeName = "";
	try {
	    Scanner reader = new Scanner(new File(fileName));
	    String line = reader.nextLine();
	    String[] parse_line = line.split("\\|");
	    if (parse_line.length > 4) {  // Appropriately formatted FASTA header
		String[] parse_name = parse_line[4].split(",");
		genomeName = parse_name[0].trim();
	    } else {  // Unexpected format for FASTA header
		genomeName = line.trim();
	    }
	    reader.close();
	} catch (IOException e) {
	    output("Error - could not read in from file " + fileName + "\n\n");
	    System.exit(0);
	}
	return genomeName;
    }

    /**
     * Checks if quality scores specified by parameter String
     * need to be adjusted by a Phred offset of 33. If not,
     * the original String is returned. If so, a new
     * adjusted String is returned.
     */
    private static String mapPhred(String qualities) {
	int offset = 33;
	int maxPhred = 104;
	boolean needsAdjusting = false;
	boolean needsResetting = false;
	for (int i=0; i<qualities.length(); i++) {
	    if ((int)(qualities.charAt(i)) < offset) needsAdjusting = true;
	    if ((int)(qualities.charAt(i)) > maxPhred) needsResetting = true;
	}
	if (!needsAdjusting && !needsResetting) return qualities;
	StringBuilder sb = new StringBuilder(qualities.length());
	for (int i=0; i<qualities.length(); i++) {
	    if (needsAdjusting) sb.append((char)(offset + qualities.charAt(i)));
	    else if (needsResetting) sb.append('(');
	}
	return sb.toString();
    }

    private static void commandLineArguments(String[] args) {
	if (args.length < 2) {
            System.err.println("\nThe Peregrine application takes one or more files of genomic sequences (such as FASTA genome files) and a file of RNA-seq reads and it aligns the reads to the genomic sequences.");
	    System.err.println("\nThe Peregrine application has two required command line arguments.");
            System.err.println("\nREQUIRED ARGUMENTS\n");
            System.err.println("\t" + "a comma-separated list of files containing replicon sequences in FASTA format");
            System.err.println("\t" + "a file <FASTQ,QSEQ,FASTA> of RNA-seq reads");

            System.err.println("\nOPTIONAL ARGUMENTS\n");
	    System.err.println("\t" + "-2 <file>" + "\t" + "a file <FASTQ,QSEQ,FASTA> of mate-paired RNA-seq reads, (default is single-end rather than paired-end reads)");
	    System.err.println("\t" + "-c <boolean>" + "\t" + "reverse complement single-end reads (default is false)");
	    System.err.println("\t" + "-ff/fr/rf/rr" + "\t" + "orientation of two mate reads for each paired-end read, f=forward and r=reverse_complement (default is fr)");
            System.err.println("\t" + "-m <number>" + "\t" + "allowed mismatches as percent of read length (default is 0.15)");
            System.err.println("\t" + "-a <boolean>" + "\t" + "report 1 alignment (true) or report all optimal alignments (false), (default is true)");
            System.err.println("\t" + "-d <integer>" + "\t" + "maximum number of bases between mate pairs for paired-end reads (default is 500)");
            System.err.println("\t" + "-l <number>" + "\t" + "minimum seed as percent of read length (default is 0.33)");
            System.err.println("\t" + "-p <integer>" + "\t" + "number of processors (default is self-identification of processors)");
	    System.err.println("\t" + "-o <DIR>" + "\t" + "directory where output files are written (default is Peregrine_Results/)");
            System.err.println("\t" + "-time      " + "\t" + "output time taken to execute program");

            System.err.println("\njava Peregrine <options> NC_002505.fna,NC_002506.fna reads.fastq\n");
            System.exit(0);
	}

	// Initially set the number of threads
	Peregrine.numThreads = Runtime.getRuntime().availableProcessors();
	if (Peregrine.numThreads > 4) Peregrine.numThreads *= 0.75;

	int i=0;
	while (i < args.length - 2) {
	    if (args[i].startsWith("-2")) {
                if ((i == args.length-3) || (args[i+1].startsWith("-"))) {
                    System.err.println("Error - command line argument -2 must be followed by a mate-paired sequencing read file.");
                    System.exit(0);
                }
		Peregrine.isPairedEnd = true;
		Peregrine.pairedEndFile = args[i+1];
                i += 2;
	    } else if (args[i].startsWith("-c")) {
                if ((i == args.length-3) || (args[i+1].startsWith("-"))) {
                    System.err.println("Error - command line argument -c must be followed by a boolean.");
                    System.exit(0);
                }
		Peregrine.singleEndOrientationReverseComplement = Boolean.valueOf(args[i+1]);
                i += 2;
	    } else if (args[i].equals("-ff")) {
		Peregrine.pairedEndOrientation = "ff";
                i += 1;
	    } else if (args[i].equals("-fr")) {
		Peregrine.pairedEndOrientation = "fr";
                i += 1;
	    } else if (args[i].equals("-rf")) {
		Peregrine.pairedEndOrientation = "rf";
                i += 1;
	    } else if (args[i].equals("-rr")) {
		Peregrine.pairedEndOrientation = "rr";
                i += 1;
	    } else if (args[i].startsWith("-m")) {
                if ((i == args.length-3) || (args[i+1].startsWith("-"))) {
                    System.err.println("Error - command line argument -m must be followed by a decimal number.");
                    System.exit(0);
                }
		Peregrine.percentMismatches = Double.parseDouble(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-a")) {
                if ((i == args.length-3) || (args[i+1].startsWith("-"))) {
                    System.err.println("Error - command line argument -a must be followed by a boolean.");
                    System.exit(0);
                }
		Peregrine.stopAfterOneHit = Boolean.valueOf(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-d")) {
                if ((i == args.length-3) || (args[i+1].startsWith("-"))) {
                    System.err.println("Error - command line argument -d must be followed by an integer.");
                    System.exit(0);
                }
		Peregrine.maxPairedEndLength = Integer.parseInt(args[i+1]);
		SamOps.maxPairedEndLength = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-l")) {
                if ((i == args.length-3) || (args[i+1].startsWith("-"))) {
                    System.err.println("Error - command line argument -l must be followed by a decimal number.");
                    System.exit(0);
                }
		Peregrine.percentSeedLength = Double.parseDouble(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-p")) {
                if ((i == args.length-3) || (args[i+1].startsWith("-"))) {
                    System.err.println("Error - command line argument -p must be followed by an integer.");
                    System.exit(0);
                }
		Peregrine.numThreads = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-o")) {
                if ((i == args.length-3) || (args[i+1].startsWith("-"))) {
                    System.err.println("Error - command line argument -o must be followed by a directory path.");
                    System.exit(0);
                }
		Peregrine.outputDIR = args[i+1];
		if (Peregrine.outputDIR.charAt(outputDIR.length()-1) != '/') outputDIR += '/';
                i += 2;
	    } else if (args[i].startsWith("-r")) {
                if ((i == args.length-3) || (args[i+1].startsWith("-"))) {
                    System.err.println("Error - command line argument -r must be followed by a long integer.");
                    System.exit(0);
                }
		Peregrine.randomSeed = Long.parseLong(args[i+1]);
                i += 2;
	    } else if (args[i].equals("-time")) {
		Peregrine.time = true;
                i += 1;
	    } else {
		System.err.println("Error - could not recognize command line argument: " + args[i]);
		System.exit(0);
	    }
	}

	// Determine required input files
	Peregrine.sequenceFiles = args[i].split(",");
	Peregrine.readsFile = args[i+1];
	Peregrine.numSequences = Peregrine.sequenceFiles.length;

	// Handle erroneous command line arguments
	for (int j=0; j<sequenceFiles.length; j++) {
	    if (!(new File(sequenceFiles[j]).exists())) {
		System.err.println("Error - file " + sequenceFiles[j] + " does not exist.");
		System.exit(0);
	    }
        }
	if (!(new File(readsFile).exists())) {
	    System.err.println("Error - file " + readsFile + " does not exist.");
	    System.exit(0);
	}
	if (isPairedEnd && !(new File(pairedEndFile).exists())) {
	    System.err.println("Error - file " + pairedEndFile + " does not exist.");
	    System.exit(0);
	}
    }



    /*************************************
     **********   Main Method   **********
     *************************************/

    public static void main(String[] args) {

	commandLineArguments(args);
	Peregrine p = new Peregrine();

    }

}



/*********************************************
 **********   ProcessReads_Thread   **********
 *********************************************/

class ProcessReads_Thread extends Thread {

    private Peregrine p;

    public ProcessReads_Thread(Peregrine p) {
	this.p = p;
    }

    public void run() {
	p.processReads();
    }
}
