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
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.util.zip.GZIPInputStream;
import java.util.concurrent.atomic.AtomicInteger;
import javax.swing.JTextArea;

public class Assembler {

    /********************************************
     **********   Instance Variables   **********
     ********************************************/

    private Dictionary kDict;
    private int PHRED_offset;
    private BufferedReader reader;
    private BufferedReader reader2;  // Mate-pair file if using paired-end reads
    private PrintWriter samWriter;  // Output reads to SAM file
    private boolean isCompressedFile;  // True if gzip file. False otherwise.
    private SamOps sam;  // Used if reads file is in SAM/BAM format
    private AtomicInteger nameCount;  // Used when outputting SAM file and reads have no name
    private boolean firstLineDone;  // Used when reading in the first line of FASTA files
    private long totalTime;  // Time taken to read in sequencing reads and assemble transcripts
    private boolean mapFullReadsToTranscriptome;  // First pass (false), second pass (true).
    private DeNovoTranscripts transcripts;

    // Results
    public AtomicInteger totalReads;
    private AtomicInteger invalidQualityReads;      // Reads containing an invalid quality score
    private AtomicInteger ambiguousReads;           // Reads with ambiguous nucleotides
    private AtomicInteger lowQualityReads;          // Reads with quality scores too low for mapping



    /*****************************************
     **********   Class Variables   **********
     *****************************************/

    // Parameters
    public static JTextArea output;  // Output to GUI or, if null, to System.out
    public static int k = 25;  // Size of k-mer
    private static int readsInputFormat;  // 0=FASTQ, 1=QSEQ, 2=FASTA, 3=SAM(single), 4=SAM(paired), 5=BAM(single), 6=BAM(paired)
    public static ArrayList<String> conditionFiles = new ArrayList<String>();
    public static boolean singleEndOrientationReverseComplement = false;  // RevComp reads
    public static int numThreads = 1;  // Number of threads
    public static int minReadLength = 35;  // Minimum allowed sequencing read length
    public static boolean stopAfterOneHit = true;
    public static boolean computeExpression = true;  // PARAM: compute differential expression
    public static int CAPACITY_POWER = 25;  // PARAM: capacity of hashtable is ~2^n
    public static int MIN_READS_MAPPING = 20;  // PARAM: Min # of full length reads mapping to do novo transcripts
    public static int minTranscriptLength;  // PARAM: Min transcript length: 2*k
    public static int minSeedExpression = 50;  // PARAM: Min expression to initiate a transcript
    public static int minExpression = 5;  // PARAM: Min number of occurrences to be part of transcript
    public static boolean unstranded = false;  // PARAM: Is RNA-seq data strand specific or ambiguous?
    public static String pairedEndOrientation = "fr";  // PARAM: ff, fr, rf, rr
    public static boolean isPairedEnd = false;  // Do we have paired-end reads?
    private static String pairedEndFile;  // Name of mate-pair file if using paired-end reads
    public static int maxPairedEndLength = 500;
    public static String[] labels;
    public static int readsFileIndex;  // Keep track of current sequencing reads file
    public static PrintWriter summaryWriter = null;  // For outputting summary file
    private static String summaryFile = "summary.txt";
    public static String expressionFile = "transcripts.txt";
    public static boolean outputSAM = false;  // Output reads to SAM file
    public static String output_DIR = "Assembler_Results/";  // PARAM: write output files to this directory
    public static boolean time = false;  // PARAM: output execution time of program
    public static boolean verbose = false;  // PARAM: verbose output



    /**************************************
     **********   Constructors   **********
     **************************************/

    public Assembler() {

	totalTime = System.currentTimeMillis();
	FileOps foo = new FileOps(getCompressedFileName(), conditionFiles);
	if (foo.isValid()) {
	    output("\n");
	    for (int i=0; i<conditionFiles.size(); i++) {
		String[] files = conditionFiles.get(i).split(",");
		for (int j=0; j<files.length; j++) {
		    String[] parse_readFiles = files[j].split("%");
		    if (parse_readFiles.length == 1) {  // Using single-end reads (except maybe SAM/BAM)
			String file1 = parse_readFiles[0];
			output("Assembling transcripts from reads in file:             \t" + removePath(file1) + "\n");
			Rockhopper.updateProgress();
			Rockhopper.updateProgress();
			Assembler.isPairedEnd = false;
		    } else {  // Using paired-end reads
			String file1 = parse_readFiles[0];
			String file2 = parse_readFiles[1];
			output("Assembling transcripts from reads in files:\n");
			output("\t" + removePath(file1) + "\n");
			output("\t" + removePath(file2) + "\n");
			Rockhopper.updateProgress();
			Rockhopper.updateProgress();
			Assembler.isPairedEnd = true;
		    }
		}
	    }
	    transcripts = new DeNovoTranscripts(foo.getTranscripts());
	    readsFileIndex = 0;
	    output("\n");
	    for (int i=0; i<conditionFiles.size(); i++) {
		String[] files = conditionFiles.get(i).split(",");
		for (int j=0; j<files.length; j++) {
		    String[] parse_readFiles = files[j].split("%");
		    if (parse_readFiles.length == 1) {  // Using single-end reads (except maybe SAM/BAM)
			String file1 = parse_readFiles[0];
			output("Aligning reads to assembled transcripts using file:\t" + removePath(file1) + "\n");
			Assembler.isPairedEnd = false;
			output("\n\tTotal reads in file:        \t" + DeNovoTranscripts.numReads[readsFileIndex] + "\n");
			output("\tPerfectly aligned reads:    \t" + DeNovoTranscripts.numMappingReads[readsFileIndex] + "\t" + Math.round((100.0*DeNovoTranscripts.numMappingReads[readsFileIndex])/DeNovoTranscripts.numReads[readsFileIndex]) + "%" + "\n\n");
			Rockhopper.updateProgress();
		    } else {  // Using paired-end reads
			String file1 = parse_readFiles[0];
			String file2 = parse_readFiles[1];
			output("Aligning reads to assembled transcripts using files:\n");
			output("\t" + removePath(file1) + "\n");
			output("\t" + removePath(file2) + "\n");
			Assembler.isPairedEnd = true;
			output("\n\tTotal reads in files:        \t" + DeNovoTranscripts.numReads[readsFileIndex] + "\n");
			output("\tPerfectly aligned reads:     \t" + DeNovoTranscripts.numMappingReads[readsFileIndex] + "\t" + Math.round((100.0*DeNovoTranscripts.numMappingReads[readsFileIndex])/DeNovoTranscripts.numReads[readsFileIndex]) + "%" + "\n\n");
			Rockhopper.updateProgress();
		    }
		    readsFileIndex++;
		}
	    }
	    if (outputSAM) output("Rockhopper did not output a SAM file since sequencing read alignments were previously cached. If you want a SAM file to be output, you must first clear the cache." + "\n\n");
	} else {

	    // Create k-mer dictionary and assemble transcripts de novo
	    kDict = new Dictionary(k);
	    totalReads = new AtomicInteger(0);
	    invalidQualityReads = new AtomicInteger(0);
	    ambiguousReads = new AtomicInteger(0);
	    lowQualityReads = new AtomicInteger(0);
	    mapFullReadsToTranscriptome = false;
	    readsFileIndex = 0;
	    output("\n");
	    for (int i=0; i<conditionFiles.size(); i++) {
		String[] files = conditionFiles.get(i).split(",");
		for (int j=0; j<files.length; j++) {
		    String[] parse_readFiles = files[j].split("%");
		    if (parse_readFiles.length == 1) {  // Using single-end reads (except maybe SAM/BAM)
			String file1 = parse_readFiles[0];
			output("Assembling transcripts from reads in file:             \t" + removePath(file1) + "\n");
			Assembler.isPairedEnd = false;
			Assembler.pairedEndFile = null;
			readsInputFormat = getReadsInputFileType(file1);
			this.PHRED_offset = getPhredOffset(file1);
			if (outputSAM) outputSamHeader(file1);
			int numThreads_OLD = numThreads;
			if ((readsInputFormat == 5) || (readsInputFormat == 6)) numThreads = 1;
			firstLineDone = false;
			readInSequencingReads(file1);  // Read in the sequencing reads
			numThreads = numThreads_OLD;
			if (samWriter != null) {
			    samWriter.close();
			    output("Sequencing reads written to SAM file:                  \t" + output_DIR + getFileNameBase(file1) + ".sam" + "\n");
			}
			Rockhopper.updateProgress();
			Rockhopper.updateProgress();
		    } else {  // Using paired-end reads
			String file1 = parse_readFiles[0];
			String file2 = parse_readFiles[1];
			output("Assembling transcripts from reads in files:\n");
			output("\t" + removePath(file1) + "\n");
			output("\t" + removePath(file2) + "\n");
			Assembler.isPairedEnd = true;
			Assembler.pairedEndFile = file2;
			readsInputFormat = getReadsInputFileType(file1);
			this.PHRED_offset = getPhredOffset(file1);
			if (outputSAM) outputSamHeader(file1);
			int numThreads_OLD = numThreads;
			if ((readsInputFormat == 5) || (readsInputFormat == 6)) numThreads = 1;
			firstLineDone = false;
			readInSequencingReads(file1);  // Read in the sequencing reads
			numThreads = numThreads_OLD;
			if (samWriter != null) {
			    samWriter.close();
			    output("Sequencing reads written to SAM file:                  \t" + output_DIR + getFileNameBase(file1) + ".sam" + "\n");
			}
			Rockhopper.updateProgress();
			Rockhopper.updateProgress();
		    }
		    readsFileIndex++;
		}
	    }
	    kDict.assembleTranscripts();  // Assemble transcripts for final time

	    // Map full length reads to transcripts
	    totalReads = new AtomicInteger(0);
	    invalidQualityReads = new AtomicInteger(0);
	    ambiguousReads = new AtomicInteger(0);
	    lowQualityReads = new AtomicInteger(0);
	    mapFullReadsToTranscriptome = true;
	    kDict.initializeReadMapping(stopAfterOneHit);
	    readsFileIndex = 0;
	    output("\n");
	    for (int i=0; i<conditionFiles.size(); i++) {
		String[] files = conditionFiles.get(i).split(",");
		for (int j=0; j<files.length; j++) {
		    String[] parse_readFiles = files[j].split("%");
		    if (parse_readFiles.length == 1) {  // Using single-end reads (except maybe SAM/BAM)
			String file1 = parse_readFiles[0];
			output("Aligning reads to assembled transcripts using file:\t" + removePath(file1) + "\n");
			Assembler.isPairedEnd = false;
			Assembler.pairedEndFile = null;
			readsInputFormat = getReadsInputFileType(file1);
			this.PHRED_offset = getPhredOffset(file1);
			int numThreads_OLD = numThreads;
			if ((readsInputFormat == 5) || (readsInputFormat == 6)) numThreads = 1;
			firstLineDone = false;
			readInSequencingReads(file1);  // Read in the sequencing reads
			numThreads = numThreads_OLD;
			if (unstranded) kDict.bwtIndex.halveReads(readsFileIndex);  // Fix read double counting
			output("\n\tTotal reads in file:        \t" + kDict.bwtIndex.getNumReads(readsFileIndex) + "\n");
			output("\tPerfectly aligned reads:    \t" + kDict.bwtIndex.getNumMappingReads(readsFileIndex) + "\t" + Math.round((100.0*kDict.bwtIndex.getNumMappingReads(readsFileIndex))/kDict.bwtIndex.getNumReads(readsFileIndex)) + "%" + "\n\n");
			Rockhopper.updateProgress();
		    } else {  // Using paired-end reads
			String file1 = parse_readFiles[0];
			String file2 = parse_readFiles[1];
			output("Aligning reads to assembled transcripts using files:\n");
			output("\t" + removePath(file1) + "\n");
			output("\t" + removePath(file2) + "\n");
			Assembler.isPairedEnd = true;
			Assembler.pairedEndFile = file2;
			readsInputFormat = getReadsInputFileType(file1);
			this.PHRED_offset = getPhredOffset(file1);
			int numThreads_OLD = numThreads;
			if ((readsInputFormat == 5) || (readsInputFormat == 6)) numThreads = 1;
			firstLineDone = false;
			readInSequencingReads(file1);  // Read in the sequencing reads
			numThreads = numThreads_OLD;
			output("\n\tTotal reads in files:        \t" + kDict.bwtIndex.getNumReads(readsFileIndex) + "\n");
			output("\tPerfectly aligned reads:     \t" + kDict.bwtIndex.getNumMappingReads(readsFileIndex) + "\t" + Math.round((100.0*kDict.bwtIndex.getNumMappingReads(readsFileIndex))/kDict.bwtIndex.getNumReads(readsFileIndex)) + "%" + "\n\n");
			Rockhopper.updateProgress();
		    }
		    readsFileIndex++;
		}
	    }
	    transcripts = new DeNovoTranscripts(kDict.bwtIndex);
	    if (computeExpression) transcripts.setMinDiffExpressionLevels(kDict.bwtIndex);
	    kDict.bwtIndex = null;  // Clear Index and Dictionary
	    kDict = null;           // Clear Index and Dictionary
	    System.gc();            // Clear Index and Dictionary
	    if (computeExpression) transcripts.computeDifferentialExpression();
	    FileOps.writeCompressedFile(getCompressedFileName(), conditionFiles, DeNovoTranscripts.numReads, DeNovoTranscripts.numMappingReads, transcripts);
	}

	// Incorporate results from SAM file
	if (sam != null) {
	    totalReads.addAndGet(sam.badReads.get());
	    invalidQualityReads.addAndGet(sam.badReads.get());
	}

	outputToFile(output_DIR + expressionFile, transcripts.toString(labels));
	totalTime = System.currentTimeMillis() - totalTime;
	outputResults();
	summaryWriter.close();
	Rockhopper.updateProgress();
    }



    /*************************************************
     **********   Public Instance Methods   **********
     *************************************************/

    /**
     * Read in sequencing reads from file.
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

	// Using paired-end reads
	String line2 = "";;
	String name2 = "";
	String nextName2 = "";  // Special case (FASTA)
	String read2 = "";
	String quality2 = "";
	String quality_fasta2 = "";  // Special case (FASTA)
	int[] qualities2 = new int[0];
	String[] parse_line2 = null;

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
			    int flag = 4;
			    synchronized(samWriter) {
				samWriter.println(name + "\t" + flag + "\t" + "*" + "\t" + "0" + "\t" + "255" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + parse_line[9] + "\t" + mapPhred(parse_line[10]));
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
			    int flag = 77;    // 01001101
			    int flag2 = 141;  // 10001101
			    synchronized(samWriter) {
				samWriter.println(name + "\t" + flag + "\t" + "*" + "\t" + "0" + "\t" + "255" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + parse_line[9] + "\t" + mapPhred(parse_line[10]));
				samWriter.println(name2 + "\t" + flag2 + "\t" + "*" + "\t" + "0" + "\t" + "255" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + parse_line2[9] + "\t" + mapPhred(parse_line2[10]));
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
			int flag = 4;
			synchronized(samWriter) {
			    samWriter.println(name + "\t" + flag + "\t" + "*" + "\t" + "0" + "\t" + "255" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + parse_line[9] + "\t" + parse_line[10]);
			}
			continue;
		    } else continue;
		} else if (readsInputFormat == 4) {  // SAM file (paired-end)
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
			int flag = 77;    // 01001101
			int flag2 = 141;  // 10001101
			synchronized(samWriter) {
			    samWriter.println(name + "\t" + flag + "\t" + "*" + "\t" + "0" + "\t" + "255" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + parse_line[9] + "\t" + parse_line[10]);
			    samWriter.println(name2 + "\t" + flag2 + "\t" + "*" + "\t" + "0" + "\t" + "255" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + parse_line2[9] + "\t" + parse_line2[10]);
			}
			continue;
		    } else continue;
		}

		// Output to SAM file
		if ((samWriter != null) && (!mapFullReadsToTranscriptome)) {
		    while (name.startsWith("@")) name = name.substring(1);
		    if (isPairedEnd || ((sam != null) && (sam.isPairedEnd()))) {  // Paired-end
			while (name2.startsWith("@")) name2 = name2.substring(1);
			if (name.length() == 0) {
			    name = "r" + nameCount.getAndIncrement();
			    name2 = name;
			}
			int flag = 77;    // 01001101
			int flag2 = 141;  // 10001101
			synchronized(samWriter) {
			    samWriter.println(name + "\t" + flag + "\t" + "*" + "\t" + "0" + "\t" + "255" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + read + "\t" + quality);
			    samWriter.println(name2 + "\t" + flag2 + "\t" + "*" + "\t" + "0" + "\t" + "255" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + read2 + "\t" + quality2);
			}
		    } else {  // Single-end
			if (name.length() == 0) name = "r" + nameCount.getAndIncrement();
			int flag = 4;
			synchronized(samWriter) {
			    samWriter.println(name + "\t" + flag + "\t" + "*" + "\t" + "0" + "\t" + "255" + "\t" + "*" + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + read + "\t" + quality);
			}
		    }
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

		// Process sequencing reads.
		totalReads.getAndIncrement();
		if (ambiguous) ambiguousReads.getAndIncrement();
		else if ((quality.length() < minReadLength) || (isPairedEnd && (quality2.length() < minReadLength)) || ((sam != null) && (sam.isPairedEnd()) && (quality2.length() < minReadLength))) lowQualityReads.getAndIncrement();  // Minimum read length
		else if (invalid) invalidQualityReads.getAndIncrement();  // Invalid quality score
		else {
		    if (!mapFullReadsToTranscriptome) {  // First pass. Assemble transcriptome.
			kDict.add(read);  // Add read to dictionary!
			if (isPairedEnd || ((sam != null) && (sam.isPairedEnd()))) kDict.add(read2);  // Using paired-end reads
		    } else {  // Second pass. Map full length reads to assembled transcripts.
			if (isPairedEnd || ((sam != null) && (sam.isPairedEnd()))) {  // Using paired-end reads
			    kDict.mapFullLengthRead(read, read2);
			} else {  // Using single-end reads
			    kDict.mapFullLengthRead(read);
			    if (unstranded) kDict.mapFullLengthRead(reverseComplement(read));
			}
		    }
		}
	    }
	} catch (IOException e) {
	    output("Error - could not read in from sequencing reads file.\n\n");
	    System.exit(0);
	}
    }

    /**
     * Return the base file name for the compressed transcripts file.
     * For example, if A,B,C are replicates in conditions 1 and 
     * X,Y,Z are replicates in condition 2, then the base file name
     * would be ABC_XYZ__parameters
     */
    public String getCompressedFileName() {
	StringBuilder sb = new StringBuilder();
	for (int i=0; i<conditionFiles.size(); i++) {
	    String[] files = conditionFiles.get(i).split(",");
	    for (int j=0; j<files.length; j++) sb.append(encode(files[j]));
	    sb.append("_");
	}
	sb.append(parametersToString());
	return sb.toString();
    }



    /**********************************************
     **********   Public Class Methods   **********
     **********************************************/

    public static void output(String s) {
	if (summaryWriter != null) summaryWriter.print(s);
	if (output == null) System.out.print(s);
	else {
	    output.append(s);
	    output.setCaretPosition(output.getDocument().getLength());
	}
    }

    /**
     * Returns a reversed version of String s.
     */
    public static String reverse(String s) {
	return (new StringBuilder(s)).reverse().toString();
    }

    /**
     * Returns a reverse complemented version of String s.
     */
    public static String reverseComplement(String s) {
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



    /**************************************************
     **********   Private Instance Methods   **********
     **************************************************/

    /**
     * Read in sequencing reads from file.
     */
    private void readInSequencingReads(String readsFile) {

        try {
	    if (isCompressedFile) reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readsFile))));
	    else reader = new BufferedReader(new FileReader(readsFile));
	    if (isPairedEnd) {
		if (FileOps.isGZIP(pairedEndFile)) reader2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(pairedEndFile))));
		else reader2 = new BufferedReader(new FileReader(pairedEndFile));
	    }
	    Thread[] threads = new Thread[numThreads];
	    for (int i=0; i<numThreads; i++) {
		threads[i] = new ProcessReadsDeNovo_Thread(this);
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
     * Outputs to STDOUT stats on sequencing reads.
     */
    private void outputResults() {
	output("Total number of assembled transcripts:\t\t" + transcripts.getNumTranscripts() + "\n");
	output("\tAverage transcript length:             \t" + transcripts.getAverageTranscriptLength() + "\n");
	output("\tMedian transcript length:              \t" + transcripts.getMedianTranscriptLength() + "\n");
	output("\tTotal number of assembled bases:       \t" + transcripts.getTotalAssembledBases() + "\n");
	output("\n");
	output("Summary of results written to file:               \t" + output_DIR + summaryFile + "\n");
	output("Details of assembled transcripts written to file: \t" + output_DIR + expressionFile + "\n");
	output("\n");
	output("FINSIHED.\n\n");
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
	    sam = new SamOps(readsFile, MAX_LINES, true);
	    if (sam.isValidSam() || sam.isValidBam()) {
		String[] replicons = sam.getReplicons();
		if (sam.isValidSam() && !sam.isPairedEnd()) return 3;  // Single-end SAM file
		else if (sam.isValidSam()) return 4;  // Paired-end SAM file
		else if (sam.isValidBam() && !sam.isPairedEnd()) return 5;  // Single-end BAM file
		else if (sam.isValidBam()) return 6;  // Paired-end BAM file
		else return -1;  // This case should never be reached
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
	    output("\nError - could not open sequencing reads file: " + readsFile + "\n\n");
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
     * Open SAM file for writing.
     * Output header information to SAM file.
     */
    private void outputSamHeader(String readsFile) {
	String samFileName = output_DIR + getFileNameBase(readsFile) + ".sam";
	if (readsFile.equals(samFileName)) return;  // Do not read and write same file at same time
	nameCount = new AtomicInteger(1);
	try {
	    samWriter = new PrintWriter(new File(samFileName));
	    samWriter.println("@HD" + "\t" + "VN:1.0" + "\t" + "SO:unsorted");
	    samWriter.println("@PG" + "\t" + "ID:Rockhopper" + "\t" + "PN:Rockhopper" + "\t" + "VN:" + Rockhopper.version);
	} catch (IOException f) {
	    samWriter = null;
	}
    }

    /**
     * Return a String representation of the parameter settings.
     */
    private String parametersToString() {
	return "_k" + k + "_e" + ("" + computeExpression).toUpperCase().charAt(0) + "_a" + ("" + stopAfterOneHit).toUpperCase().charAt(0) + "_d" + maxPairedEndLength + "_j" + minReadLength + "_n" + CAPACITY_POWER + "_b" + MIN_READS_MAPPING + "_u" + minTranscriptLength + "_w" + minSeedExpression + "_x" + minExpression + "_" + pairedEndOrientation + "_s" + ("" + unstranded).toUpperCase().charAt(0) + "_c" + ("" + singleEndOrientationReverseComplement).toUpperCase().charAt(0);
    }



    /***********************************************
     **********   Private Class Methods   **********
     ***********************************************/

    /**
     * Returns the base of the specified file name, i.e.,
     * the file extension and path to the file are excluded.
     */
    private static String getFileNameBase(String fileName) {
	int slashIndex = fileName.lastIndexOf(separator());
	if (slashIndex == -1) slashIndex = fileName.lastIndexOf("/");
	int periodIndex = fileName.indexOf('.', slashIndex+1);
	if (periodIndex == -1) periodIndex = fileName.length();
	return fileName.substring(slashIndex+1, periodIndex);
    }

    /**
     * Returns the system-dependent file separator that
     * can be used as a RegEx.
     */
    private static String separator() {
	if (File.separatorChar == '\\') return "\\\\";
	else return File.separator;
    }

    /**
     * Returns the file name with the leading path removed.
     */
    private static String removePath(String fileName) {
	int slashIndex = fileName.lastIndexOf(separator());
	if (slashIndex == -1) slashIndex = fileName.lastIndexOf("/");
	return fileName.substring(slashIndex+1);
    }

    /**
     * Returns a 2 character encoding of the String s.
     * First, s is converted to an integer via hashing.
     * Then, the integer is converted into a 2 character
     * encoding, where each character takes on one of
     * 52 different values (a-z, A-Z).
     */
    private static String encode(String s) {
	int encoding = s.hashCode() % 2704;  // 52*52=2704
	if (encoding < 0) encoding += 2704;
	int digit1 = encoding / 52;
	int digit2 = encoding % 52;
	if (digit2 < 0) digit2 += 52;
	char char1 = '?';
	char char2 = '?';
	if (digit1 < 26) char1 = (char)('A' + digit1);
	else char1 = (char)('a' - 26 + digit1);
	if (digit2 < 26) char2 = (char)('A' + digit2);
	else char2 = (char)('a' - 26 + digit2);
	return char1 + "" + char2;
    }

    /**
     * Output the specified string s to the specified file.
     */
    private static void outputToFile(String fileName, String s) {
	try {
	    PrintWriter writer = new PrintWriter(new File(fileName));
	    writer.print(s);
	    writer.close();
	} catch (FileNotFoundException e) {
	    output("\nError - could not open file " + fileName + "\n\n");
	}
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
	if (args.length < 1) {
            output("\nThe Assembler application takes one or more files of RNA-seq reads and it assembles the reads (de novo) into a transcriptome map.\n");
	    output("\nThe Assembler application has one required command line argument.\n");
            output("\nREQUIRED ARGUMENT\n\n");
            output("\t" + "a FASTQ file of RNA-seq reads\n");

            output("\nOPTIONAL ARGUMENTS\n\n");
	    output("\t" + "-k <integer>" + "\t" + "size of k-mer, range of values is 15 to 31 (default is 25)\n");
	    output("\t" + "-c <boolean>" + "\t" + "reverse complement single-end reads (default is false)\n");
            output("\t" + "-ff/fr/rf/rr" + "\t" + "orientation of two mate reads for paired-end read, f=forward and r=reverse_complement (default is fr)\n");
	    output("\t" + "-s <boolean>" + "\t" + "RNA-seq experiments are strand specific (true) or strand ambiguous (false) (default is true)\n");
            output("\t" + "-p <integer>" + "\t" + "number of processors (default is self-identify)\n");
            output("\t" + "-e <boolean>" + "\t" + "compute differential expression for transcripts in pairs of experimental conditions (default is true)\n");
            output("\t" + "-a <boolean>" + "\t" + "report 1 alignment (true) or report all optimal alignments (false), (default is true)\n");
            output("\t" + "-d <integer>" + "\t" + "maximum number of bases between mate pairs for paired-end reads (default is 500)\n");
            output("\t" + "-j <integer>" + "\t" + "minimum length required to use a sequencing read after trimming/processing (default is 35)\n");
            output("\t" + "-n <integer>" + "\t" + "size of k-mer hashtable is ~ 2^n (default is 25). HINT: should normally be 25 or, if more memory is available, 26. WARNING: if increased above 25 then more than 1.2M of memory must be allocated\n");
            output("\t" + "-b <integer>" + "\t" + "minimum number of full length reads required to map to a de novo assembled trancript (default is 20)\n");
           output("\t" + "-u <integer>" + "\t" + "minimum length of de novo assembled transcripts (default is 2*k)\n");
           output("\t" + "-w <integer>" + "\t" + "minimum count of k-mer to use it to seed a new de novo assembled transcript (default is 50)\n");
           output("\t" + "-x <integer>" + "\t" + "minimum count of k-mer to use it to extend an existing de novo assembled transcript (default is 5)\n");
	    output("\t" + "-L <comma separated list>" + "\t" + "labels for each condition\n");
            output("\t" + "-v <boolean>" + "\t" + "verbose output including raw/normalized counts aligning to each transcript (default is false)\n");
	    output("\t" + "-TIME       " + "\t" + "output time taken to execute program\n");

	    output("\nEXAMPLE EXECUTION: SINGLE-END READS\n");
	    output("\njava Assembler <options> aerobic_replicate1.fastq,aerobic_replicate2.fastq anaerobic_replicate1.fastq,anaerobic_replicate2.fastq\n");
	    output("\nEXAMPLE EXECUTION: PAIRED-END READS\n");
	    output("\njava Assembler <options> aerobic_replicate1_pairedend1.fastq%aerobic_replicate1_pairedend2.fastq,aerobic_replicate2_pairedend1.fastq%aerobic_replicate2_pairedend2.fastq anaerobic_replicate1_pairedend1.fastq%anaerobic_replicate1_pairedend2.fastq,anaerobic_replicate2_pairedend1.fastq%anaerobic_replicate2_pairedend2.fastq\n");
	    output("\n");
            System.exit(0);
	}

	// Initially set the number of threads
	Assembler.numThreads = Runtime.getRuntime().availableProcessors();
	if (Assembler.numThreads > 4) Assembler.numThreads = (int)Math.min(Assembler.numThreads * 0.75, 8.0);

	int i=0;
	while (i < args.length) {
	    if (args[i].startsWith("-k")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -k must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.k = Integer.parseInt(args[i+1]);
		if (Assembler.minTranscriptLength == 0) Assembler.minTranscriptLength = 2*Assembler.k;
                i += 2;
	    } else if (args[i].startsWith("-c")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -c must be followed by a boolean.\n");
                    System.exit(0);
                }
		Assembler.singleEndOrientationReverseComplement = Boolean.valueOf(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-s")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -s must be followed by a boolean.\n");
                    System.exit(0);
                }
		Assembler.unstranded = !Boolean.valueOf(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-p")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -p must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.numThreads = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-e")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -e must be followed by a boolean.\n");
                    System.exit(0);
                }
		Assembler.computeExpression = Boolean.valueOf(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-a")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -a must be followed by a boolean.\n");
                    System.exit(0);
                }
		Assembler.stopAfterOneHit = Boolean.valueOf(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-d")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -d must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.maxPairedEndLength = Integer.parseInt(args[i+1]);
		SamOps.maxPairedEndLength = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-j")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -j must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.minReadLength = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-n")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -n must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.CAPACITY_POWER = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-b")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -b must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.MIN_READS_MAPPING = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-u")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -u must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.minTranscriptLength = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-w")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -w must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.minSeedExpression = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-x")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -x must be followed by an integer.\n");
                    System.exit(0);
                }
		Assembler.minExpression = Integer.parseInt(args[i+1]);
                i += 2;
	    } else if (args[i].startsWith("-L")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -L must be followed by a comma separated list of names for the conditions.\n");
                    System.exit(0);
                }
		Assembler.labels = args[i+1].split(",");
		i += 2;
	    } else if (args[i].startsWith("-v")) {
                if ((i == args.length-2) || (args[i+1].startsWith("-"))) {
                    output("Error - command line argument -v must be followed by a boolean.\n");
                    System.exit(0);
                }
		Assembler.verbose = Boolean.valueOf(args[i+1]);
                i += 2;
	    } else if (args[i].equals("-ff")) {
		Assembler.pairedEndOrientation = "ff";
		i++;
	    } else if (args[i].equals("-fr")) {
		Assembler.pairedEndOrientation = "fr";
		i++;
	    } else if (args[i].equals("-rf")) {
		Assembler.pairedEndOrientation = "rf";
		i++;
	    } else if (args[i].equals("-rr")) {
		Assembler.pairedEndOrientation = "rr";
		i++;
	    } else if (args[i].equals("-TIME")) {
		time = true;
		i++;
	    } else {  // Sequencing reads files
		Assembler.conditionFiles.add(args[i]);
		i++;
	    }
	}

	// Handle erroneous command line arguments
	if (conditionFiles.size() == 0) {
	    output("Error - sequencing reads files (in FASTQ, QSEQ, or FASTA format) are required as command line arguments.\n");
	    System.exit(0);
	}
	for (int j=0; j<conditionFiles.size(); j++) {
	    String[] files = conditionFiles.get(j).split(",");
	    for (int k=0; k<files.length; k++) {
		String[] pairedEnd_files = files[k].split("%");
		if (!(new File(pairedEnd_files[0]).exists())) {
		    output("Error - file " + pairedEnd_files[0] + " does not exist.\n");
		    System.exit(0);
		}
		if ((pairedEnd_files.length > 1) && (!(new File(pairedEnd_files[1]).exists()))) {
		    output("Error - file " + pairedEnd_files[1] + " does not exist.\n");
		    System.exit(0);
		}
		if (pairedEnd_files.length > 1) Assembler.unstranded = false;  // Assume paired-end reads are strand specific
	    }
	}

	// If not set by command line arguments, set de novo assembled transcript length
	if (Assembler.minTranscriptLength == 0) Assembler.minTranscriptLength = 2*Assembler.k;

	// Output directory does not exist. Create it.
	if (!(new File(output_DIR).exists())) {
	    if (!(new File(output_DIR).mkdir())) {
		output("Error - could not create directory " + output_DIR + "." + "\n");
		System.exit(0);
	    }
	}

	try {
	    summaryWriter = new PrintWriter(new File(output_DIR + summaryFile));
	} catch (FileNotFoundException e) {
	    output("Error - could not create file " + (output_DIR+summaryFile) + "\n");
	    System.exit(0);
	}
    }



    /*************************************
     **********   Main Method   **********
     *************************************/

    public static void main(String[] args) {

	commandLineArguments(args);
	Assembler a = new Assembler();

    }

}



/***************************************************
 **********   ProcessReadsDeNovo_Thread   **********
 ***************************************************/

class ProcessReadsDeNovo_Thread extends Thread {

    private Assembler a;

    public ProcessReadsDeNovo_Thread(Assembler a) {
	this.a = a;
    }

    public void run() {
	a.processReads();
    }
}
