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

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.zip.GZIPOutputStream;
import java.util.zip.GZIPInputStream;

/**
 * Supports reading and writing of compressed files.
 */
public class FileOps {

    /********************************************
     **********   Instance Variables   **********
     ********************************************/

    private String readFileName;
    private long readFileLength;
    private long readFileDate;
    private String repliconFileName;
    private long repliconFileLength;
    private long repliconFileDate;
    private String repliconHeader;
    private int repliconLength;
    private long avgLengthReads;
    private int[] coordinates_plus;
    private int[] coordinates_minus;
    private ArrayList<DeNovoTranscript> transcripts;  // For de novo assembly
    private boolean isValid;

    // Results
    private int totalReads;
    private int invalidQualityReads;      // Reads containing an invalid quality score
    private int ambiguousReads;           // Reads with ambiguous nucleotides
    private int lowQualityReads;          // Reads with quality scores too low for mapping
    private int mappedOnceReads;          // Reads mapping to one location
    private int mappedMoreThanOnceReads;  // Reads mapping to multiple locations
    private int exactMappedReads;         // Reads mapping without mismatches/errors
    private int inexactMappedReads;       // Reads mapping with mismatches/errors



    /*****************************************
     **********   Class Variables   **********
     *****************************************/

    private static String DIR = Rockhopper.output_DIR + "intermediary/";
    private static byte[] buffer1 = new byte[1];  // Temp array used when reading in GZIP file



    /**************************************
     **********   Constructors   **********
     **************************************/

    public FileOps(String inFileName) {
	isValid = readCompressedFile(inFileName);
    }

    public FileOps(String inFileName, String readFileName, String repliconFileName) {
	isValid = readCompressedFile(inFileName);
	File readFile = new File(readFileName);
	File repliconFile = new File(repliconFileName);
	isValid = isValid && (readFile.getName().equals(this.readFileName)) && (readFile.length() == readFileLength) && (readFile.lastModified() == readFileDate) && (repliconFile.getName().equals(this.repliconFileName)) && (repliconFile.length() == repliconFileLength) && (repliconFile.lastModified() == repliconFileDate);
    }

    /**
     * Constructor used for de novo transcript assembly.
     */
    public FileOps(String inFileName, ArrayList<String> conditionFiles) {
	isValid = readCompressedFile_DeNovo(inFileName, conditionFiles);
    }



    /*************************************************
     **********   Public Instance Methods   **********
     *************************************************/

    /**
     * Returns true if the input file is valid, i.e.,
     * the stats of the read-file and replicon-file file 
     * correspond to the stats from the input file.
     * Returns false otherwise.
     */
    public boolean isValid() {
	return this.isValid;
    }

    /**
     * Returns the name of the sequencing read file.
     */
    public String getReadFileName() {
	return this.readFileName;
    }

    /**
     * Returns the average length of sequencing reads
     * in the sequencing reads file.
     */
    public long getAvgLengthReads() {
	return avgLengthReads;
    }

    /**
     * Returns plus strand coordinates.
     */
    public int[] getCoordinates_plus() {
	return this.coordinates_plus;
    }

    /**
     * Returns minus strand coordinates.
     */
    public int[] getCoordinates_minus() {
	return this.coordinates_minus;
    }

    /**
     * Returns the total reads
     */
    public int getTotalReads() {
	return this.totalReads;
    }

    /**
     * Returns the invalide quality reads
     */
    public int getInvalidQualityReads() {
	return this.invalidQualityReads;
    }

    /**
     * Returns the ambiguous reads
     */
    public int getAmbiguousReads() {
	return this.ambiguousReads;
    }

    /**
     * Returns the low quality reads
     */
    public int getLowQualityReads() {
	return this.lowQualityReads;
    }

    /**
     * Returns the mapped once reads
     */
    public int getMappedOnceReads() {
	return this.mappedOnceReads;
    }

    /**
     * Returns the mapped more than once reads
     */
    public int getMappedMoreThanOnceReads() {
	return this.mappedMoreThanOnceReads;
    }

    /**
     * Returns the exact mapped reads
     */
    public int getExactMappedReads() {
	return this.exactMappedReads;
    }

    /**
     * Returns the inexact mapped reads
     */
    public int getInexactMappedReads() {
	return this.inexactMappedReads;
    }

    /**
     * Returns the de novo assembled transcripts
     */
    public ArrayList<DeNovoTranscript> getTranscripts() {
	return this.transcripts;
    }



    /**********************************************
     **********   Public Class Methods   **********
     **********************************************/

    /**
     * Converts the specified pair of int arrays to byte arrays
     * and writes these byte arrays to a compressed file. The following
     * header information is also written to the compressed output file:
     *      version
     *      read file name
     *      read file length
     *      read file last modified date
     *      replicon file name
     *      replicon file length
     *      replicon file last modified date
     *      replicon FASTA header line
     *      replicon length
     *      average length of reads
     *      total reads
     *      invalid quality reads
     *      ambiguous reads
     *      low quality reads
     *      mapped once reads
     *      mapped more than once reads
     *      exact mapped reads
     *      inexact mapped reads
     */
    public static void writeCompressedFile(String outFileName, String readFileName, String repliconFileName, String repliconHeader, int repliconLength, int totalReads, int invalidQualityReads, int ambiguousReads, int lowQualityReads, int mappedOnceReads, int mappedMoreThanOnceReads, int exactMappedReads, int inexactMappedReads, long avgLengthReads, int[] a1, int[] a2) {
	File readFile = new File(readFileName);
	File repliconFile = new File(repliconFileName);
	byte[] header = (Rockhopper.version + "\n" + readFile.getName() + "\n" + readFile.length() + "\n" + readFile.lastModified() + "\n" + repliconFile.getName() + "\n" + repliconFile.length() + "\n" + repliconFile.lastModified() + "\n" + repliconHeader + "\n" + repliconLength + "\n" + avgLengthReads + "\n" + totalReads + "\n" + invalidQualityReads + "\n" + ambiguousReads + "\n" + lowQualityReads + "\n" + mappedOnceReads + "\n" + mappedMoreThanOnceReads + "\n" + exactMappedReads + "\n" + inexactMappedReads + "\n").getBytes();
	byte[] b1 = intArrayToByteArray(a1);
	byte[] b2 = intArrayToByteArray(a2);
	try {
	    // Set up directory
	    File d = new File(DIR);
	    if (!d.isDirectory()) d.mkdir();
	    GZIPOutputStream gzos = new GZIPOutputStream(new FileOutputStream(DIR + outFileName + ".gz"));
	    gzos.write(header, 0, header.length);
	    gzos.write(b1, 0, b1.length);
	    gzos.write(b2, 0, b2.length);
	    gzos.close();
	} catch (IOException e) {
	    Rockhopper.output("Could not write to file " + outFileName + ".gz" + "\n");
	}
    }

    /**
     * Writes a list of de novo assembled transcripts to a compressed file.
     * The following header information is also written to the compressed file:
     *      version
     *      number of conditions
     *      number of replicates in each condition
     *      for each read file: its name, its length, and its last modified time
     *      number of files
     *      for each file: total reads, total mapping reads
     *      number of transcripts in list
     *      number of p-values per transcript
     *      for each transcript:
     *            sequence
     *            for each condition: rawCounts_nts, rawCounts_reads, normalizedCounts
     *            RPKMs, means, variances, lowess, pValues, qValues
     */
    public static void writeCompressedFile(String outFileName, ArrayList<String> conditionFiles, int[] numReads, int[] numMappingReads, DeNovoTranscripts transcripts) {

	// Version, number of conditions, and number of replicates in each condition
	StringBuilder header = new StringBuilder();
	header.append(Rockhopper.version + "\n");
	header.append(conditionFiles.size() + "\n");
	for (int i=0; i<conditionFiles.size(); i++) {
	    header.append(conditionFiles.get(i).split(",").length + "\n");
	}

	// Name, length, and last modified time of each file (including mate-pair files)
	for (int i=0; i<conditionFiles.size(); i++) {
	    String[] files = conditionFiles.get(i).split(",");
	    for (int j=0; j<files.length; j++) {
		String[] parse_readFiles = files[j].split("%");
		if (parse_readFiles.length == 1) {  // Using single-end reads
		    String file1 = parse_readFiles[0];
		    File readFile = new File(file1);
		    header.append(readFile.getName() + "\n" + readFile.length() + "\n" + readFile.lastModified() + "\n");
		} else {  // Using paired-end reads
		    String file1 = parse_readFiles[0];
		    File readFile1 = new File(file1);
		    header.append(readFile1.getName() + "\n" + readFile1.length() + "\n" + readFile1.lastModified() + "\n");
		    String file2 = parse_readFiles[1];
		    File readFile2 = new File(file2);
		    header.append(readFile2.getName() + "\n" + readFile2.length() + "\n" + readFile2.lastModified() + "\n");
		}
	    }
	}

	// Number of files. Total reads per file and total number of mapping reads per file
	header.append(numReads.length + "\n");
	for (int i=0; i<numReads.length; i++) header.append(numReads[i] + "\n" + numMappingReads[i] + "\n");

	// Number of transcripts and number of p-values per transcript
	header.append(transcripts.getNumTranscripts() + "\n");
	if (transcripts.getNumTranscripts() == 0) header.append("0" + "\n");
	else header.append(transcripts.getTranscript(0).getNumPvalues() + "\n");

	// Write header information and information for each transcript
	try {
	    // Set up directory
	    File d = new File(DIR);
	    if (!d.isDirectory()) d.mkdir();
	    GZIPOutputStream gzos = new GZIPOutputStream(new FileOutputStream(DIR + outFileName + ".gz"));
	    byte[] header_bytes = header.toString().getBytes();
	    gzos.write(header_bytes, 0, header_bytes.length);
	    for (int i=0; i<transcripts.getNumTranscripts(); i++) {
		DeNovoTranscript transcript = transcripts.getTranscript(i);
		byte[] b = (transcript.sequence() + "\n").getBytes();  // Add new line to sequence
		gzos.write(b, 0, b.length);
		for (int j=0; j<conditionFiles.size(); j++) {
		    b = longArrayToByteArray(transcript.getRawCounts_nts_array(j));
		    gzos.write(b, 0, b.length);
		    b = longArrayToByteArray(transcript.getRawCounts_reads_array(j));
		    gzos.write(b, 0, b.length);
		    b = longArrayToByteArray(transcript.getNormalizedCounts_array(j));
		    gzos.write(b, 0, b.length);
		}
		b = longArrayToByteArray(transcript.getRPKMs());
		gzos.write(b, 0, b.length);
		b = longArrayToByteArray(transcript.getMeans());
		gzos.write(b, 0, b.length);
		b = longArrayToByteArray(transcript.getVariances());
		gzos.write(b, 0, b.length);
		b = longArrayToByteArray(transcript.getLowess());
		gzos.write(b, 0, b.length);
		b = doubleArrayToByteArray(transcript.getPvalues());
		gzos.write(b, 0, b.length);
		b = doubleArrayToByteArray(transcript.getQvalues());
		gzos.write(b, 0, b.length);
	    }
	    gzos.close();
	} catch (IOException e) {
	    Rockhopper.output("Could not write to file " + outFileName + ".gz" + "\n");
	}
    }

    /**
     * Returns true is file is successfully read. False otherwise.
     */
    public boolean readCompressedFile(String inFileName) {
	try {
	    GZIPInputStream in = new GZIPInputStream(new FileInputStream(DIR + inFileName + ".gz"));

	    // Read in header information
	    String version = getLine(in);
	    if (!version.equals(Rockhopper.version)) return false;
	    this.readFileName = getLine(in);
	    this.readFileLength = Long.parseLong(getLine(in));
	    this.readFileDate = Long.parseLong(getLine(in));
	    this.repliconFileName = getLine(in);
	    this.repliconFileLength = Long.parseLong(getLine(in));
	    this.repliconFileDate = Long.parseLong(getLine(in));
	    this.repliconHeader = getLine(in);
	    this.repliconLength = Integer.parseInt(getLine(in));
	    this.avgLengthReads = Long.parseLong(getLine(in));
	    this.totalReads = Integer.parseInt(getLine(in));
	    this.invalidQualityReads = Integer.parseInt(getLine(in));
	    this.ambiguousReads = Integer.parseInt(getLine(in));
	    this.lowQualityReads = Integer.parseInt(getLine(in));
	    this.mappedOnceReads = Integer.parseInt(getLine(in));
	    this.mappedMoreThanOnceReads = Integer.parseInt(getLine(in));
	    this.exactMappedReads = Integer.parseInt(getLine(in));
	    this.inexactMappedReads = Integer.parseInt(getLine(in));

	    // Retrieve data
	    this.coordinates_plus = new int[repliconLength];
	    for (int i=0; i<repliconLength; i++) coordinates_plus[i] = getInt(in);
	    this.coordinates_minus = new int[repliconLength];
	    for (int i=0; i<repliconLength; i++) coordinates_minus[i] = getInt(in);
	    in.close();
	    return true;
	} catch (IOException e) {
	    return false;
	}
    }

    /**
     * Returns true is file is successfully read and is valid. False otherwise.
     */
    public boolean readCompressedFile_DeNovo(String inFileName, ArrayList<String> conditionFiles) {
	try {

	    // Read in file
	    GZIPInputStream in = new GZIPInputStream(new FileInputStream(DIR + inFileName + ".gz"));
	    ArrayList<Byte> bytes = new ArrayList<Byte>(1000000);
	    byte[] buffer = new byte[1024];
	    int len;
	    while ((len = in.read(buffer)) >= 0) {
		for (int i=0; i<len; i++)
		    bytes.add(buffer[i]);
	    }
	    in.close();

	    // First header line: version
	    int startIndex = 0;
	    int stopIndex = getStopIndex(bytes, startIndex);
	    String version = getLine(bytes, startIndex, stopIndex);
	    if (!version.equals(Rockhopper.version)) return false;

	    // Second header line: number of conditions
	    startIndex = stopIndex + 1;
	    stopIndex = getStopIndex(bytes, startIndex);
	    int numConditions = Integer.parseInt(getLine(bytes, startIndex, stopIndex));
	    startIndex = stopIndex + 1;

	    // One line for each condition indicating number of replicates in condition
	    int[] numReplicates = new int[numConditions];
	    for (int i=0; i<numConditions; i++) {
		stopIndex = getStopIndex(bytes, startIndex);
		numReplicates[i] = Integer.parseInt(getLine(bytes, startIndex, stopIndex));
		startIndex = stopIndex + 1;
	    }

	    // For each read file, we have its name, length, and last modified time
	    String[][][] names = new String[2][numConditions][];  // First coordinate is mate-pair
	    long[][][] lengths = new long[2][numConditions][];    // First coordinate is mate-pair
	    long[][][] times = new long[2][numConditions][];      // First coordinate is mate-pair
	    for (int i=0; i<numConditions; i++) {
		names[0][i] = new String[numReplicates[i]];
		lengths[0][i] = new long[numReplicates[i]];
		times[0][i] = new long[numReplicates[i]];
		names[1][i] = new String[numReplicates[i]];
		lengths[1][i] = new long[numReplicates[i]];
		times[1][i] = new long[numReplicates[i]];
		String[] files = conditionFiles.get(i).split(",");
		for (int j=0; j<numReplicates[i]; j++) {
		    stopIndex = getStopIndex(bytes, startIndex);
		    names[0][i][j] = getLine(bytes, startIndex, stopIndex);
		    startIndex = stopIndex + 1;
		    stopIndex = getStopIndex(bytes, startIndex);
		    lengths[0][i][j] = Long.parseLong(getLine(bytes, startIndex, stopIndex));
		    startIndex = stopIndex + 1;
		    stopIndex = getStopIndex(bytes, startIndex);
		    times[0][i][j] = Long.parseLong(getLine(bytes, startIndex, stopIndex));
		    startIndex = stopIndex + 1;

		    // Handle mate-pair files, if any
		    String[] parse_readFiles = files[j].split("%");
		    if (parse_readFiles.length > 1) {  // Using paired-end reads
			stopIndex = getStopIndex(bytes, startIndex);
			names[1][i][j] = getLine(bytes, startIndex, stopIndex);
			startIndex = stopIndex + 1;
			stopIndex = getStopIndex(bytes, startIndex);
			lengths[1][i][j] = Long.parseLong(getLine(bytes, startIndex, stopIndex));
			startIndex = stopIndex + 1;
			stopIndex = getStopIndex(bytes, startIndex);
			times[1][i][j] = Long.parseLong(getLine(bytes, startIndex, stopIndex));
			startIndex = stopIndex + 1;
		    }
		}
	    }

	    // Number of files. For each file, number of reads and number of mapping reads
	    stopIndex = getStopIndex(bytes, startIndex);
	    int numFiles = Integer.parseInt(getLine(bytes, startIndex, stopIndex));
	    startIndex = stopIndex + 1;
	    int[] numReads = new int[numFiles];
	    int[] numMappingReads = new int[numFiles];
	    for (int i=0; i<numFiles; i++) {
		stopIndex = getStopIndex(bytes, startIndex);
		numReads[i] = Integer.parseInt(getLine(bytes, startIndex, stopIndex));
		startIndex = stopIndex + 1;
		stopIndex = getStopIndex(bytes, startIndex);
		numMappingReads[i] = Integer.parseInt(getLine(bytes, startIndex, stopIndex));
		startIndex = stopIndex + 1;
	    }

	    // Final two header lines
	    stopIndex = getStopIndex(bytes, startIndex);
	    int numTranscripts = Integer.parseInt(getLine(bytes, startIndex, stopIndex));
	    startIndex = stopIndex + 1;
	    stopIndex = getStopIndex(bytes, startIndex);
	    int numPvalues = Integer.parseInt(getLine(bytes, startIndex, stopIndex));
	    startIndex = stopIndex + 1;

	    // Transcripts
	    ArrayList<DeNovoTranscript> transcripts = new ArrayList<DeNovoTranscript>(numTranscripts);
	    for (int i=0; i<numTranscripts; i++) {  // For each transcript
		stopIndex = getStopIndex(bytes, startIndex);
		String sequence = getLine(bytes, startIndex, stopIndex);
		startIndex = stopIndex + 1;
		ArrayList<ArrayList<Long>> rawCounts_nts = new ArrayList<ArrayList<Long>>(numConditions);
		ArrayList<ArrayList<Long>> rawCounts_reads = new ArrayList<ArrayList<Long>>(numConditions);
		ArrayList<ArrayList<Long>> normalizedCounts = new ArrayList<ArrayList<Long>>(numConditions);
		for (int j=0; j<numConditions; j++) {
		    byte[] tempData = new byte[numReplicates[j]*8];
		    for (int k=startIndex; k<startIndex+numReplicates[j]*8; k++) tempData[k-startIndex] = bytes.get(k).byteValue();
		    rawCounts_nts.add(byteToLongList(tempData));
		    startIndex += numReplicates[j]*8;
		    for (int k=startIndex; k<startIndex+numReplicates[j]*8; k++) tempData[k-startIndex] = bytes.get(k).byteValue();
		    rawCounts_reads.add(byteToLongList(tempData));
		    startIndex += numReplicates[j]*8;
		    for (int k=startIndex; k<startIndex+numReplicates[j]*8; k++) tempData[k-startIndex] = bytes.get(k).byteValue();
		    normalizedCounts.add(byteToLongList(tempData));
		    startIndex += numReplicates[j]*8;
		}
		byte[] tempData = new byte[numConditions*8];
		for (int k=startIndex; k<startIndex+numConditions*8; k++) tempData[k-startIndex] = bytes.get(k).byteValue();
		ArrayList<Long> RPKMs = byteToLongList(tempData);
		startIndex += numConditions*8;
		for (int k=startIndex; k<startIndex+numConditions*8; k++) tempData[k-startIndex] = bytes.get(k).byteValue();
		ArrayList<Long> means = byteToLongList(tempData);
		startIndex += numConditions*8;
		tempData = new byte[numConditions*numConditions*8];
		for (int k=startIndex; k<startIndex+numConditions*numConditions*8; k++) tempData[k-startIndex] = bytes.get(k).byteValue();
		long[][] variances = byteToLong2Darray(tempData);
		startIndex += numConditions*numConditions*8;
		for (int k=startIndex; k<startIndex+numConditions*numConditions*8; k++) tempData[k-startIndex] = bytes.get(k).byteValue();
		long[][] lowess = byteToLong2Darray(tempData);
		startIndex += numConditions*numConditions*8;
		tempData = new byte[numPvalues*8];
		for (int k=startIndex; k<startIndex+numPvalues*8; k++) tempData[k-startIndex] = bytes.get(k).byteValue();
		ArrayList<Double> pValues = byteToDoubleList(tempData);
		startIndex += numPvalues*8;
		for (int k=startIndex; k<startIndex+numPvalues*8; k++) tempData[k-startIndex] = bytes.get(k).byteValue();
		ArrayList<Double> qValues = byteToDoubleList(tempData);
		startIndex += numPvalues*8;
		transcripts.add(new DeNovoTranscript(sequence, rawCounts_nts, rawCounts_reads, normalizedCounts, RPKMs, means, variances, lowess, pValues, qValues));
	    }

	    // Check validity
	    if (numConditions != conditionFiles.size()) return false;
	    for (int i=0; i<numConditions; i++) {
		String[] files = conditionFiles.get(i).split(",");
		if (numReplicates[i] != files.length) return false;
		for (int j=0; j<numReplicates[i]; j++) {
		    String[] parse_readFiles = files[j].split("%");
		    if (parse_readFiles.length == 1) {  // Using single-end reads
			String file1 = parse_readFiles[0];
			File f = new File(file1);
			if (!names[0][i][j].equals(f.getName())) return false;
			if (lengths[0][i][j] != f.length()) return false;
			if (times[0][i][j] != f.lastModified()) return false;
		    } else {  // Using paired-end reads
			String file1 = parse_readFiles[0];
			String file2 = parse_readFiles[1];
			File f1 = new File(file1);
			if (!names[0][i][j].equals(f1.getName())) return false;
			if (lengths[0][i][j] != f1.length()) return false;
			if (times[0][i][j] != f1.lastModified()) return false;
			File f2 = new File(file2);
			if (!names[1][i][j].equals(f2.getName())) return false;
			if (lengths[1][i][j] != f2.length()) return false;
			if (times[1][i][j] != f2.lastModified()) return false;
		    }
		}
	    }
	    this.transcripts = transcripts;
	    DeNovoTranscripts.numReads = numReads;
	    DeNovoTranscripts.numMappingReads = numMappingReads;
	    return true;
	} catch (Exception e) {
	    return false;
	}
    }

    /**
     * Returns true if the specified file is in gzip format.
     * Returns false otherwise.
     */
    public static boolean isGZIP(String fileName) {
	try {
	    GZIPInputStream in = new GZIPInputStream(new FileInputStream(fileName));
	    in.close();
	    return true;
	} catch (Exception e) {
	    return false;
	}
    }



    /**************************************************
     **********   Private Instance Methods   **********
     **************************************************/

    /**
     * Returns the first index after the specified startIndex in the Byte list 
     * corresponding to a new line (end of line).
     */
    private int getStopIndex(ArrayList<Byte> bytes, int startIndex) {
	int stopIndex = startIndex;
	while ((stopIndex < bytes.size()) && (bytes.get(stopIndex).byteValue() != 10)) stopIndex++;
	return stopIndex;
    }

    /**
     * Returns a String corresponding to one line, i.e., the bytes from
     * startIndex to stopIndex.
     */
    private String getLine(ArrayList<Byte> bytes, int startIndex, int stopIndex) {
	// Convert sequence of bytes to a String representing a line
	byte[] bytesInLine = new byte[stopIndex - startIndex];
	for (int i=startIndex; i<stopIndex; i++) bytesInLine[i - startIndex] = bytes.get(i).byteValue();
	return (new String(bytesInLine));
    }



    /***********************************************
     **********   Private Class Methods   **********
     ***********************************************/

    /**
     * Returns one line from the specified input stream.
     */
    private static String getLine(GZIPInputStream in) {
	try {
	    int len = -1;
	    ArrayList<Byte> buffer = new ArrayList<Byte>();
	    while (((len = in.read(buffer1)) >= 0) && (buffer1[0] != 10)) {
		buffer.add(buffer1[0]);
	    }
	    byte[] buffer_bytes = new byte[buffer.size()];
	    for (int i=0; i<buffer.size(); i++) buffer_bytes[i] = buffer.get(i).byteValue();
	    return (new String(buffer_bytes));
	} catch (IOException e) {
	    return "";
	}
    }

    /**
     * Returns one integer from the specified input stream.
     */
    private static int getInt(GZIPInputStream in) {
	try {
	    int result = 0;
	    in.read(buffer1);
	    result += (buffer1[0] & 0xFF) << 24;
	    in.read(buffer1);
	    result += (buffer1[0] & 0xFF) << 16;
	    in.read(buffer1);
	    result += (buffer1[0] & 0xFF) << 8;
	    in.read(buffer1);
	    result += (buffer1[0] & 0xFF);
	    return result;
	} catch (IOException e) {
	    return -1;
	}
    }

    /**
     * Returns a byte array representation of the specified int array.
     */
    private static byte[] intArrayToByteArray(int[] a) {
	byte[] b = new byte[a.length*4];
	for (int i=0; i<a.length; i++) {
	    b[4*i+0] = (byte)(a[i] >>> 24);
	    b[4*i+1] = (byte)(a[i] >>> 16);
	    b[4*i+2] = (byte)(a[i] >>> 8);
	    b[4*i+3] = (byte)(a[i] >>> 0);
	}
	return b;
    }

    /**
     * Returns a byte array representation of the specified long array.
     */
    private static byte[] longArrayToByteArray(long[] a) {
	byte[] b = new byte[a.length*8];
	for (int i=0; i<a.length; i++) {
	    b[8*i+0] = (byte)(a[i] >>> 56);
	    b[8*i+1] = (byte)(a[i] >>> 48);
	    b[8*i+2] = (byte)(a[i] >>> 40);
	    b[8*i+3] = (byte)(a[i] >>> 32);
	    b[8*i+4] = (byte)(a[i] >>> 24);
	    b[8*i+5] = (byte)(a[i] >>> 16);
	    b[8*i+6] = (byte)(a[i] >>> 8);
	    b[8*i+7] = (byte)(a[i] >>> 0);
	}
	return b;
    }

    /**
     * Returns a byte array representation of the specified double array.
     */
    private static byte[] doubleArrayToByteArray(double[] a) {
	byte[] b = new byte[a.length*8];
	for (int i=0; i<a.length; i++) {
	    b[8*i+0] = (byte)(Double.doubleToRawLongBits(a[i]) >>> 56);
	    b[8*i+1] = (byte)(Double.doubleToRawLongBits(a[i]) >>> 48);
	    b[8*i+2] = (byte)(Double.doubleToRawLongBits(a[i]) >>> 40);
	    b[8*i+3] = (byte)(Double.doubleToRawLongBits(a[i]) >>> 32);
	    b[8*i+4] = (byte)(Double.doubleToRawLongBits(a[i]) >>> 24);
	    b[8*i+5] = (byte)(Double.doubleToRawLongBits(a[i]) >>> 16);
	    b[8*i+6] = (byte)(Double.doubleToRawLongBits(a[i]) >>> 8);
	    b[8*i+7] = (byte)(Double.doubleToRawLongBits(a[i]) >>> 0);
	}
	return b;
    }

    /**
     * Returns an int array representation of the specified byte array.
     */
    private static int[] byteToIntArray(byte[] b) {
	int[] a = new int[b.length/4];
	for (int i=0; i<b.length; i+=4) {
	    a[i/4] = (b[i+0] << 24) + ((b[i+1] & 0xFF) << 16) + ((b[i+2] & 0xFF) << 8) + (b[i+3] & 0xFF);
	}
	return a;
    }

    /**
     * Returns a long list representation of the specified byte array.
     */
    private static ArrayList<Long> byteToLongList(byte[] b) {
	ArrayList<Long> a = new ArrayList<Long>(b.length/8);
	for (int i=0; i<b.length; i+=8) {
	    a.add(((long)(b[i+0] & 0xFF) << 56) | ((long)(b[i+1] & 0xFF) << 48) | ((long)(b[i+2] & 0xFF) << 40) | ((long)(b[i+3] & 0xFF) << 32) | ((long)(b[i+4] & 0xFF) << 24) | ((long)(b[i+5] & 0xFF) << 16) | ((long)(b[i+6] & 0xFF) << 8) | ((long)(b[i+7] & 0xFF) << 0));
	}
	return a;
    }

    /**
     * Returns a double list representation of the specified byte array.
     */
    private static ArrayList<Double> byteToDoubleList(byte[] b) {
	ArrayList<Double> a = new ArrayList<Double>(b.length/8);
	for (int i=0; i<b.length; i+=8) {
	    a.add(Double.longBitsToDouble(((long)(b[i+0] & 0xFF) << 56) | ((long)(b[i+1] & 0xFF) << 48) | ((long)(b[i+2] & 0xFF) << 40) | ((long)(b[i+3] & 0xFF) << 32) | ((long)(b[i+4] & 0xFF) << 24) | ((long)(b[i+5] & 0xFF) << 16) | ((long)(b[i+6] & 0xFF) << 8) | ((long)(b[i+7] & 0xFF) << 0)));
	}
	return a;
    }

    /**
     * Returns a 2D long array representation of the specified byte array.
     */
    private static long[][] byteToLong2Darray(byte[] b) {
	ArrayList<Long> a = byteToLongList(b);

	// Convert 1D long list into 2D long array
	int size = (int)(Math.sqrt(a.size()));
	long[][] c = new long[size][size];
	int index = 0;
	for (int x=0; x<size; x++) {
	    for (int y=0; y<size; y++) {
		c[x][y] = a.get(index);
		index++;
	    }
	}
	return c;
    }



    /*************************************
     **********   Main Method   **********
     *************************************/

    public static void main(String[] args) {

    }

}

