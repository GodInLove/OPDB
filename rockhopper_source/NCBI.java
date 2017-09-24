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
import java.util.HashMap;
import java.util.Scanner;
import java.net.URL;
import java.net.URLConnection;
import java.io.IOException;
import java.io.File;
import java.io.PrintWriter;

/**
 * Maintains a list of bacterial genomes and their annotations available from NCBI.
 */
public class NCBI {

    /*****************************************
     **********   CLASS VARIABLES   **********
     *****************************************/

    public static int HASH_SIZE = 3;
    private static String fileName = "replicons.txt";
    private static String Genome_url = "http://cs.wellesley.edu/~btjaden/genomes/";
    private static String Rockhopper_url = "http://cs.wellesley.edu/~btjaden/Rockhopper/";
    public static String genome_DIR = Rockhopper.output_DIR + "genomes/";



    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private boolean Genome_valid;  // We are able to establish a connection to the server.
    private boolean Rockhopper_valid;  // We are able to establish a connection to the server.
    private ArrayList<String> repliconList;
    private HashMap<String,HashMap<Integer,Boolean>> repliconDictionary;



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    public NCBI() {
	Rockhopper_valid = establishServerConnection(Rockhopper_url + fileName);
	repliconList = readInRepliconList();
	repliconDictionary = createRepliconDictionary(repliconList);
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    /**
     * Returns true if a connection to the Genome server has been successfully established.
     * Returns false otherwise.
     */
    public boolean isValid_Genome() {
	return this.Genome_valid;
    }

    /**
     * Returns true if a connection to the NCBI server has been successfully established.
     * Returns false otherwise.
     */
    public boolean isValid_Rockhopper() {
	return this.Rockhopper_valid;
    }

    /**
     * Return a list of all replicon names that contain
     * the substring s.
     */
    public ArrayList<String> getRepliconNames(String s) {
	s = s.toLowerCase();
	ArrayList<String> reps = new ArrayList<String>();
	if (s.length() < HASH_SIZE) return reps;
	String seed = s.substring(0, HASH_SIZE);
	if (!repliconDictionary.containsKey(seed)) return reps;
	HashMap<Integer,Boolean> list = repliconDictionary.get(seed);
	for (Integer i : list.keySet()) {
	    String repliconName = repliconList.get(i.intValue()).split("\t")[0];
	    if (repliconName.toLowerCase().indexOf(s) >= 0) reps.add(repliconName);
	}
	Collections.sort(reps);
	return reps;
    }

    /**
     * Return the name of a FASTA file on the Genome server corresponding
     * to the given replicon name.
     */
    public String getRepliconFileName(String repliconName) {
	repliconName = repliconName.toLowerCase();
	if (repliconName.length() < HASH_SIZE) return "";
	String seed = repliconName.substring(0, HASH_SIZE);
	if (!repliconDictionary.containsKey(seed)) return "";
	HashMap<Integer,Boolean> list = repliconDictionary.get(seed);
	for (Integer i : list.keySet()) {
	    String[] parse_repliconInfo = repliconList.get(i.intValue()).split("\t");
	    if (repliconName.equalsIgnoreCase(parse_repliconInfo[0])) return parse_repliconInfo[1];
	}
	return "";
    }

    /**
     * Ensures the appropriate genomes files (*.fna, *.ptt, *.rnt) are
     * available locally for the specified replicon. If not, the
     * files are copied from the server and created locally.
     * Returns true if successful, false otherwise (in which
     * case the user must specifify the files locally).
     */
    public boolean prepGenomeFiles(String repliconName) {
	boolean successful = true;
	String serverFileName = getRepliconFileName(repliconName);
	if (serverFileName.length() == 0) return false;
	int slashIndex = serverFileName.lastIndexOf("/");
	int periodIndex = serverFileName.lastIndexOf(".");
	String baseFileName = serverFileName.substring(slashIndex+1, periodIndex);
	String replicon_DIR = getRepliconDir(repliconName);
	File resultsDIR = new File(Rockhopper.output_DIR);
	if (!resultsDIR.exists()) resultsDIR.mkdir();
	File genDIR = new File(genome_DIR);
	if (!genDIR.exists()) genDIR.mkdir();
	File repliconDIR = new File(genome_DIR + replicon_DIR);
	if (!repliconDIR.exists()) repliconDIR.mkdir();
	File genomeFile = new File(genome_DIR + replicon_DIR + baseFileName + ".fna");
	if (!genomeFile.exists()) successful = downloadFile(serverFileName, genome_DIR + replicon_DIR + baseFileName + ".fna");
	File geneFile = new File(genome_DIR + replicon_DIR + baseFileName + ".ptt");
	if (!geneFile.exists()) downloadFile(serverFileName.substring(0, serverFileName.length()-4) + ".ptt", genome_DIR + replicon_DIR + baseFileName + ".ptt");
	File rnaFile = new File(genome_DIR + replicon_DIR + baseFileName + ".rnt");
	if (!rnaFile.exists()) downloadFile(serverFileName.substring(0, serverFileName.length()-4) + ".rnt", genome_DIR + replicon_DIR + baseFileName + ".rnt");
	return successful;
    }



    /**************************************************
     **********   PRIVATE INSTANCE METHODS   **********
     **************************************************/

    /**
     * Returns true if we can connect to the servers, false otherwise.
     */
    private boolean establishServerConnection(String url) {
	try {
	    Scanner reader = new Scanner(new URL(url).openStream());
	    return true;
	} catch (IOException e) {
	    return false;
	}
    }

    /**
     * Read in the list of replicons either from a local file or from
     * the Rockhopper server. If there is no local file or its file
     * size differs from that of the server file, then a new local
     * file is created. Returns an empty list if neither a local
     * file nor the server file could be found/accessed.
     */
    private ArrayList<String> readInRepliconList() {
	ArrayList<String> repliconList = new ArrayList<String>();

	try {
	    // Local file
	    File resultsDIR = new File(Rockhopper.output_DIR);
	    if (!resultsDIR.exists()) resultsDIR.mkdir();
	    File genDIR = new File(genome_DIR);
	    if (!genDIR.exists()) genDIR.mkdir();
	    File repFile = new File(genome_DIR + fileName);
	    long fileLength = 0;
	    if (repFile.exists()) fileLength = repFile.length();

	    // Server file
	    long serverLength = 0;
	    if (isValid_Rockhopper())
		serverLength = new URL(Rockhopper_url + fileName).openConnection().getContentLength();

	    // Neither local file nor server file is available.
	    if (!repFile.exists() && !isValid_Rockhopper()) return repliconList;

	    // Use appropriate file (local or server) to read in replicon information.
	    Scanner reader;
	    PrintWriter writer = null;
	    if (serverLength == fileLength) reader = new Scanner(repFile);
	    else {  // Read from server rather than local file. Create new local file.
		reader = new Scanner(new URL(Rockhopper_url + fileName).openStream());
		writer = new PrintWriter(repFile);
	    }
	    while (reader.hasNextLine()) {
		String repliconInfo = reader.nextLine();
		repliconList.add(repliconInfo);
		if (serverLength != fileLength) writer.println(repliconInfo);
	    }
	    reader.close();
	    if (serverLength != fileLength) writer.close();
	} catch (IOException e) {
	    Rockhopper.output("Could not process file " + genome_DIR + fileName + " and file " + Rockhopper_url + fileName + "\n");
	}
	return repliconList;
    }

    /**
     * Create dictionary of replicon tokens. Values are indices in ArrayList.
     */
    private HashMap<String,HashMap<Integer,Boolean>> createRepliconDictionary(ArrayList<String> repliconList) {
	HashMap<String,HashMap<Integer,Boolean>> repliconDictionary = new HashMap<String,HashMap<Integer,Boolean>>();
	for (int i=0; i<repliconList.size(); i++) {
	    String name = repliconList.get(i).split("\t")[0].toLowerCase();
	    for (int j=0; j<name.length()-HASH_SIZE+1; j++) {
		String key = name.substring(j, j+HASH_SIZE);
		if (!repliconDictionary.containsKey(key)) repliconDictionary.put(key, new HashMap<Integer,Boolean>());
		repliconDictionary.get(key).put(i, true);
	    }
	}
	return repliconDictionary;
    }

    /**
     * Retrieves a file from a server and stores it locally.
     * Returns true on success, false otherwise.
     */
    private boolean downloadFile(String serverFileName, String localFileName) {
	try {
	    //System.setProperty("java.net.preferIPv4Stack", "true");  // Windows firewall bug
	    Scanner reader = new Scanner(new URL(serverFileName).openStream());
	    PrintWriter writer = new PrintWriter(new File(localFileName));
	    while (reader.hasNextLine()) writer.println(reader.nextLine());
	    reader.close();
	    writer.close();
	    return true;
	} catch (IOException e) {
	    return false;
	}
    }



    /**********************************************
     **********   PUBLIC CLASS METHODS   **********
     **********************************************/

    public static String getRepliconDir(String repliconName) {
	for (int i=0; i<repliconName.length(); i++) {
	    char ch = repliconName.charAt(i);
	    if (!Character.isLetterOrDigit(ch)) repliconName = repliconName.replace(ch, '_');
	}
	return repliconName + "/";
    }



    /*************************************
     **********   MAIN METHOD   **********
     *************************************/

    public static void main(String[] args) {

    }

}

