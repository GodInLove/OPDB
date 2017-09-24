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
import java.util.Scanner;
import java.io.PrintWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * An instance of the Genome class represents a genome and its
 * annotation, including the genome sequence, protein-coding
 * genes, and RNA genes in the genome.
 */
public class Genome {

    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private String ID;
    private String name;
    private String genome;  // We add a '?' character so it is 1-indexed
    private String formalGenomeName = "";  // First token of FASTA header line, e.g., gi|49175990|ref|NC_000913.2|
    private String baseFileName = "";  // Path and name of genome file sans .fna extension
    private ArrayList<Gene> codingGenes = new ArrayList<Gene>();
    private ArrayList<Gene> rnas = new ArrayList<Gene>();
    private ArrayList<Gene> genes = new ArrayList<Gene>();  // Merged coding genes and RNAs



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    /**
     * Constructs a new Genome object based on files in the specified
     * directory (e.g., genome.fna, genes.ptt, rna.rnt).
     */
    public Genome(String directory) {
	String genomeFileName = null;
	String geneFileName = null;
	String rnaFileName = null;

	if (!directory.endsWith("/")) directory += "/";
	File dir = new File(directory);
	String[] files = dir.list();
	for (int i=0; i<files.length; i++) {
	    if (files[i].endsWith(".fna")) genomeFileName = files[i];
	    if (files[i].endsWith(".ptt")) geneFileName = files[i];
	    if (files[i].endsWith(".rnt")) rnaFileName = files[i];
	}

	if (genomeFileName != null) readInGenome(directory + genomeFileName);
	if (geneFileName != null) this.codingGenes = readInGenes(directory + geneFileName, "ORF");
	if (rnaFileName != null) this.rnas = readInGenes(directory + rnaFileName, "RNA");
	this.genes = mergeGenes(this.codingGenes, this.rnas);
	if (genomeFileName != null) create_GFF_file(this.genes, directory + genomeFileName, this.genome.length()-1);
	if (genomeFileName != null) baseFileName = directory + genomeFileName.substring(0, genomeFileName.length()-4);
    }

   /**
     * Constructs a new Genome object based on the specified files 
     * (e.g., genome.fna, genes.ptt, rna.rnt). If any of the file
     * names are empty, then they are simply ignored.
     */
    public Genome(String genomeFileName, String geneFileName, String rnaFileName) {
	if (genomeFileName != null) readInGenome(genomeFileName);
	if (geneFileName != null) this.codingGenes = readInGenes(geneFileName, "ORF");
	if (rnaFileName != null) this.rnas = readInGenes(rnaFileName, "RNA");
	this.genes = mergeGenes(this.codingGenes, this.rnas);
	if (genomeFileName != null) create_GFF_file(this.genes, genomeFileName, this.genome.length()-1);
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    /**
     * Returns the size in nucleotides of the genome
     * sequence. Note: the returned size will be 1
     * bigger than the actual genome size since we
     * add a character to the beginning of the 
     * genome sequence to ensure 1-indexing.
     */
    public int size() {
	return genome.length();
    }

    /**
     * Returns the ID of the genome.
     */
    public String getID() {
	return ID;
    }

    /**
     * Returns the name of the genome.
     */
    public String getName() {
	return name;
    }

    /**
     * Returns the path and name of the genome file sans .fna extension.
     */
    public String getBaseFileName() {
	return baseFileName;
    }

    /**
     * Returns a subsequence of the (1-indexed) genome
     * as specified by the two coordinates, inclusive.
     */
    public String getSeq(int start, int stop) {
	return genome.substring(Math.max(start, 1), Math.min(stop+1, genome.length()));
    }

    /**
     * Return the number of genes in the genome.
     */
    public int numGenes() {
	return genes.size();
    }

    /**
     * Returns the formal genome name, e.g., gi|49175990|ref|NC_000913.2|
     */
    public String getFormalGenomeName() {
	return this.formalGenomeName;
    }

    /**
     * Return the Gene at the specified index.
     */
    public Gene getGene(int i) {
	if ((i >= 0) && (i < genes.size()))
	    return genes.get(i);
	return null;
    }

    /**
     * Return the collection of protein coding genes.
     */
    public ArrayList<Gene> getCodingGenes() {
	return this.codingGenes;
    }

    public void addPredictedRNAs(ArrayList<Gene> predictedRNAs) {
	this.genes = mergeGenes(this.genes, predictedRNAs);
    }

    /**
     * Return the collection of all genes.
     */
    public ArrayList<Gene> getGenes() {
	return this.genes;
    }

    /**
     * Returns a String[] representing the annotation
     * (gene, rRNA, tRNA, RNA, "") for each
     * coordinate in the genome on the specified strand.
     */
    public String[] getAnnotations(char strand) {
	String[] annotation = new String[genome.length()];
	for (int i=0; i<annotation.length; i++) annotation[i] = "";  // Initialize
	for (int i=0; i<codingGenes.size(); i++) {
	    Gene g = codingGenes.get(i);
	    if (g.getStrand() == strand) {
		for (int j=g.getFirst(); j<=g.getLast(); j++) annotation[j] = g.getName();
	    }
	}
	for (int i=0; i<rnas.size(); i++) {
	    Gene g = rnas.get(i);
	    if (g.getStrand() == strand) {
		for (int j=g.getFirst(); j<=g.getLast(); j++) {
		    if (g.getProduct().contains("tRNA")) annotation[j] = "tRNA";
		    else if ((g.getProduct().contains("rRNA")) || (g.getProduct().contains("ribosomal"))) annotation[j] = "rRNA";
		    else annotation[j] = "RNA";  // Miscellaneous RNA
		}
	    }
	}
	return annotation;
    }

    /**
     * Return a String representation of all genes in the genome.
     */
    public String genesToString(ArrayList<Condition> conditions, String[] labels) {
	StringBuilder sb = new StringBuilder();
	sb.append("Transcription Start" + "\t" + "Translation Start" + "\t" + "Translation Stop" + "\t" + "Transcription Stop" + "\t" + "Strand" + "\t" + "Name" + "\t" + "Synonym" + "\t" + "Product");
	for (int j=0; j<conditions.size(); j++) {
	    String conditionName = "" + (j+1);
	    if ((labels != null) && (labels.length == conditions.size())) conditionName = labels[j];
	    if (Rockhopper.verbose) {  // Verbose output
		if (conditions.get(j).numReplicates() == 1) {
		    sb.append("\t" + "Raw Counts " + conditionName);
		    sb.append("\t" + "Normalized Counts " + conditionName);
		} else {
		    for (int k=0; k<conditions.get(j).numReplicates(); k++)
			sb.append("\t" + "Raw Counts " + conditionName + " Replicate " + (k+1));
		    for (int k=0; k<conditions.get(j).numReplicates(); k++)
			sb.append("\t" + "Normalized Counts " + conditionName + " Replicate " + (k+1));
		}
		sb.append("\t" + "RPKM " + conditionName);
	    }
	    sb.append("\t" + "Expression " + conditionName);  // Label experiment conditions
	}
	int numQvalues = 0;
	for (int x=0; x<conditions.size()-1; x++) {  // First in pair
	    for (int y=x+1; y<conditions.size(); y++) {  // Second in pair
		if (Rockhopper.verbose) {  // Verbose output
		    if ((labels != null) && (labels.length == conditions.size()))  // Use labels
			sb.append("\t" + "pValue " + labels[x] + " vs " + labels[y]);
		    else  // Do not use labels
			sb.append("\t" + "pValue " + (x+1) + " vs " + (y+1));
		}
		if ((numGenes() > 0) && (genes.get(0).hasQvalue(numQvalues))) {
		    if ((labels != null) && (labels.length == conditions.size()))  // Use labels
			sb.append("\t" + "qValue " + labels[x] + " vs " + labels[y]);
		    else  // Do not use labels
			sb.append("\t" + "qValue " + (x+1) + " vs " + (y+1));
		}
		numQvalues++;
	    }
	}
	sb.append("\n");
	for (int i=0; i<numGenes(); i++) {
	    sb.append(genes.get(i).toString() + genes.get(i).expressionToString() + "\n");
	}
	return sb.toString();
    }



    /**************************************************
     **********   PRIVATE INSTANCE METHODS   **********
     **************************************************/

    /* Sets the "genome" instance variable with the FASTA formatted
     * sequence found in the specified file. Also, sets the "ID"
     * instance variable and the "name" instance variable based
     * on the FASTA header line.
     */
    private void readInGenome(String fileName) {
	this.ID = "";
	this.name = "";
	this.genome = "";

        StringBuilder sequence = new StringBuilder();
        try {
            Scanner reader = new Scanner(new File(fileName));

	    // Extract header info
            String header = reader.nextLine();  // Header line of FASTA file
            if ((header.length() == 0) || (header.charAt(0) != '>')) {
                Rockhopper.output("Error - first line of file " + fileName + " is not in FASTA format.\n");
                return;
            }
	    String[] parse_header = header.split("\\|");
	    if (parse_header.length >= 5) {
		String[] parse_ID = parse_header[3].split("\\.");
		this.ID = parse_ID[0];
		String[] parse_name = parse_header[4].split(",");
		this.name = parse_name[0].trim();
	    }

	    // Read in sequence
	    sequence.append('?');
            while (reader.hasNext()) {  // continue until we reach end of file
                sequence.append(reader.nextLine());
            }
            reader.close();
	    this.genome = sequence.toString();
        } catch (FileNotFoundException e) {
            Rockhopper.output("Error - the file " + fileName + " could not be found and opened.\n");
        }
    }

    /**
     * Reads in a file of genes (either *.ptt or *.rnt) and returns
     * an ArrayList of gene objects.
     */
    private ArrayList<Gene> readInGenes(String fileName, String type) {
	ArrayList<Gene> listOfGenes = new ArrayList<Gene>();
        try {
            Scanner reader = new Scanner(new File(fileName));
	    for (int i=0; i<3; i++) reader.nextLine();  // Ignore 3 header lines
	    while (reader.hasNext()) {  // Continue until end of file
		listOfGenes.add(new Gene(reader.nextLine(), type));  // Create new gene
	    }
            reader.close();
        } catch (FileNotFoundException e) {
            Rockhopper.output("Error - the file " + fileName + " could not be found and opened.\n");
        }
	return listOfGenes;
    }

    /**
     * Given two lists of Genes, creates a new list of Genes containing
     * all the Genes from the two specified lists combined and sorted
     * by their coordinates.
     */
    private ArrayList<Gene> mergeGenes(ArrayList<Gene> genes1, ArrayList<Gene> genes2) {
	ArrayList<Gene> listOfGenes = new ArrayList<Gene>();
	int i=0;  // Index for genes1 
	int j=0;  // Index for genes2
	while ((i < genes1.size()) && (j < genes2.size())) {
	    if (genes1.get(i).getMinCoordinate() < genes2.get(j).getMinCoordinate()) {
		listOfGenes.add(genes1.get(i));
		i++;
	    } else {
		listOfGenes.add(genes2.get(j));
		j++;
	    }
	}
	while (i < genes1.size()) {
	    listOfGenes.add(genes1.get(i));
	    i++;
	}
	while (j < genes2.size()) {
	    listOfGenes.add(genes2.get(j));
	    j++;
	}
	return listOfGenes;
    }

    /**
     * Based on a list of genes and a FASTA genome file, 
     * create a GFF file that can be read by genome browser.
     * The newly created file is output to the same directory
     * as the FASTA genome file and has the same name except
     * with the extension *.gff.
     */
    private void create_GFF_file(ArrayList<Gene> genes, String genomeFileName, int genomeSize) {
	String GFF_fileName = genomeFileName.substring(0, genomeFileName.length()-4) + ".gff";
	try {

	    // Get info from FASTA genome file
	    Scanner reader = new Scanner(new File(genomeFileName));
	    String line = reader.nextLine();  // Header line
	    this.formalGenomeName = line.split("\\s+")[0].substring(1);
	    reader.close();

	    // Create GFF file
	    PrintWriter writer = new PrintWriter(new File(GFF_fileName));
	    writer.println("track name=Genes color=255,0,255");
	    writer.println("##gff-version 3");  // Header line
	    writer.println("##sequence-region " + formalGenomeName + " 1 " + genomeSize);
	    for (int i=0; i<genes.size(); i++) {
		Gene g = genes.get(i);
		String type = "RNA";
		int start = g.getMinCoordinate();
		int stop = g.getMaxCoordinate();
		if (g.isORF()) {
		    type = "Coding gene";
		    start = Math.min(g.getStart(), g.getStop());
		    stop = Math.max(g.getStart(), g.getStop());
		}
		writer.println(formalGenomeName + "\t" + "RefSeq" + "\t" + type + "\t" + start + "\t" + stop + "\t" + "." + "\t" + g.getStrand() + "\t" + "." + "\t" + "name=" + g.getName() + ";" + "product=" + "\"" + g.getProduct() + "\"");
	    }
	    writer.close();
	} catch (IOException e) {
	    Rockhopper.output("Error - could not create GFF file " + GFF_fileName + "\n");
	    return;
	}
    }



    /*************************************
     **********   MAIN METHOD   **********
     *************************************/

    /**
     * The main method is used to test the methods of this class.
     */
    public static void main(String[] args) {

	if (args.length == 0) {
	    System.err.println("\nThe Genome application must be invoked with a command line argument corresponding to the name of a directory containing genome files (e.g., *.fna, *.ptt, *.rnt).\n");
	    System.exit(0);
	}

	Genome g = new Genome(args[0]);
	System.out.println();
	System.out.println("Genome ID:\t\t" + g.ID);
	System.out.println("Genome name:\t\t" + g.name);
	System.out.println("Size of genome is:\t" + g.genome.length());
	System.out.println("First 10 nucleotides:\t" + g.genome.substring(0, 10));
	System.out.println("Final 10 nucleotides:\t" + g.genome.substring(g.genome.length()-10));
	System.out.println("Number of ORFs:\t\t" + g.codingGenes.size());
	System.out.println("Number of RNAs:\t\t" + g.rnas.size());
	System.out.println("Number of all genes:\t" + g.genes.size());
	System.out.println();
    }

}

