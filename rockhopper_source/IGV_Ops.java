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
import java.io.PrintWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import org.broad.igv.feature.genome.GenomeImporter;

public class IGV_Ops {

    /**
     * Creates a genome registry for IGV. File is named genomes.txt.
     */
    public static boolean createGenomeRegistry(File genomeRegistryFile, ArrayList<Genome> genomes) {
	try {
	    PrintWriter writer = new PrintWriter(genomeRegistryFile);
	    writer.println("<Server-Side Genome List>");
	    for (Genome genome : genomes)
		writer.println(getGenomeDisplayName(genome.getName()) + "\t" + genome.getBaseFileName() + ".genome" + "\t" + genome.getName().replace(' ', '_'));
	    writer.close();
	    return true;
	} catch (FileNotFoundException ex) {
	    return false;
	}
    }

    /**
     * Creates a genome archive file for use by IGV. This method
     * utilizes IGV code.
     */
    public static boolean createGenomeArchive(String genomeArchiveFileName, String genomeName, String fastaFileName, String geneFileName) {
	try {
	    (new GenomeImporter()).createGenomeArchive(new File(genomeArchiveFileName), genomeName.replace(' ', '_'), getGenomeDisplayName(genomeName), fastaFileName, new File(geneFileName), null, null);
	    return true;
	} catch (IOException ex) {
	    return false;
	}
    }

    /**
     * Create a batch script to execute when IGV launches to load
     * the sequencing reads data.
     */
    public static boolean createBatchFileToLoadData(String batchFileName, String genomeName, String dataFiles, String rockhopperFiles) {
	try {
	    PrintWriter writer = new PrintWriter(new File(batchFileName));
	    writer.println("new");
	    writer.println("genome " + genomeName.replace(' ', '_'));
	    if (dataFiles.length() > 0) writer.println("load " + dataFiles);
	    if (rockhopperFiles.length() > 0) writer.println("load " + rockhopperFiles);
	    writer.close();
	    return true;
	} catch (FileNotFoundException ex) {
	    return false;
	}
    }

    /**
     * Returns the name of the genome to display in IGV's drop-down menu.
     */
    private static String getGenomeDisplayName(String genomeName) {
	genomeName = genomeName.trim();
	int indexOfSpace = genomeName.indexOf(' ');
	if (indexOfSpace >= 0) genomeName = genomeName.charAt(0) + "." + genomeName.substring(indexOfSpace);
	return genomeName;
    }

}

