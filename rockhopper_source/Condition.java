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
 * A Condition keeps track of the RNA-seq data in each
 * replicate experiment for a given condition.
 */
public class Condition {

    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private ArrayList<Replicate> replicates;
    private int partner;  // If there are no replicates, index of most similar other condition.
    private int minDiffExpressionLevel;  // Used for differentially expressed p-value correction



    /*****************************************
     **********   CLASS VARIABLES   **********
     *****************************************/

    private static long avgUpperQuartile;  // Avg upper quartile over all replicates in all conditions



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    public Condition() {
	this.replicates = new ArrayList<Replicate>();
	this.partner = -1;
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    /**
     * Return the Replicate at the specified index.
     */
    public Replicate getReplicate(int i) {
	if ((i >= 0) && (i < replicates.size()))
	    return replicates.get(i);
	return null;
    }

    /**
     * Adds a replicate to this Condition.
     */
    public void addReplicate(Replicate r) {
	this.replicates.add(r);
    }

    /**
     * Return the number of Replicates performed for this Condition.
     */
    public int numReplicates() {
	return replicates.size();
    }

    /**
     * Sets the partner of this Condition, i.e., the index of the other
     * condition that is most similar to this Condition.
     */
    public void setPartner(int partner) {
	this.partner = partner;
    }

    /**
     * Return the partner of this Condition, i.e., the index of the other
     * condition that is most similar to this Condition.
     */
    public int getPartner() {
	return this.partner;
    }

    /**
     * Estimates the minimum expression level, across all replicates,
     * necessary to compute a p-value for differential expression.
     */
    public int getMinDiffExpressionLevel() {
	return minDiffExpressionLevel;
    }

    /**
     * The background probability decreases as the number of reads increases.
     * We use this fact to estimate a minimum level of expression necessary
     * for us to be able to compute a p-value of differential expression.
     * For a given Condition, we ensure that all Replicates meet the
     * expression threshold. 
     * Currently, the expression threshold is set (based on anecdote and
     * experience) to 0.005.
     */
    public void setMinDiffExpressionLevel() {
	double THRESHOLD = 0.005;
	minDiffExpressionLevel = 0;
	while (true) {
	    boolean foundThreshold = true;
	    for (int i=0; i<replicates.size(); i++) {
		if (replicates.get(i).getBackgroundProb(minDiffExpressionLevel) >= THRESHOLD)
		    foundThreshold = false;
	    }
	    if (foundThreshold) return;
	    minDiffExpressionLevel++;
	}
    }

    /**
     * Returns a String representation of a this object.
     */
    public String toString() {
	StringBuilder sb = new StringBuilder();
	for (int i=0; i<replicates.size(); i++) {
	    sb.append(replicates.get(i).toString() + "\n");
	}
	return sb.toString();
    }



    /**********************************************
     **********   PUBLIC CLASS METHODS   **********
     **********************************************/

    /**
     * Set the average upper quartile gene expression over
     * all replicates in all all conditions.
     */
    public static void setAvgUpperQuartile(long avgUpperQuartile) {
	Condition.avgUpperQuartile = avgUpperQuartile;
    }

    /**
     * Returns the average upper quartile gene expression over 
     * all replicates in all conditions.
     */
    public static long getAvgUpperQuartile() {
	return Condition.avgUpperQuartile;
    }



    /**************************************************
     **********   PRIVATE INSTANCE METHODS   **********
     **************************************************/



    /*************************************
     **********   MAIN METHOD   **********
     *************************************/

    /**
     * The main method is used to test the methods of this class.
     */
    public static void main(String[] args) {

    }

}

