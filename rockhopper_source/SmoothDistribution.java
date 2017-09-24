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
import java.util.Collections;

public class SmoothDistribution {

    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private ArrayList<Integer> data;
    private double bandwidth;
    private int BIN_SIZE;
    private int minimum;
    private int maximum;
    private double pseudocount;

    private ArrayList<Integer> histogram;
    private ArrayList<Double> histogramNormalized;
    private ArrayList<Integer> values;
    private ArrayList<Double> smoothed;
    private ArrayList<Double> smoothedNormalized;



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    public SmoothDistribution(ArrayList<Integer> data) {
	this(data, 1.0, 1);
    }

    public SmoothDistribution(ArrayList<Integer> data, double bandwidth, int BIN_SIZE) {
	this(data, bandwidth, BIN_SIZE, Collections.min(data), Collections.max(data));
    }

    public SmoothDistribution(ArrayList<Integer> data, double bandwidth, int BIN_SIZE, int minimum, int maximum) {
	this.data = data;
	this.minimum = minimum;
	this.maximum = maximum;
	this.bandwidth = bandwidth;
	this.BIN_SIZE = BIN_SIZE;
	this.pseudocount = 0.0;
	generateHistogram();
	smoothData();
	normalize();
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    /**
     * Get smoothed value at point "x".
     */
    public double getSmoothedValue(double x) {
	if (x < this.minimum) return this.pseudocount;
	if (x > this.maximum) return this.pseudocount;
	int index = binarySearch(values, x, 0, values.size()-1);
	// If "x" is not in our discrete distribution, the index may be off by 1.
	// For example, if our distribution includes 4,5,6,7 and x=5.3, we
	// may return index 1 (corresponding to 5) or index 2 (corresponding to 6).
	if (index > 0) {
	    if (Math.abs(x - values.get(index)) > Math.abs(x - values.get(index-1)))
		return Math.max(smoothedNormalized.get(index-1), this.pseudocount);
	}
	if (index < smoothedNormalized.size()-1) {
	    if (Math.abs(x - values.get(index)) > Math.abs(x - values.get(index+1)))
		return Math.max(smoothedNormalized.get(index+1), this.pseudocount);
	}
	return Math.max(smoothedNormalized.get(index), this.pseudocount);
    }

    /**
     * Set the pseudocount for this distribution.
     */
    public void setPseudocount(double pseudocount) {
	this.pseudocount = pseudocount;
    }

    /**
     * Returns a String representation of this smoothed distribution.
     */
    public String toString() {
	StringBuilder sb = new StringBuilder();
	sb.append("VALUES" + "\t" + "SMOOTHED" + "\t" + "NORMALIZED" + "\n");
	for (int i=0; i<smoothedNormalized.size(); i++) {
	    sb.append(values.get(i) + "\t" + histogram.get(i) + "\t" + histogramNormalized.get(i) + "\t" + smoothedNormalized.get(i) + "\n");
	}
	return sb.toString();
    }



    /**************************************************
     **********   PRIVATE INSTANCE METHODS   **********
     **************************************************/

    /**
     * Generate histogram of data.
     */
    private void generateHistogram() {
	this.histogram = new ArrayList<Integer>();
	this.histogramNormalized = new ArrayList<Double>();
	for (int i=this.minimum; i<=maximum; i+=BIN_SIZE) {
	    histogram.add(0);
	    histogramNormalized.add(0.0);
	}
	for (int i=0; i<data.size(); i++) {
	    if (data.get(i) < minimum) histogram.set(0, histogram.get(0) + 1);
	    else if (data.get(i) > maximum) histogram.set(histogram.size()-1, histogram.get(histogram.size()-1) + 1);
	    else histogram.set((data.get(i)-minimum)/BIN_SIZE, histogram.get((data.get(i)-minimum)/BIN_SIZE) + 1);
	}
	for (int i=0; i<histogram.size(); i++) histogramNormalized.set(i, histogram.get(i) / (double)data.size());
    }

    /**
     * Generate smooth distibution of data (based on Epanechnikov kernel).
     */
    private void smoothData() {
	values = new ArrayList<Integer>();
	smoothed = new ArrayList<Double>();
	int i = this.minimum;
	while (i <= this.maximum) {  // Compute weighted value for index i
	    double sum = 0.0;
	    for (int j=0; j<data.size(); j++) {
		if (Math.abs(i - data.get(j)) <= this.bandwidth) {
		    double u = (i - data.get(j)) / this.bandwidth;
		    sum += (3.0/4.0) * (1.0 - u*u);
		}
	    }
	    values.add(i);
	    smoothed.add(sum / (data.size() * this.bandwidth));
	    i += this.BIN_SIZE;
	}
    }

    /**
     * Determine pseudocounts and normalize distribution so it sums to 1.0.
     */
    private void normalize() {
	smoothedNormalized = new ArrayList<Double>();
	double sum = 0.0;
	for (int i=0; i<smoothed.size(); i++) sum += smoothed.get(i);
	double min = Double.MAX_VALUE;
	for (int i=0; i<smoothed.size(); i++) {
	    double normalizedValue = smoothed.get(i) / sum;
	    smoothedNormalized.add(normalizedValue);
	    if ((normalizedValue > 0) && (normalizedValue < min)) min = normalizedValue;
	}
	this.pseudocount = min / 10.0;
    }

    /**
     * Perform a binary search for the specified value.
     */
    private int binarySearch(ArrayList<Integer> a, double value, int lo, int hi) {
	if (hi <= lo) return lo;
	int mid = lo + (hi-lo)/2;
	if (a.get(mid) > value) return binarySearch(a, value, lo, mid-1);
	if (a.get(mid) < value) return binarySearch(a, value, mid+1, hi);
	return mid;
    }



    /*************************************
     **********   MAIN METHOD   **********
     *************************************/

    public static void main(String[] args) {
	System.err.println("\nThe SmoothDistribution application cannot be executed from the command line. It must be instantiated from another Java application. It takes a set of data and generates a smoothed version of the data distribution based on the Epanechnikov kernel. Smoothing occurs at each point over an interval [-bandwidth,bandwidth]. The BIN_SIZE indicates how big each BIN_SIZE should be.\n");
    }

}
