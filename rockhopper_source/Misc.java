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
import java.util.Random;
import java.util.Iterator;

public class Misc {

    /*************************************************
     **********   PRIVATE CLASS VARIABLES   **********
     *************************************************/

    private static Random rand;  // Random number generator



    /***********************************************
     **********   PUBLIC STATIC METHODS   **********
     ***********************************************/

    /**
     * Return the element in the ArrayList corresponding to 
     * the specified order statistic.
     */
    public static long select_Long(ArrayList<Long> a, int orderStatistic) {
	return select_Long(a, 0, a.size()-1, orderStatistic);
    }

    /**
     * Return the element in the ArrayList corresponding to 
     * the specified order statistic.
     */
    public static double select_Double(ArrayList<Double> a, int orderStatistic) {
	return select_Double(a, 0, a.size()-1, orderStatistic);
    }

    /**
     * Returns the correlation coefficient of two ArrayLists.
     */
    public static double correlation(ArrayList<Long> x, ArrayList<Long> y) {
	if (x.size() != y.size()) return 0.0;  // Lists have different lengths
	double meanX = mean(x);
	double meanY = mean(y);
	double sumXY = 0.0;
	double sumXX = 0.0;
	double sumYY = 0.0;
	for (int i=0; i<x.size(); i++) {
	    sumXY += (x.get(i) - meanX) * (y.get(i) - meanY);
	    sumXX += (x.get(i) - meanX) * (x.get(i) - meanX);
	    sumYY += (y.get(i) - meanY) * (y.get(i) - meanY);
	}
	if ((sumXX == 0.0) || (sumYY == 0.0)) return 0.0;
	return sumXY / (Math.sqrt(sumXX) * Math.sqrt(sumYY));
    }

    /**
     * Returns the mean of an ArrayList.
     */
    public static double mean(ArrayList<Long> a) {
	double sum = 0;
	for (Iterator<Long> iter = a.iterator(); iter.hasNext(); )
	    sum += iter.next();
	return sum / a.size();
    }



    /************************************************
     **********   PRIVATE STATIC METHODS   **********
     ************************************************/

    /**
     * Return the element in the ArrayList corresponding to
     * the specified order statistic.
     */
    private static long select_Long(ArrayList<Long> a, int lo, int hi, int orderStatistic) {
	if (lo >= hi) {
	    if (a.size() == 0) return 0;
	    return a.get(lo);
	}
	int p = partition(a, lo, hi);
	if (orderStatistic == p-lo+1) return a.get(p);
	else if (orderStatistic < p-lo+1) return select_Long(a, lo, p-1, orderStatistic);
	else return select_Long(a, p+1, hi, orderStatistic-(p-lo+1));
    }

    /**
     * Partitions "a" around random pivot element. 
     * Returns the index of the pivot.
     */
    private static int partition(ArrayList<Long> a, int lo, int hi) {
	if (rand == null) rand = new Random();
	swap(a, rand.nextInt(hi-lo+1)+lo, hi);  // Choose pivot randomly and place at end of "a"
	long x = a.get(hi);
	int i = lo-1;
	for (int j=lo; j<hi; j++) {
	    if (a.get(j) <= x) {
		i++;
		swap(a, i, j);
	    }
	}
	swap(a, i+1, hi);
	return i+1;
    }

    /**
     * Return the element in the ArrayList corresponding to
     * the specified order statistic.
     */
    private static double select_Double(ArrayList<Double> a, int lo, int hi, int orderStatistic) {
	if (lo == hi) return a.get(lo);
	int p = partitionDouble(a, lo, hi);
	if (orderStatistic == p-lo+1) return a.get(p);
	else if (orderStatistic < p-lo+1) return select_Double(a, lo, p-1, orderStatistic);
	else return select_Double(a, p+1, hi, orderStatistic-(p-lo+1));
    }

    /**
     * Partitions "a" around random pivot element. 
     * Returns the index of the pivot.
     */
    private static int partitionDouble(ArrayList<Double> a, int lo, int hi) {
	if (rand == null) rand = new Random();
	swap(a, rand.nextInt(hi-lo+1)+lo, hi);  // Choose pivot randomly and place at end of "a"
	double x = a.get(hi);
	int i = lo-1;
	for (int j=lo; j<hi; j++) {
	    if (a.get(j) <= x) {
		i++;
		swap(a, i, j);
	    }
	}
	swap(a, i+1, hi);
	return i+1;
    }

    /**
     * Swaps two elements in an ArrayList at the specified indices.
     */
    private static <E> void swap(ArrayList<E> a, int i, int j) {
	E temp = a.get(i);
	a.set(i, a.get(j));
	a.set(j, temp);
    }


    /*************************************
     **********   MAIN METHOD   **********
     *************************************/

    public static void main(String[] args) {

	System.err.println("\nClass Methods:\n");
	System.err.println("\t" + "select(ArrayList<Integer or Double> a, int orderStatistic)" + "\t" + "Computes the specified order statistic in the list.\n");
	System.err.println("\t" + "correlation(ArrayList<Integer> x, ArrayList<Integer> y)" + "\t" + "Computes the correlation coefficient of two lists.\n");
	System.err.println("\t" + "mean(ArrayList<Integer> a)" + "\t" + "Computes the mean of a list.\n");

    }

}

