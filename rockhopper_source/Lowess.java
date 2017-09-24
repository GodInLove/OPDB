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

public class Lowess {

    /**
     * Lowess smoother. Robust locally weighted regression.
     * The lowess function fits a nonparametric regression curve to a 
     * scatterplot. x and y must be of equal length. Returns the
     * estimatad smooth values of y.
     */
    public static ArrayList<Long> lowess(ArrayList<Long> x, ArrayList<Long> y) {
	return lowess(x, y, 0.1, 1);  // Default parameters
    }

    /**
     * Lowess smoother. Robust locally weighted regression.
     * The lowess function fits a nonparametric regression curve to a 
     * scatterplot. x and y must be of equal length. Returns the
     * estimatad smooth values of y.
     *
     * A larger smoothing span results in a smoother curve.
     * The number of iterations can be increased to improve robustness
     * at a cost to running time.
     */
    public static ArrayList<Long> lowess(ArrayList<Long> x, ArrayList<Long> y, double smoothingSpan, int iterations) {
	if (x.size() != y.size()) {
	    Rockhopper.output("Error - cannot compute lowess of two matrices with differing sizes.\n");
	    return new ArrayList<Long>();
	}

	// Create lists of nonzero elements
	ArrayList<Long> X1 = new ArrayList<Long>();
	ArrayList<Long> Y1 = new ArrayList<Long>();
	ArrayList<Long> X1_sorted = new ArrayList<Long>();
	for (int i=0; i<x.size(); i++) {
	    if (x.get(i) > 0) {
		X1.add(x.get(i));
		Y1.add(y.get(i));
		X1_sorted.add(x.get(i));
	    }
	}

	// "process" method creates a 2D array based on the size of X1.
	// If X1 is too large, then we run out of memory and crash.
	// Here we attempt to limit the size of the dimensions, i.e., the X1 array.
	int power = 0;
	int threshold = 0;
	while ((X1.size() >= 10000) && (Runtime.getRuntime().maxMemory() < 1300000000) && (power < 31)) {
	    threshold = (int)Math.pow(2, power);
	    X1.clear();
	    Y1.clear();
	    X1_sorted.clear();
	    for (int i=0; i<x.size(); i++) {
		if (x.get(i) > threshold) {
		    X1.add(x.get(i));
		    Y1.add(y.get(i));
		    X1_sorted.add(x.get(i));
		}
	    }
	    power++;
	}

	// Sort X1 in reverse order
	Collections.sort(X1_sorted);
	Collections.reverse(X1_sorted);

	ArrayList<Long> h = new ArrayList<Long>();
	ArrayList<Long> dist = new ArrayList<Long>();
	for (int i=0; i<X1.size(); i++) {
	    Long p = X1.get(i);
	    dist.clear();
	    double fraction = (X1_sorted.size() - X1_sorted.indexOf(p)) / (double)(X1_sorted.size());
	    int r = (int)(Math.ceil(Math.min(2.0*fraction, 0.4)*X1.size()));
	    for (int j=0; j<X1.size(); j++) dist.add(Math.abs(p.longValue() - X1.get(j).longValue()));
	    Collections.sort(dist);
	    if (r < dist.size()) h.add(dist.get(r));
	    else h.add(dist.get(dist.size()-1));
	}

	double[][] w = process(X1, h);
	cube(w);
	oneMinus(w);
	cube(w);

	double[] yest = new double[X1.size()];
	double[] delta = new double[X1.size()];

	for (int z=0; z<iterations; z++) {
	    for (int i=0; i<X1.size(); i++) {
		ArrayList<Double> weights = new ArrayList<Double>();
		for (int j=0; j<w.length; j++) weights.add(w[j][i]);
		double b1 = sum(product(weights, Y1));
		double b2 = sum(product(product(weights, Y1), X1));
		double A1 = sum(weights);
		double A2 = sum(product(weights, X1));
		double A3 = A2;
		double A4 = sum(product(product(weights, X1), X1));
		double[][] b = {{b1}, {b2}};
		double[][] A = {{A1, A2}, {A3, A4}};
		Matrix b_Matrix = new Matrix(b);
		Matrix A_Matrix = new Matrix(A);
		double[][] beta = A_Matrix.solve(b_Matrix).getArray();  // Linear algebra solver
		yest[i] = beta[0][0] + (beta[1][0] * X1.get(i));
	    }
	    double[] residuals = new double[Y1.size()];
	    ArrayList<Double> residuals_absValue = new ArrayList<Double>();
	    for (int i=0; i<Y1.size(); i++) {
		residuals[i] = Y1.get(i) - yest[i];
		residuals_absValue.add(Math.abs(Y1.get(i) - yest[i]));
	    }
	    double s = 1.0;
	    if (residuals.length > 0) s = Misc.select_Double(residuals_absValue, 1+(residuals.length/2));

	    for (int i=0; i<delta.length; i++) {
		delta[i] = residuals[i] / (6.0 * s);
		if (delta[i] < -1) delta[i] = -1.0;
		if (delta[i] > 1) delta[i] = 1.0;
	    }
	    oneMinus(product(delta, delta));
	    delta = product(delta, delta);
	}

	double minPositive = Double.MAX_VALUE;
	for (int i=0; i<yest.length; i++) {
	    if (yest[i] > threshold) {
		if (yest[i] < minPositive)
		    minPositive = yest[i];
	    }
	}

	ArrayList<Long> yest_Long = new ArrayList<Long>();
	int j=0;
	for (int i=0; i<x.size(); i++) {
	    if (x.get(i) <= threshold)
		yest_Long.add((long)minPositive);
	    else {
		yest_Long.add((long)yest[j]);
		j++;
	    }
	}
	return yest_Long;
    }



    /************************************************
     **********   PRIVATE STATIC METHODS   **********
     ************************************************/

    /**
     * clip(abs(([x]-transpose([x]))/h),0.0,1.0)
     */
    private static double[][] process(ArrayList<Long> x, ArrayList<Long> h) {
	double[][] w = new double[x.size()][x.size()];
	for (int i=0; i<x.size(); i++) {
	    for (int j=0; j<x.size(); j++) {
		double value = Math.abs(((double)(x.get(i) - x.get(j)))/h.get(i));
		if (value < 0.0) value = 0.0;
		if (value > 1.0) value = 1.0;
		w[j][i] = value;
	    }
	}
	return w;
    }

    /**
     * Cube every element in a 2D array.
     */
    private static void cube(double[][] w) {
	for (int i=0; i<w.length; i++)
	    for (int j=0; j<w[0].length; j++)
		w[i][j] = w[i][j] * w[i][j] * w[i][j];
    }

    /**
     * w = 1 - w
     */
    private static void oneMinus(double[][] w) {
	for (int i=0; i<w.length; i++)
	    for (int j=0; j<w[0].length; j++)
		w[i][j] = 1.0 - w[i][j];
    }

    /**
     * w = 1 - w
     */
    private static void oneMinus(double[] w) {
	for (int i=0; i<w.length; i++)
	    w[i] = 1.0 - w[i];
    }

    /**
     * Return a new ArrayList where each element is the product of the corresponding
     * elements in the two specified ArrayLists.
     */
    private static ArrayList<Double> product(ArrayList<Double> x, ArrayList<Long> y) {
	ArrayList<Double> result = new ArrayList<Double>();
	for (int i=0; i<x.size(); i++)
	    result.add(x.get(i) * y.get(i));
	return result;
    }

    /**
     * Returns a new array where each element is the product of the corresponding
     * elements in the two specified arrays.
     */
    private static double[] product(double[] x, double[] y) {
	double[] result = new double[x.length];
	for (int i=0; i<x.length; i++)
	    result[i] = x[i] * y[i];
	return result;
    }

    /**
     * Returns the sum of all elements in an ArrayList.
     */
    private static double sum(ArrayList<Double> x) {
	double sum = 0.0;
	for (int i=0; i<x.size(); i++) sum += x.get(i);
	return sum;
    }

    /**
     * Output to Stdout a 2D array.
     */
    private static void print(double[][] w) {
	Rockhopper.output("\n");
	for (int i=0; i<w.length; i++) {
	    for (int j=0; j<w[0].length; j++) {
		Rockhopper.output("\t" + w[i][j]);
	    }
	    Rockhopper.output("\n");
	}
	Rockhopper.output("\n");
    }



    /*************************************
     **********   MAIN METHOD   **********
     *************************************/
    /**
     * The main method is used to test the methods of this class.
     */
    public static void main(String[] args) {
	ArrayList<Long> a = new ArrayList<Long>();
	a.add((long)1);
	a.add((long)2);
	a.add((long)5);
	a.add((long)8);
	a.add((long)7);
	a.add((long)3);
	a.add((long)20);
	a.add((long)13);
	a.add((long)10);
	ArrayList<Long> b = new ArrayList<Long>();
	b.add((long)4);
	b.add((long)9);
	b.add((long)6);
	b.add((long)0);
	b.add((long)1);
	b.add((long)5);
	b.add((long)15);
	b.add((long)11);
	b.add((long)18);
	System.out.println("\n" + lowess(a, b) + "\n");


    }

}

