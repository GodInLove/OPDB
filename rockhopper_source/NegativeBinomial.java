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

/**
 * Probability mass function of negative binomial distribution.
 * Based on saddle point algorithm found in "Fast and Accurate
 * Computation of Binomial Probabilities" by Catherine Loader, 2000.
 */
public class NegativeBinomial {

    /*****************************************
     **********   CLASS VARIABLES   **********
     *****************************************/

    private static double S0 = 0.08333333333333333333;
    private static double S1 = 0.00277777777777777778;
    private static double S2 = 0.00079365079365079365;
    private static double S3 = 0.00059523809523809524;
    private static double S4 = 0.00084175084175084175;
    private static double[] sfe = {0.0, 0.081061466795327258219670264, 0.041340695955409294093822081, 0.0276779256849983391487892927, 0.020790672103765093111522771, 0.0166446911898211921631948653, 0.013876128823070747998945727, 0.0118967099458917700950557241, 0.010411265261972096497478567, 0.0092554621827127329177286366, 0.008330563433362871256469318, 0.0075736754879518407949720242, 0.006942840107209529865664152, 0.0064089941880042070684396370, 0.005951370112758847735624416, 0.0055547335519628013710386899};



    /**********************************************
     **********   PUBLIC CLASS METHODS   **********
     **********************************************/

    public static double pmf(double k, double n, double p, boolean b) {
	if (p == 0.0) {
	    if (k == 0) return 1.0;
	    else return 0.0;
	} else if (p == 1.0) {
	    if (k == n) return 1.0;
	    else return 0.0;
	} else if (k == 0) {
	    return Math.exp(n * Math.log(1.0 - p));
	} else if (k == n) {
	    return Math.exp(n * Math.log(p));
	} else {
	    double lc = stirlerr(n) - stirlerr(k) - stirlerr(n - k) - bd0(k, n*p) - bd0(n-k, n*(1.0-p));
	    return p * Math.exp(lc) * Math.sqrt(n / (2.0*Math.PI*k*(n-k)));  // We multiply by "p" here
	}
    }

    public static double pmf(double k, double n, double p) {
	if (p == 0.0) {
	    if (k == 0) return 1.0;
	    else return 0.0;
	} else if (p == 1.0) {
	    if (k == n) return 1.0;
	    else return 0.0;
	} else if (k == 0) {
	    return Math.exp(n * Math.log(1.0 - p));
	} else if (k == n) {
	    return Math.exp(n * Math.log(p));
	} else {
	    double lc = stirlerr(n) - stirlerr(k) - stirlerr(n - k) - bd0(k, n*p) - bd0(n-k, n*(1.0-p));
	    return p * Math.exp(lc) * Math.sqrt(n / (2.0*Math.PI*k*(n-k)));  // We multiply by "p" here
	}
    }



    /***********************************************
     **********   PRIVATE CLASS METHODS   **********
     ***********************************************/

    /**
     * log(n!) - log(sqrt(2*pi*n)*(n/e)^n)
     */
    private static double stirlerr(double n) {
	if (n < 16) return sfe[(int)n];
	double nn = n*n;
	if (n > 500) return (S0 - S1/nn)/n;
	if (n > 80) return (S0 - (S1/S2/nn)/nn)/n;
	if (n > 35) return (S0 - (S1 - (S2-S3/nn)/nn)/nn)/n;
	return (S0 - (S1 - (S2 - (S3 - S4/nn)/nn)/nn)/nn)/n;
    }

    /**
     * Deviance term: k*lg(k/np) + np - k
     */
    private static double bd0(double k, double np) {
	if (Math.abs(k-np) < 0.1*(k+np)) {
	    double s = (k-np)*(k-np)/(k+np);
	    double v = (k-np)/(k+np);
	    double ej = 2*k*v;
	    int j = 1;
	    while (true) {
		ej = ej*v*v;
		double s1 = s+ej/(2*j+1);
		if (s1 == s) return s;
		s = s1;
		j += 1;
	    }

	}
	return (k*Math.log(k/np)+np-k);
    }



    /*************************************
     **********   MAIN METHOD   **********
     *************************************/

    public static void main(String[] args) {

	if (args.length < 3) {
	    System.err.println("\nUSAGE: java NegativeBinomial <k> <r> <p>" + "\n");
	    System.err.println("NegativeBinomial computes the probability mass function for the negative binomial distribution with parameters (r,p) at value k.\n");
	    System.exit(0);
	}

	int k = Integer.parseInt(args[0]);
	int r = Integer.parseInt(args[1]);
	double p = Double.parseDouble(args[2]);
	System.out.println(pmf(r-1, k+r-1, p));
    }

}
