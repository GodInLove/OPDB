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

public class Matrix {

    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private double[][] M;
    private int rows, cols;  // Number or rows and columns



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    public Matrix (int rows, int cols) {
	M = new double[rows][cols];
	this.rows = rows;
	this.cols = cols;
    }

    public Matrix (double[][] M) {
	this.M = M;
	rows = M.length;
	cols = M[0].length;
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    /**
     * Returns the number of rows in this Matrix.
     */
    public int getRows () {
	return rows;
    }

    /**
     * Returns the number of columns in this Matrix.
     */
    public int getCols () {
	return cols;
    }

    /**
     * Returns the 2D array representing this Matrix.
     */
    public double[][] getArray () {
	return M;
    }

    /**
     * Returns a copy of the 2D array representing this Matrix.
     */

    public double[][] getArrayCopy () {
	double[][] M2 = new double[rows][cols];
	for (int i = 0; i < rows; i++)
	    for (int j = 0; j < cols; j++)
		M2[i][j] = M[i][j];
	return M2;
    }

    /**
     * Return a submatrix of this Matrix.
     */
    public Matrix getMatrix (int[] r, int x, int y) {
	Matrix mtrx = new Matrix(r.length,y-x+1);
	double[][] M2 = mtrx.getArray();
	for (int i = 0; i < r.length; i++)
	    for (int j = x; j <= y; j++)
		M2[i][j-x] = M[r[i]][j];
	return mtrx;
    }

    /**
     * Return a submatrix of this Matrix.
     */
    public Matrix getMatrix (int x1, int x2, int y1, int y2) {
	Matrix mtrx = new Matrix(x2-x1+1,y2-y1+1);
	double[][] M2 = mtrx.getArray();
	for (int i = x1; i <= x2; i++)
	    for (int j = y1; j <= y2; j++)
		M2[i-x1][j-y1] = M[i][j];
	return mtrx;
    }

    public Matrix solve (Matrix B) {
	if (rows == cols) return (new LUDecomposition(this)).solve(B);
	else return (new QRDecomposition(this)).solve(B);
    }



    /***************************************
     **********   CLASS METHODS   **********
     ***************************************/

    public static double hypot(double x, double y) {
	double r;
	if (Math.abs(x) > Math.abs(y)) {
	    r = y/x;
	    r = Math.abs(x)*Math.sqrt(1+r*r);
	} else if (y != 0) {
	    r = x/y;
	    r = Math.abs(y)*Math.sqrt(1+r*r);
	} else {
	    r = 0.0;
	}
	return r;
    }

}



/*****************************************
 **********   LU MATRIX CLASS   **********
 *****************************************/

class LUDecomposition {

    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private double[][] LU;
    private int rows;
    private int cols;
    private int signOfPivot; 
    private int[] pivot;



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    public LUDecomposition (Matrix A) {

	LU = A.getArrayCopy();
	rows = A.getRows();
	cols = A.getCols();
	pivot = new int[rows];
	for (int i = 0; i < rows; i++) pivot[i] = i;
	signOfPivot = 1;
	double[] LUrowi;
	double[] LUcolj = new double[rows];
	
	for (int j = 0; j < cols; j++) {	    
	    for (int i = 0; i < rows; i++) LUcolj[i] = LU[i][j];
	    for (int i = 0; i < rows; i++) {
		LUrowi = LU[i];
		int kmax = Math.min(i,j);
		double s = 0.0;
		for (int k = 0; k < kmax; k++) s += LUrowi[k]*LUcolj[k];
		LUrowi[j] = LUcolj[i] -= s;
	    }
	    int p = j;
	    for (int i = j+1; i < rows; i++) {
		if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p])) p = i;
	    }
	    if (p != j) {
		for (int k = 0; k < cols; k++) {
		    double t = LU[p][k];
		    LU[p][k] = LU[j][k];
		    LU[j][k] = t;
		}
		int k = pivot[p]; pivot[p] = pivot[j]; pivot[j] = k;
		signOfPivot = -signOfPivot;
	    }
	    if (j < rows & LU[j][j] != 0.0) {
		for (int i = j+1; i < rows; i++) LU[i][j] /= LU[j][j];
	    }
	}
    }

    

    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    public Matrix solve (Matrix B) {

	int cols2 = B.getCols();
	Matrix mtrx = B.getMatrix(pivot,0,cols2-1);
	double[][] M2 = mtrx.getArray();
	
	for (int k = 0; k < cols; k++)
	    for (int i = k+1; i < cols; i++)
		for (int j = 0; j < cols2; j++)
		    M2[i][j] -= M2[k][j]*LU[i][k];

	for (int k = cols-1; k >= 0; k--) {
	    for (int j = 0; j < cols2; j++)
		M2[k][j] /= LU[k][k];
	    for (int i = 0; i < k; i++)
		for (int j = 0; j < cols2; j++)
		    M2[i][j] -= M2[k][j]*LU[i][k];
	}
	return mtrx;
    }

}



/*****************************************
 **********   QR MATRIX CLASS   **********
 *****************************************/

class QRDecomposition {

    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private double[][] QR;
    private int rows;
    private int cols;
    private double[] Rdiag;



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    public QRDecomposition (Matrix A) {
	QR = A.getArrayCopy();
	rows = A.getRows();
	cols = A.getCols();
	Rdiag = new double[cols];
	for (int k = 0; k < cols; k++) {
	    double nrm = 0;
	    for (int i = k; i < rows; i++) nrm = Matrix.hypot(nrm,QR[i][k]);
	    if (nrm != 0.0) {
		if (QR[k][k] < 0) nrm = -nrm;
		for (int i = k; i < rows; i++) QR[i][k] /= nrm;
		QR[k][k] += 1.0;
		for (int j = k+1; j < cols; j++) {
		    double s = 0.0; 
		    for (int i = k; i < rows; i++) s += QR[i][k]*QR[i][j];
		    s = -s/QR[k][k];
		    for (int i = k; i < rows; i++) QR[i][j] += s*QR[i][k];
		}
	    }
	    Rdiag[k] = -nrm;
	}
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    public boolean isFullRank () {
	for (int j = 0; j < cols; j++) {
	    if (Rdiag[j] == 0) return false;
	}
	return true;
    }

    public Matrix solve (Matrix B) {
	int cols2 = B.getCols();
	double[][] M2 = B.getArrayCopy();
	for (int k = 0; k < cols; k++) {
	    for (int j = 0; j < cols2; j++) {
		double s = 0.0; 
		for (int i = k; i < rows; i++) s += QR[i][k]*M2[i][j];
		s = -s/QR[k][k];
		for (int i = k; i < rows; i++) M2[i][j] += s*QR[i][k];
	    }
	}
	for (int k = cols-1; k >= 0; k--) {
	    for (int j = 0; j < cols2; j++) M2[k][j] /= Rdiag[k];
	    for (int i = 0; i < k; i++)
		for (int j = 0; j < cols2; j++)
		    M2[i][j] -= M2[k][j]*QR[i][k];
	}
	return (new Matrix(M2).getMatrix(0,cols-1,0,cols2-1));
    }
}


