
/*
 * File SnAPLikelihoodCore.java
 *
 * Copyright (C) 2010 Remco Bouckaert, David Bryant remco@cs.auckland.ac.nz
 *
 * This file is part of SnAP.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * SnAP is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  SnAP is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SnAP; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */
package snapper;


import java.util.Arrays;

import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.BeerLikelihoodCore;
import beast.evolution.tree.Node;

public class SnapperLikelihoodCore extends BeerLikelihoodCore {
    boolean m_bReuseCache = false;
    final double MIN_STEP = 0.01;
    final static boolean debug = false;
    
    
    // 2 x #nodes x #patterns at bottom of branch 
    ChebyshevPolynomial[][][] chebPoly;
    double[][]time;
    int N;

	double [] v1;
	double [] v2;
	double [] Q1;
	double [] Q2;
	double [] delta;

    public SnapperLikelihoodCore(Node root, Alignment data, int N) {
    	super(N);
    	this.N = N;
    	int nodeCount = root.getNodeCount();
    	int patternCount = data.getPatternCount(); 
    	chebPoly = new ChebyshevPolynomial[2][nodeCount][patternCount];
    	for (int i = 0; i < nodeCount; i++) {
    		for (int j = 0; j < patternCount; j++) {
    			chebPoly[0][i][j] = new ChebyshevPolynomial(N);
    			chebPoly[1][i][j] = new ChebyshevPolynomial(N);
    		}
    	}
    	v1 = new double[N];
    	v2 = new double[N];
    	Q1 = new double[N*N];
    	Q2 = new double[N*N];
    	
    	delta = ChebyshevPolynomial.getInterValContributions(N);
    	
    	states = new int[nodeCount][];
    }
    
    void exponentiate(double time, double [] matrix, double [] a) {
    	 MatrixExponentiator.expmvRK4(time, matrix, a, MIN_STEP);

    	 ChebyshevPolynomial c = new ChebyshevPolynomial(N);
    	 c.a = a;
    	 c.aToF();
    	 System.arraycopy(c.f, 0, a, 0, N);
    	
//		double [] a2 = new double[N+1];
//		System.arraycopy(a, 0, a2, 1, N);
//		QMatrix Q = new QMatrix(matrix);
//		double[] b = null;
//		try {
//			b = MatrixExponentiator.expmv(time, Q, a2);
//		} catch (Exception e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		System.arraycopy(b, 1, a, 0, N);

    }
    
    @Override
    public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories,
    		boolean useAmbiguities) {
    	super.initialize(nodeCount, patternCount, matrixCount, integrateCategories, useAmbiguities);
    	time = new double[nodeCount][matrixCount];
    }
    
    public void setLeafPolyFactors(int nodeIndex, int patternIndex, int r, int n) {
    	chebPoly[0][nodeIndex][patternIndex].init(r, n);
    	chebPoly[1][nodeIndex][patternIndex].init(r, n);
    }

  
    
    public void clearCache(int nNodeNrMax, int nRedsMax) {
        // m_cache = new FCacheT(nNodeNrMax, nRedsMax + 1);
    }

    
	public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3) {
        if (states[nodeIndex1] != null) {
            if (states[nodeIndex2] != null) {
                calculateStatesStatesPruning(
                        nodeIndex1, matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                        nodeIndex2, matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                        chebPoly[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
            } else {
                calculateStatesPartialsPruning(nodeIndex1, matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                		nodeIndex2,chebPoly[currentPartialsIndex[nodeIndex2]][nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                		chebPoly[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
            }
        } else {
            if (states[nodeIndex2] != null) {
                calculateStatesPartialsPruning(nodeIndex2, matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                		nodeIndex1,chebPoly[currentPartialsIndex[nodeIndex1]][nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                		chebPoly[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
            } else {
                calculatePartialsPartialsPruning(nodeIndex1, chebPoly[currentPartialsIndex[nodeIndex1]][nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                		nodeIndex2, chebPoly[currentPartialsIndex[nodeIndex2]][nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                		chebPoly[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
            }
        }
        if (debug) {
        	System.out.println(nodeIndex1 + " x " + nodeIndex2 + " => " + nodeIndex3);
        	System.out.println(Arrays.toString(chebPoly[currentPartialsIndex[nodeIndex1]][nodeIndex1][0].f));
        	System.out.println(Arrays.toString(chebPoly[currentPartialsIndex[nodeIndex2]][nodeIndex2][0].f));
        	System.out.println(Arrays.toString(chebPoly[currentPartialsIndex[nodeIndex3]][nodeIndex3][0].f));
        }
        
        if (useScaling) {
            scalePartials(nodeIndex3);
        }
    }

    
    /**
     * Calculates partial likelihoods at a node when both children have states.
     */
    protected void calculateStatesStatesPruning(int nodeIndex1, double[] matrices1,
                                                int nodeIndex2, double[] matrices2,
                                                ChebyshevPolynomial[] partials3) {
        for (int l = 0; l < nrOfMatrices; l++) {
            int w = l * matrixSize;
            System.arraycopy(matrices1, w, Q1, 0, matrixSize);
            System.arraycopy(matrices2, w, Q2, 0, matrixSize);
            for (int k = 0; k < nrOfPatterns; k++) {

                if (debug) System.out.println(nodeIndex1 + ":exp(" + Arrays.toString(chebPoly[0][nodeIndex1][0].f)+")");
                System.arraycopy(chebPoly[0][nodeIndex1][k].a, 0, v1, 0, N);
                exponentiate(time[nodeIndex1][l], Q1, v1);
                if (debug) System.out.println("="+Arrays.toString(v1));
                
                if (debug) System.out.println(nodeIndex2 + ":exp(" + Arrays.toString(chebPoly[0][nodeIndex2][0].f)+")");
                System.arraycopy(chebPoly[0][nodeIndex2][k].a, 0, v2, 0, N);
                exponentiate(time[nodeIndex2][l], Q2, v2);
                if (debug) System.out.println("="+Arrays.toString(v2));
                
                double [] f = partials3[k].f; 
                for (int i = 0; i < N; i++) {
                	f[i] = Math.max(v1[i] * v2[i], 0);
                }
            }
        }
    }

    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     */
    protected void calculateStatesPartialsPruning(int nodeIndex1, double[] matrices1,
                                                  int nodeIndex2, ChebyshevPolynomial[] partials2, double[] matrices2,
                                                  ChebyshevPolynomial[] partials3) {
        for (int l = 0; l < nrOfMatrices; l++) {
            int w = l * matrixSize;
            System.arraycopy(matrices1, w, Q1, 0, matrixSize);
            System.arraycopy(matrices2, w, Q2, 0, matrixSize);
            for (int k = 0; k < nrOfPatterns; k++) {

                System.arraycopy(chebPoly[0][nodeIndex1][k].a, 0, v1, 0, N);
                exponentiate(time[nodeIndex1][l], Q1, v1);
                
                partials2[k].fToA();
                System.arraycopy(partials2[k].a, 0, v2, 0, N);
                exponentiate(time[nodeIndex2][l], Q2, v2);
                
                double [] f = partials3[k].f; 
                for (int i = 0; i < N; i++) {
                	f[i] = Math.max(v1[i] * v2[i], 0);
                }
            }
        }

    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected void calculatePartialsPartialsPruning(int nodeIndex1, ChebyshevPolynomial[] partials1, double[] matrices1,
    		int nodeIndex2, ChebyshevPolynomial[] partials2, double[] matrices2,
    		ChebyshevPolynomial[] partials3) {
        for (int l = 0; l < nrOfMatrices; l++) {
            int w = l * matrixSize;
            System.arraycopy(matrices1, w, Q1, 0, matrixSize);
            System.arraycopy(matrices2, w, Q2, 0, matrixSize);
            for (int k = 0; k < nrOfPatterns; k++) {

            	partials1[k].fToA();
                System.arraycopy(partials1[k].a, 0, v1, 0, N);
                exponentiate(time[nodeIndex1][l], Q1, v1);
                
                partials2[k].fToA();
                System.arraycopy(partials2[k].a, 0, v2, 0, N);
                exponentiate(time[nodeIndex2][l], Q2, v2);
                
                double [] f = partials3[k].f; 
                for (int i = 0; i < N; i++) {
                	f[i] = Math.max(v1[i] * v2[i], 0);
                }
            }
        }
    }

    @Override
    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {
        calculateIntegratePartials(chebPoly[currentPartialsIndex[nodeIndex]][nodeIndex], proportions, outPartials);
    }
    
    /**
     * Integrates partials across categories.
     *
     * @param inPartials  the array of partials to be integrated
     * @param proportions the proportions of sites in each category
     * @param outPartials an array into which the partials will go
     */
	protected void calculateIntegratePartials(ChebyshevPolynomial[] inPartials, double[] proportions, double[] outPartials) {

        int u = 0;
//        int v = 0;
        double [] f;
        for (int k = 0; k < nrOfPatterns; k++) {
        	f = inPartials[k].f;
            for (int i = 0; i < nrOfStates; i++) {
                outPartials[u] = f[i] * proportions[0];
                u++;
            }
        }


        for (int l = 1; l < nrOfMatrices; l++) {
            u = 0;

            for (int k = 0; k < nrOfPatterns; k++) {
            	f = inPartials[k].f;

                for (int i = 0; i < nrOfStates; i++) {

                    outPartials[u] += f[i] * proportions[l];
                    u++;
                }
            }
        }
    }

    /**
     * Calculates pattern log likelihoods at a node.
     *
     * @param partials          the partials used to calculate the likelihoods
     * @param frequencies       an array of state frequencies
     * @param outLogLikelihoods an array into which the likelihoods will go
     */
    @Override
	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
        int v = 0;
        for (int k = 0; k < nrOfPatterns; k++) {

            double sum = 0.0;
            for (int i = 0; i < nrOfStates; i++) {

                sum += frequencies[i] * partials[v] * delta[i];
                v++;
            }
            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
        }
    }


} // class SnapperLikelihoodCore
