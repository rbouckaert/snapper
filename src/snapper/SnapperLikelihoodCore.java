
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


import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.BeerLikelihoodCore;
import beast.evolution.tree.Node;

public class SnapperLikelihoodCore extends BeerLikelihoodCore {
    boolean m_bReuseCache = false;
    final double MIN_STEP = 0.001;
    
    
    // 2 x #nodes x #patterns at bottom of branch 
    ChebyshevPolynomial[][][] chebPoly;
    double[][]time;
    MatrixExponentiator exponentiator;
    int N;

	double [] v1;
	double [] v2;
	double [] Q1;
	double [] Q2;

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
    	exponentiator = new MatrixExponentiator();
    	v1 = new double[N];
    	v2 = new double[N];
    	Q1 = new double[N*N];
    	Q2 = new double[N*N];
    }
    
    void exponentiate(double time, double [] matrix, double [] a) {
    	exponentiator.expmvRK4(time, matrix, a, MIN_STEP);
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

                System.arraycopy(chebPoly[0][nodeIndex1][k].a, 0, v1, 0, N);
                exponentiate(time[nodeIndex1][l], Q1, v1);
                
                System.arraycopy(chebPoly[0][nodeIndex2][k].a, 0, v2, 0, N);
                exponentiate(time[nodeIndex2][l], Q2, v2);
                
                double [] f = partials3[k].f; 
                for (int i = 0; i < N; i++) {
                	f[i] = v1[i] * v2[i];
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
                	f[i] = v1[i] * v2[i];
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
                	f[i] = v1[i] * v2[i];
                }
            }
        }
    }

    /**
     * Integrates partials across categories.
     *
     * @param inPartials  the array of partials to be integrated
     * @param proportions the proportions of sites in each category
     * @param outPartials an array into which the partials will go
     */
    @Override
	protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials) {

        int u = 0;
        int v = 0;
        for (int k = 0; k < nrOfPatterns; k++) {

            for (int i = 0; i < nrOfStates; i++) {

                outPartials[u] = inPartials[v] * proportions[0];
                u++;
                v++;
            }
        }


        for (int l = 1; l < nrOfMatrices; l++) {
            u = 0;

            for (int k = 0; k < nrOfPatterns; k++) {

                for (int i = 0; i < nrOfStates; i++) {

                    outPartials[u] += inPartials[v] * proportions[l];
                    u++;
                    v++;
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

                sum += frequencies[i] * partials[v];
                v++;
            }
            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
        }
    }


} // class SnAPLikelihoodCore
