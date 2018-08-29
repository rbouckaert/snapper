
/*
 * File SnAPTreeLikelihood.java
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





import java.io.PrintStream;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.distribution.BetaDistribution;
import org.apache.commons.math.distribution.BetaDistributionImpl;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import snap.Data;
import snap.NodeData;


@Description("Implements a tree Likelihood Function for Single Site Sorted-sequences on a tree.") 

public class SnapperTreeLikelihood extends TreeLikelihood {
//	public Input<Data> m_pData = new Input<Data>("data", "set of alignments");
//	public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations");

	public Input<Integer> NInput = new Input<Integer>("N", "dimension of Chebyshef functions used, should be power of 2 plus 1. Higher values are more accurate but slower", 33);
	public Input<Boolean> m_bInitFromTree = new Input<Boolean>("initFromTree", "whether to initialize coalescenceRate from starting tree values (if true), or vice versa (if false)", false);
	public Input<String> m_pPattern = new Input<String>("pattern", "pattern of metadata element associated with this parameter in the tree");


	public Input<Boolean> m_usenNonPolymorphic = new Input<Boolean>("non-polymorphic", 
			"Check box only if constant sites have been left in the data and are to be included in the likelihood calculation. " +
			"Leave unchecked if all but the variable sites have been removed.",
			//"Whether to use non-polymorphic data in the sequences. " +
			//"If true, constant-sites in the data will be used as part of the likelihood calculation. " +
			//"If false (the default) constant sites will be removed from the sequence data and a normalization factor is " +
			//"calculated for the likelihood.", 
			true);
	public Input<IntegerParameter> ascSiteCountInput = new Input<IntegerParameter>("ascSiteCount", "Counts for number of ascertained sites");
	public Input<IntegerParameter> totalSiteCountInput = new Input<IntegerParameter>("totalSiteCount","Number of sites including those removed by ascertainment");
	
	public Input<Boolean> mutationOnlyAtRoot = new Input<Boolean>("mutationOnlyAtRoot", "Emulate the likelihood calculation of RoyChoudhury et al (2008) which assumes that mutations occur only in the ancestral (root) population", false);
	public Input<Boolean> hasDominantMarkers = new Input<Boolean>("dominant", "indicate that alleles are dominant (default false)", false);
	public Input<Boolean> showPatternLikelihoodsAndQuit = new Input<Boolean>("showPatternLikelihoodsAndQuit", "print out likelihoods for all patterns for the starting state, then quit", false);
	public Input<Boolean> useLogLikelihoodCorrection = new Input<Boolean>("useLogLikelihoodCorrection", "use correction of log likelihood for the purpose of calculating " +
			"Bayes factors for different species assignments. There is (almost) no computational cost involved for the MCMC chain, but the log likelihood " +
			"might be reported as positive number with this correction since the likelihood is not a proper likelihood any more.", true);
	
	public SnapperTreeLikelihood() throws Exception {
		// suppress some validation rules
		siteModelInput.setRule(Validate.OPTIONAL);
	}
	
	/** shadow variable of m_pData input */
	Data m_data2;
	/** SampleSizes = #lineages per taxon **/
	int [] m_nSampleSizes;
	/** likelihood core, doing the actual hard work of calculating the likelihood **/
	SnapperLikelihoodCore m_core;
	
	/** some variable for shadowing inputs **/
	boolean m_bUsenNonPolymorphic;
	boolean m_bMutationOnlyAtRoot;
	boolean m_bHasDominantMarkers;
//	double [] fSiteProbs;
//	double [] fStoredSiteProbs;
	
	double m_fP0 = 0.0, m_fP1 = 0.0;
	double m_fStoredP0 = 0.0, m_fStoredP1 = 0.0;
	double ascLogP = Double.NaN, storedAscLogP = Double.NaN;
	
	int N;
	
	SnapSubstitutionModel m_substitutionmodel;
	
	// Correction so that the returned value is a likelihood instead
	// of a sufficient statistic for the likelihood
	double m_fLogLikelihoodCorrection = 0;
	
	// Sampled parameter equal to the number of sites which have been removed from the data during ascertainment
	IntegerParameter ascSiteCount;
	QMatrix Q;

	@Override
	public void initAndValidate() {
		N =  NInput.get();
		Q = new QMatrix(N);
		// check N = 2^k + 1 for some k
		int k = N - 1;
		while (k > 1) {
			if (k % 2 == 1) {
				throw new IllegalArgumentException("N should be a power of 2 plus 1");
			}
			k = k / 2;
		}
		
		ascSiteCount = ascSiteCountInput.get();
		// check that alignment has same taxa as tree
    	if (!(dataInput.get() instanceof Data)) {
    		throw new IllegalArgumentException("The data input should be a snap.Data object");
    	}
    	if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
    		throw new IllegalArgumentException("The number of leaves in the tree does not match the number of sequences");
    	}

    	m_bUsenNonPolymorphic = m_usenNonPolymorphic.get();
    	m_bMutationOnlyAtRoot = mutationOnlyAtRoot.get();
    	m_bHasDominantMarkers = hasDominantMarkers.get();
    	
    	m_siteModel = (SiteModel.Base) siteModelInput.get();
    	
    	TreeInterface tree = treeInput.get();
    	m_substitutionmodel = ((SnapSubstitutionModel)m_siteModel.substModelInput.get());
    	Input<RealParameter> thetaInput = m_substitutionmodel.thetaInput;
		
		Double [] values = new Double[tree.getNodeCount()];
		String sTheta = "";
		if (m_bInitFromTree.get() == true) {
			tree.getMetaData(tree.getRoot(), values, m_pPattern.get());
			for (Double d : values) {
				sTheta += d + " ";
			}
		} else {
	    	List<Double> sValues = thetaInput.get().valuesInput.get();
	        for (int i = 0; i < values.length; i++) {
	            values[i] = new Double(sValues.get(i % sValues.size()));
	            sTheta += values[i] + " ";
	        }
			tree.setMetaData(tree.getRoot(), values, m_pPattern.get());
		}
		RealParameter pTheta = thetaInput.get();
		RealParameter theta = new RealParameter();
		theta.initByName("value", sTheta, "upper", pTheta.getUpper(), "lower", pTheta.getLower(), "dimension", values.length);
		theta.setID(pTheta.getID());
		thetaInput.get().assignFrom(theta);
	
    	
    	
    	
    	m_data2 = (Data) dataInput.get();
    	if ( BeastMCMC.m_nThreads == 1) {
    		// single threaded likelihood core
    		m_core = new SnapperLikelihoodCore(treeInput.get().getRoot(), dataInput.get(), N);
    	} else {
    		// multi-threaded likelihood core
    		throw new RuntimeException("Threading not implemented yet -- use single tread");
    		//m_core = new SnapperLikelihoodCoreT(treeInput.get().getRoot(), dataInput.get());
    	}
    	Integer [] nSampleSizes = m_data2.getStateCounts().toArray(new Integer[0]);
    	m_nSampleSizes = new int[nSampleSizes.length];
    	for (int i = 0; i < nSampleSizes.length; i++) {
    		m_nSampleSizes[i] = nSampleSizes[i];
    	}
    	if (!(treeInput.get().getRoot() instanceof NodeData)) {
    		throw new IllegalArgumentException("Tree has nodes of the wrong type. NodeData expected, but found " + 
    				treeInput.get().getRoot().getClass().getName());
    	}

		int numPatterns = m_data2.getPatternCount();
		
		
		
		// calculate Likelihood Correction. 
		// When the assignment of individuals to populations/species is fixed, the allele counts in each population are sufficient 
		// statistics for the species tree parameters. However when testing species assignments this is no longer the case.
		// To address this we multiply the likelihood computed from allele counts by the probability of observing
		// the given sequences given those allele counts (and the species assignments).
		m_fLogLikelihoodCorrection = 0;
		if (useLogLikelihoodCorrection.get()) {
			// RRB: note that increasing the number of constant sites
			// does not change the m_fLogLikelihoodCorrection since the
			// contribution of constant sites is zero. This means,
			// m_fLogLikelihoodCorrection does not need to be recalculated
			// when ascSiteCount changes.
			// DJB: This is true, but only until we start looking at non-constant sites being ascertained.
	    	for (int i = 0; i < numPatterns; i++) {
	            int [] thisSite = m_data2.getPattern(i);  //count of red alleles for this site
	            int [] lineageCounts = m_data2.getPatternLineagCounts(i); //count of total lineages for this site
	            for (int j = 0; j < thisSite.length; j++) {
	            	m_fLogLikelihoodCorrection -= logBinom(thisSite[j], lineageCounts[j]) * m_data2.getPatternWeight(i);
	            }
	    	}
    	}
		System.err.println("Log Likelihood Correction = " + m_fLogLikelihoodCorrection);
		
		branchRateModel = branchRateModelInput.get();
		if (branchRateModel == null) {
			branchRateModel = new StrictClockModel();
		}
		if (!(branchRateModel instanceof StrictClockModel)) {
			//We assume that the mutation rate (but not theta) is constant for the species tree.
			throw new IllegalArgumentException("Only strict clock model allowed for branchRateModel, not " + branchRateModel.getClass().getName());
		}

		
		
		
    	
    	int nodeCount = tree.getNodeCount();
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        patternLogLikelihoods = new double[numPatterns];
        m_fRootPartials = new double[numPatterns * N];
        
        m_core.initialize(
                nodeCount,
                numPatterns,
                m_siteModel.getCategoryCount(),
                true, m_useAmbiguities.get()
        );

		// set up Chebyshef functions for leaf observations
		// not particularly efficient, but we only do this once
		double [] f = new double[N];
    	for (int i = 0; i < numPatterns; i++) {
            int [] thisSite = m_data2.getPattern(i);
            int [] lineageCounts = m_data2.getPatternLineagCounts(i);
            
            for (int j = 0; j < thisSite.length; j++) {
            	int r = thisSite[j];
            	int n = lineageCounts[j];
            	m_core.setLeafPolyFactors(j, i, r, n);
            }
    	}
    }

    private double logBinom(int k, int n) {
    	double f = 0;
    	for (int i = k + 1; i <= n; i++) {
    		f += Math.log(i) - Math.log(n - i + 1);
    	}
		return f;
	}

	/**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLogP() {
    	try {
    		// get current tree
	    	NodeData root = (NodeData) treeInput.get().getRoot();
	    	//Double [] theta = m_substitutionmodel.thetaInput.get().getValues();
	    	// assing gamma values to tree
//	    	if (m_pGamma.get().somethingIsDirty()) {
//	    		// sync gammas in parameter with gammas in tree, if necessary
//	    		m_pGamma.get().prepare();
//	    	}
	    	
	    	//double u = m_substitutionmodel.m_pU.get().getValue();
	    	//double v  = m_substitutionmodel.m_pV.get().getValue();
			boolean useCache = false;
			//boolean useCache = false;
			boolean dprint = showPatternLikelihoodsAndQuit.get();
			if (dprint) {
				System.out.println("Log Likelihood Correction = " + m_fLogLikelihoodCorrection);
			}
			
			
			double [] fCategoryRates = m_siteModel.getCategoryRates(null);
			if (branchRateModel != null) {
				double branchRate = branchRateModel.getRateForBranch(null);
				for (int i = 0; i < fCategoryRates.length; i++) {
					fCategoryRates[i] *= branchRate;
				}
			}
			
			traverse(root);
			
			// amalgamate site probabilities over patterns
			int numPatterns = m_data2.getPatternCount();
			// claculate log prob
			logP = 0;
			for(int id = 0; id < numPatterns - (m_bUsenNonPolymorphic ? 0 : 2); id++) {
				double freq = m_data2.getPatternWeight(id);
				double siteL = patternLogLikelihoods[id];
				if (Double.isInfinite(siteL)) {
					logP = -10e100;
					break;
				}
				logP += (double)freq * siteL;
			}
			// correction for constant sites. If we are sampling the numbers of constant sites 
			// (stored in ascSiteCount) then we include these probabilities. Otherwise we 
			// assume that we want conditional likelihood, in which case we divide 
			// by the probability that a site is not ascertained (or more correctly,
			// subtract the log probability.
			if (!m_bUsenNonPolymorphic) {
				m_fP0 =  patternLogLikelihoods[numPatterns - 2];
				m_fP1 =  patternLogLikelihoods[numPatterns - 1];
				if (ascSiteCount != null) {   
					ascLogP = (double)ascSiteCount.getValue(0) * Math.log(m_fP0) +
							  (double)ascSiteCount.getValue(1) * Math.log(m_fP1);
					logP += ascLogP;
				} else {
					logP -= (double) m_data2.getSiteCount() * Math.log(1.0 - m_fP0 - m_fP1);
				}
			}				
			
			if (useLogLikelihoodCorrection.get()) {
				logP += m_fLogLikelihoodCorrection;
			}
			
			
//			logP = m_core.computeLogLikelihood(root, u , v, 
//	    			m_nSampleSizes, 
//	    			m_data2,
//	    			coalescenceRate,
//	    			fCategoryRates, fCategoryProportions,
//	    			useCache,
//	    			m_bUsenNonPolymorphic,
//	    			dprint /*= false*/);
			return logP;
    	} catch (Exception e) {
			e.printStackTrace();
			return 0;
		}
    } // calculateLogLikelihood

    
    @Override
    protected boolean requiresRecalculation() {
    	boolean isDirty = super.requiresRecalculation();
    	if (ascSiteCount != null && ascSiteCount.somethingIsDirty()) {
    		isDirty = true;
    	}
    	return isDirty;
    }
    
    /** CalculationNode methods **/
	@Override
	public void store() {
        storedLogP = logP;
    	m_core.m_bReuseCache = true;
    	super.store();
    	m_fStoredP0 = m_fP0;
    	m_fStoredP1 = m_fP1;
    	storedAscLogP = ascLogP; 
    }

	@Override
    public void restore() {
        logP = storedLogP;
    	m_core.m_bReuseCache = false;
    	super.restore();
    	m_fP0 = m_fStoredP0;
    	m_fP1 = m_fStoredP1;
    	ascLogP = storedAscLogP; 
    }

	
	@Override public List<String> getArguments() {return null;}
	@Override public List<String> getConditions() {return null;}
	@Override public void sample(State state, Random random) {};
	
	@Override
	public void init(PrintStream out) {
		super.init(out);
		if (!m_bUsenNonPolymorphic) {
			out.append("P0\t");
			out.append("P1\t");
		}
	}
	
	@Override
	public void log(long nSample, PrintStream out) {
		super.log(nSample, out);
		if (!m_bUsenNonPolymorphic) {
			out.append(m_fP0 + "\t");
			out.append(m_fP1 + "\t");
		}
	}

	public double getProbVariableSites() {
		if (!m_bUsenNonPolymorphic) {
			return 1.0 - m_fP0 - m_fP1;
		} else {
			return 1.0;
		}
	}

	
	// return contribution of ascertained sites to log likelihood
	public double getAscSitesLogP() {
		if (ascSiteCount != null) {
			return ascLogP;
		} else {
			return 0.0;
		}
	}
	
	
    /* Assumes there IS a branch rate model as opposed to traverse() */
    int traverse(final Node node) {

        int update = (node.isDirty() | hasDirt);
        update = Tree.IS_FILTHY;

        final int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex] ||
        		m_substitutionmodel.thetaInput.get().getStoredValue(nodeIndex) != m_substitutionmodel.thetaInput.get().getValue(nodeIndex))) {
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            m_core.setNodeMatrixForUpdate(nodeIndex);
    		Double[] theta = m_substitutionmodel.thetaInput.get().getValues();
    		double u = m_substitutionmodel.m_pU.get().getValue();
    		double v = m_substitutionmodel.m_pV.get().getValue();
			double[] fCategoryRates = m_siteModel.getCategoryRates(null);
			double [] time = m_core.time[nodeIndex];
			
    		for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
    			double scaledTheta = theta[nodeIndex] / fCategoryRates[i];
            	final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
            	time[i] = jointBranchRate * branchTime;// * scaledTheta / 2;
            	Q.setQ(u, v, scaledTheta);
                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                m_core.setNodeMatrix(nodeIndex, i, Q.Q);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                m_core.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                	m_core.setNodeStatesForUpdate(nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {
                	m_core.calculatePartials(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods
                    final double[] frequencies = getFrequencies();

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    m_core.integratePartials(node.getNr(), proportions, m_fRootPartials);

                    m_core.calculateLogLikelihoods(m_fRootPartials, frequencies, patternLogLikelihoods);
                }

            }
        }
        return update;
    } // traverseWithBRM

    ChebyshevPolynomial freqs = null;
	private double[] getFrequencies() {
		if (freqs == null) {
			freqs = new ChebyshevPolynomial(N);
			freqs.a[0] = 1;
			freqs.aToF();
			
//			BetaDistribution beta = new BetaDistributionImpl(0.006, 0.006);
//			freqs.f[0] = 10;
//			freqs.f[N-1] = 10;
//			for (int i = 1; i < N-1; i++) {
//	    		double x = 0.5 - Math.cos(-i/(N-1.0)*Math.PI) / 2.0;
//				freqs.f[i] = beta.density(x);
//			}
		}
		
		return freqs.f;
	}

} // class SnapperTreeLikelihood