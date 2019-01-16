
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

import org.apache.commons.math.distribution.BetaDistribution;
import org.apache.commons.math.distribution.BetaDistributionImpl;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.special.Beta;

import beast.app.BeastMCMC;
import beast.app.beauti.Beauti;
import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
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
	

	public Input<Boolean> useBetaRootPriorInput = new Input<Boolean>("useBetaRootPrior", "instead of using a uniform prior for allele frequencies at the root, "
			+ "use a beta root prior", false);

	final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads", "maximum number of threads to use, if less than 1 the number of threads in BeastMCMC is used (default -1)", -1);
    final public Input<String> proportionsInput = new Input<>("proportions", "specifies proportions of patterns used per thread as space "
    		+ "delimited string. This is useful when using a mixture of BEAGLE devices that run at different speeds, e.g GPU and CPU. "
    		+ "The string is duplicated if there are more threads than proportions specified. For example, "
    		+ "'1 2' as well as '33 66' with 2 threads specifies that the first thread gets a third of the patterns and the second "
    		+ "two thirds. With 3 threads, it is interpreted as '1 2 1' = 25%, 50%, 25% and with 7 threads it is "
    		+ "'1 2 1 2 1 2 1' = 10% 20% 10% 20% 10% 20% 10%. If not specified, all threads get the same proportion of patterns.");
	
    /** private list of likelihoods, to notify framework of TreeLikelihoods being created in initAndValidate() **/
    final private Input<List<TreeLikelihood>> likelihoodsInput = new Input<>("*","",new ArrayList<>());

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

	
    /** number of threads to use, changes when threading causes problems **/
    private int threadCount;
    private ExecutorService pool = null;
    private final List<Callable<Double>> likelihoodCallers = new ArrayList<Callable<Double>>();
	// specified a set ranges of patterns assigned to each thread
	// first patternPoints contains 0, then one point for each thread
    private int [] patternPoints;
	
	SnapperTreeLikelihood [] treelikelihood;
	double [] logPByThread;
	
	
	@Override
	public void initAndValidate() {
		threadCount = BeastMCMC.m_nThreads;

		if (maxNrOfThreadsInput.get() > 0) {
			threadCount = Math.min(maxNrOfThreadsInput.get(), BeastMCMC.m_nThreads);
		}
        String instanceCount = System.getProperty("beast.instance.count");
        if (instanceCount != null && instanceCount.length() > 0) {
        	threadCount = Integer.parseInt(instanceCount);
        }

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
		

		branchRateModel = branchRateModelInput.get();
		if (branchRateModel == null) {
			branchRateModel = new StrictClockModel();
		}
		if (!(branchRateModel instanceof StrictClockModel)) {
			//We assume that the mutation rate (but not theta) is constant for the species tree.
			throw new IllegalArgumentException("Only strict clock model allowed for branchRateModel, not " + branchRateModel.getClass().getName());
		}


    	if (threadCount <= 1) {
    		m_core = new SnapperLikelihoodCore(treeInput.get().getRoot(), dataInput.get(), N);
    	} else {
    		treelikelihood = new SnapperTreeLikelihood[threadCount];
    		logPByThread = new double[threadCount];
    		
    		// multi-threaded likelihood core
        	pool = Executors.newFixedThreadPool(threadCount);
    		
        	calcPatternPoints(dataInput.get().getSiteCount());
        	for (int i = 0; i < threadCount; i++) {
        		Alignment rawdata = dataInput.get();
        		String filterSpec = (patternPoints[i] +1) + "-" + (patternPoints[i + 1]);
        		if (rawdata.isAscertained) {
        			filterSpec += rawdata.excludefromInput.get() + "-" + rawdata.excludetoInput.get() + "," + filterSpec;
        		}
        		SnapperTreeLikelihood likelihood = null;
				try {
					likelihood = new SnapperTreeLikelihood();
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
        		likelihood.setID(getID() + i);
        		likelihood.getOutputs().add(this);
        		likelihoodsInput.get().add(likelihood);

        		snap.FilteredAlignment filter = new snap.FilteredAlignment();
        		if (i == 0 && dataInput.get() instanceof FilteredAlignment && ((FilteredAlignment)dataInput.get()).constantSiteWeightsInput.get() != null) {
        			filter.initByName("data", dataInput.get()/*, "userDataType", m_data.get().getDataType()*/, 
        							"filter", filterSpec, 
        							"constantSiteWeights", ((FilteredAlignment)dataInput.get()).constantSiteWeightsInput.get()
        							);
        		} else {
        			filter.initByName("data", dataInput.get()/*, "userDataType", m_data.get().getDataType()*/, 
        							"filter", filterSpec,
        							"userDataType", dataInput.get().getDataType()
        							);
        		}
        		likelihood.initByName("data", filter, 
        				"tree", treeInput.get(), 
        				"siteModel", duplicate((BEASTInterface) siteModelInput.get(), i), 
        				"branchRateModel", duplicate(branchRateModelInput.get(), i), 
        				"initFromTree", m_bInitFromTree.get(),
                        "pattern" , m_pPattern.get() + "",
                        "threads", 1,
                        "useBetaRootPrior", useBetaRootPriorInput.get()                        
        				);
        	    treelikelihood[i] = likelihood;        		
        		likelihoodCallers.add(new TreeLikelihoodCaller(i, likelihood));
        	}
        	return;
    	}
    	

		
		
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

    
    class TreeLikelihoodCaller implements Callable<Double> {
        private final int threadNr;
        private final SnapperTreeLikelihood likelihood;
        
        public TreeLikelihoodCaller(int threadNr, SnapperTreeLikelihood likelihood) {
        	this.threadNr = threadNr;
        	this.likelihood = likelihood;
        }

        public Double call() throws Exception {
  		  	try {
	            logPByThread[threadNr] = likelihood.calculateLogP();
  		  	} catch (Exception e) {
  		  		System.err.println("Something went wrong in thread " + threadNr);
				e.printStackTrace();
				System.exit(0);
			}
            return 0.0;
        }

    }
    
    
    /** create new instance of src object, connecting all inputs from src object
     * Note if input is a SubstModel, it is duplicated as well.
     * @param src object to be copied
     * @param i index used to extend ID with.
     * @return copy of src object
     */
    private Object duplicate(BEASTInterface src, int i) {
    	if (src == null) { 
    		return null;
    	}
    	BEASTInterface copy;
		try {
			copy = src.getClass().newInstance();
        	copy.setID(src.getID() + "_" + i);
		} catch (InstantiationException | IllegalAccessException e) {
			e.printStackTrace();
			throw new RuntimeException("Programmer error: every object in the model should have a default constructor that is publicly accessible: " + src.getClass().getName());
		}
        for (Input<?> input : src.listInputs()) {
            if (input.get() != null) {
                if (input.get() instanceof List) {
                    // handle lists
                	//((List)copy.getInput(input.getName())).clear();
                    for (Object o : (List<?>) input.get()) {
                        if (o instanceof BEASTInterface) {
                        	// make sure it is not already in the list
                            copy.setInputValue(input.getName(), o);
                        }
                    }
                } else if (input.get() instanceof SubstitutionModel) {
                	// duplicate subst models
                	BEASTInterface substModel = (BEASTInterface) duplicate((BEASTInterface) input.get(), i);
            		copy.setInputValue(input.getName(), substModel);
            	} else {
                    // it is some other value
            		copy.setInputValue(input.getName(), input.get());
            	}
            }
        }
        copy.initAndValidate();
		return copy;
	}
    
	private void calcPatternPoints(int nPatterns) {
		patternPoints = new int[threadCount + 1];
		if (proportionsInput.get() == null) {
			int range = nPatterns / threadCount;
			for (int i = 0; i < threadCount - 1; i++) {
				patternPoints[i+1] = range * (i+1);
			}
			patternPoints[threadCount] = nPatterns;
		} else {
			String [] strs = proportionsInput.get().split("\\s+");
			double [] proportions = new double[threadCount];
			for (int i = 0; i < threadCount; i++) {
				proportions[i] = Double.parseDouble(strs[i % strs.length]);
			}
			// normalise
			double sum = 0;
			for (double d : proportions) {
				sum += d;
			}
			for (int i = 0; i < threadCount; i++) {
				proportions[i] /= sum;
			}
			// cummulative 
			for (int i = 1; i < threadCount; i++) {
				proportions[i] += proportions[i- 1];
			}
			// calc ranges
			for (int i = 0; i < threadCount; i++) {
				patternPoints[i+1] = (int) (proportions[i] * nPatterns + 0.5);
			}
		}
    }
	/**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLogP() {
		if (threadCount > 1) {
			try {
	            pool.invokeAll(likelihoodCallers);
			} catch (RejectedExecutionException | InterruptedException e) {
				e.printStackTrace();
				System.exit(0);
			}
			logP = 0;
	    	for (double f : logPByThread) {
	    		logP += f;
	    	}
	    	return logP;
		}

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
				// System.out.println(Arrays.toString(m_data2.getPattern(id)) + " " + id + " " + (freq * siteL));
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
    	if (threadCount > 1) {
    		boolean requiresRecalculation = false;
    		for (SnapperTreeLikelihood b : treelikelihood) {
    			requiresRecalculation |= b.requiresRecalculation();
    		}
    		return requiresRecalculation;    		
    	}
    	
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
    	if (threadCount == 1) {
	    	m_core.m_bReuseCache = true;
	    	m_fStoredP0 = m_fP0;
	    	m_fStoredP1 = m_fP1;
	    	storedAscLogP = ascLogP;
	    	super.store();
    	}
    }

	@Override
    public void restore() {
		logP = storedLogP;
    	if (threadCount == 1) {
    		m_core.m_bReuseCache = false;
    		m_fP0 = m_fStoredP0;
    		m_fP1 = m_fStoredP1;
    		ascLogP = storedAscLogP;
        	super.restore();
    	}
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
        		m_substitutionmodel.thetaInput.get().isDirty(nodeIndex))) {
        		//m_substitutionmodel.thetaInput.get().getStoredValue(nodeIndex) != m_substitutionmodel.thetaInput.get().getValue(nodeIndex)
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
                    final double[] frequencies = Tn_ito_xn_lookup_fast(); // getFrequencies();

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    m_core.integratePartials(node.getNr(), proportions, m_fRootPartials);

                    m_core.calculateLogLikelihoods(m_fRootPartials, frequencies, patternLogLikelihoods);
                }

            }
        }
        return update;
    } // traverse

    ChebyshevPolynomial freqs = null;
    
	// returns root allele frequencies
	private double[] getFrequencies() {
		if (freqs == null) {
			freqs = new ChebyshevPolynomial(N);
			if (!useBetaRootPriorInput.get()) {
				// use uniform prior for root allele frequencies
				freqs.a[0] = 1;
				freqs.aToF();
			} else {
	    		double u = m_substitutionmodel.m_pU.get().getValue();
	    		double v = m_substitutionmodel.m_pV.get().getValue();
	    		Double [] thetas = m_substitutionmodel.thetaInput.get().getValues();
	    		double theta = thetas[thetas.length - 1]; // theta for the root	    				

	    		BetaDistribution distr;
	    		if (Math.abs(u -v) < 1e-10) { 
	    			// u == v
	    			distr = new BetaDistributionImpl(theta/2.0, theta/2.0);
	    		} else {
	    			// u != v
		    		double alpha = u/v;
					double beta1 = (-alpha * theta - theta)/(2.0 * (alpha * alpha * theta + 2 * alpha * theta - 2 * alpha + theta));
					double beta2 =  beta1 * alpha;
					distr = new BetaDistributionImpl(beta1, beta2);
				}
				freqs.f[0] = distr.density(0.0001);
				freqs.f[N-1] = distr.density(0.9999);
				for (int i = 1; i < N-1; i++) {
		    		double x = 0.5 - Math.cos(-i/(N-1.0)*Math.PI) / 2.0;
					freqs.f[i] = distr.density(x);
				}
			}			
		}
		
		return freqs.f;
	}
	
	
	private static double []  D = new double [] {	2.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.000000000000000000e+00,
			8.000000000000000000e+00,-8.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			3.200000000000000000e+01,-4.800000000000000000e+01,1.800000000000000000e+01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.000000000000000000e+00,
			1.280000000000000000e+02,-2.560000000000000000e+02,1.600000000000000000e+02,-3.200000000000000000e+01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			5.120000000000000000e+02,-1.280000000000000000e+03,1.120000000000000000e+03,-4.000000000000000000e+02,5.000000000000000000e+01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.000000000000000000e+00,
			2.048000000000000000e+03,-6.144000000000000000e+03,6.912000000000000000e+03,-3.584000000000000000e+03,8.400000000000000000e+02,-7.200000000000000000e+01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			8.192000000000000000e+03,-2.867200000000000000e+04,3.942400000000000000e+04,-2.688000000000000000e+04,9.408000000000000000e+03,-1.568000000000000000e+03,9.800000000000000000e+01,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.000000000000000000e+00,
			3.276800000000000000e+04,-1.310720000000000000e+05,2.129920000000000000e+05,-1.802240000000000000e+05,8.448000000000000000e+04,-2.150400000000000000e+04,2.688000000000000000e+03,-1.280000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			1.310720000000000000e+05,-5.898240000000000000e+05,1.105920000000000000e+06,-1.118208000000000000e+06,6.589440000000000000e+05,-2.280960000000000000e+05,4.435200000000000000e+04,-4.320000000000000000e+03,1.620000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.000000000000000000e+00,
			5.242880000000000000e+05,-2.621440000000000000e+06,5.570560000000000000e+06,-6.553600000000000000e+06,4.659200000000000000e+06,-2.050048000000000000e+06,5.491200000000000000e+05,-8.448000000000000000e+04,6.600000000000000000e+03,-2.000000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			2.097152000000000000e+06,-1.153433600000000000e+07,2.739404800000000000e+07,-3.676569600000000000e+07,3.063808000000000000e+07,-1.640038400000000000e+07,5.637632000000000000e+06,-1.208064000000000000e+06,1.510080000000000000e+05,-9.680000000000000000e+03,2.420000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.000000000000000000e+00,
			8.388608000000000000e+06,-5.033164800000000000e+07,1.321205760000000000e+08,-1.992294400000000000e+08,1.905131520000000000e+08,-1.203240960000000000e+08,5.069209600000000000e+07,-1.405747200000000000e+07,2.471040000000000000e+06,-2.562560000000000000e+05,1.372800000000000000e+04,-2.880000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			3.355443200000000000e+07,-2.181038080000000000e+08,6.270484480000000000e+08,-1.049624576000000000e+09,1.133117440000000000e+09,-8.255569920000000000e+08,4.127784960000000000e+08,-1.412136960000000000e+08,3.236147200000000000e+07,-4.759040000000000000e+06,4.164160000000000000e+05,-1.892800000000000000e+04,3.380000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.000000000000000000e+00,
			1.342177280000000000e+08,-9.395240960000000000e+08,2.936012800000000000e+09,-5.402263552000000000e+09,6.499598336000000000e+09,-5.369233408000000000e+09,3.111714816000000000e+09,-1.270087680000000000e+09,3.611811840000000000e+08,-6.970163200000000000e+07,8.712704000000000000e+06,-6.522880000000000000e+05,2.548000000000000000e+04,-3.920000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			5.368709120000000000e+08,-4.026531840000000000e+09,1.358954496000000000e+10,-2.726297600000000000e+10,3.617587199900000000e+10,-3.342650572700000000e+10,2.205220863900000000e+10,-1.047822335900000000e+10,3.572121598000000000e+09,-8.599551980000000000e+08,1.418926060000000000e+08,-1.527551800000000000e+07,9.900780000000000000e+05,-3.359800000000000000e+04,4.480000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			2.147483648000000000e+09,-1.717986918400000000e+10,6.227702579200000000e+10,-1.352914698240000000e+11,1.962934272000000000e+11,-2.006555033600000000e+11,1.485622476800000000e+11,-8.064807731200000000e+10,3.213321830400000000e+10,-9.313976320000000000e+09,1.926299648000000000e+09,-2.751856640000000000e+08,2.579865600000000000e+07,-1.462272000000000000e+06,4.352000000000000000e+04,-5.120000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			8.589934592000000000e+09,-7.301444403200000000e+10,2.829309706240000000e+11,-6.616933990400000000e+11,1.042167103488000000e+12,-1.167945891840000000e+12,9.593841254410000000e+11,-5.862902988810000000e+11,2.677768192010000000e+11,-9.104411852900000000e+10,2.276102963300000000e+10,-4.093386753000000000e+09,5.116733450000000000e+08,-4.217088100000000000e+07,2.108545000000000000e+06,-5.548900000000000000e+04,5.790000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-2.000000000000000000e+00,
			3.435973836800000000e+10,-3.092376453120000000e+11,1.275605286912000000e+12,-3.195455668224000000e+12,5.429778186240000000e+12,-6.620826304512000000e+12,5.977134858240000000e+12,-4.063273943040000000e+12,2.095125626880000000e+12,-8.190820352000000000e+11,2.409991372800000000e+11,-5.258162995200000000e+10,8.307167232000000000e+09,-9.168445440000000000e+08,6.697728000000000000e+07,-2.976768000000000000e+06,6.976800000000000000e+04,-6.480000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			1.374389534720000000e+11,-1.305670057984000000e+12,5.712306503680000000e+12,-1.526001880268800000e+13,2.782709311078400000e+13,-3.668116819148800000e+13,3.610802493849600000e+13,-2.703941959680000000e+13,1.554766626816000000e+13,-6.880289095680000000e+12,2.334383800320000000e+12,-6.012806758400000000e+11,1.156308992000000000e+11,-1.618832588800000000e+10,1.589924864000000000e+09,-1.036907520000000000e+08,4.124064000000000000e+06,-8.664000000000000000e+04,7.220000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.000000000000000000e+00,
			5.497558138880000000e+11,-5.497558138880000000e+12,2.542620639232000000e+13,-7.215545057280000000e+13,1.405528047616000000e+14,-1.991834033192960000e+14,2.123646579507200000e+14,-1.737529019596800000e+14,1.102923694080000000e+14,-5.455321497600000000e+13,2.100298776576000000e+13,-6.254808268800000000e+12,1.424085811200000000e+12,-2.434334720000000000e+11,3.042918400000000000e+10,-2.677768192000000000e+09,1.569004800000000000e+08,-5.617920000000000000e+06,1.064000000000000000e+05,-8.000000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			2.199023255552000000e+12,-2.308974418329600000e+13,1.125625028935680000e+14,-3.381685450178560000e+14,7.008098136883200000e+14,-1.062579203997696000e+15,1.219998345330688000e+15,-1.083059755548672000e+15,7.525672566128640000e+14,-4.117581791232000000e+14,1.775707147468800000e+14,-6.014491951104000000e+13,1.587157598208000000e+13,-3.220624834560000000e+12,4.929527808000000000e+11,-5.538111488000000000e+10,4.393213440000000000e+09,-2.325818880000000000e+08,7.537376000000000000e+06,-1.293600000000000000e+05,8.820000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.000000000000000000e+00,
			8.796093022208000000e+12,-9.675702324428800000e+13,4.958797441269760000e+14,-1.572301627719680000e+15,3.454150138396672000e+15,-5.579780992794624000e+15,6.864598984556544000e+15,-6.573052309536768000e+15,4.964023879598080000e+15,-2.978414327758848000e+15,1.423506847825920000e+15,-5.411678925619200000e+14,1.627731551846400000e+14,-3.837084303360000000e+13,6.988974981120000000e+12,-9.639965491200000000e+11,9.790589952000000000e+10,-7.038986240000000000e+09,3.384128000000000000e+08,-9.974272000000000000e+06,1.558480000000000000e+05,-9.680000000000000000e+02,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,
			3.518437208883200000e+13,-4.046202790215680000e+14,2.174833999740928000e+15,-7.257876254949376000e+15,1.684864130613248000e+16,-2.888925570295398400e+16,3.791714811012710400e+16,-3.895882800326246400e+16,3.178220179213516800e+16,-2.075864531468288000e+16,1.089828879020851400e+16,-4.599927086776322000e+15,1.555857691115522000e+15,-4.188847629926410000e+14,8.882601000960000000e+13,-1.461331132416000000e+13,1.826663915520000000e+12,-1.685866291200000000e+11,1.103841024000000000e+10,-4.841408000000000000e+08,1.303456000000000000e+07,-1.862080000000000000e+05,1.058000000000000000e+03,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.000000000000000000e+00,
			1.407374883553280000e+14,-1.688849860263936000e+15,9.499780463984640000e+15,-3.328441599603507200e+16,8.141443799030169600e+16,-1.476820037963612160e+17,2.059929537080197120e+17,-2.260898272405094400e+17,1.981818641905090560e+17,-1.400259325334650880e+17,8.014642191060171200e+16,-3.721787159976344000e+16,1.399977809018879200e+16,-4.246086541639672000e+15,1.030300410839036000e+15,-1.977344222822380000e+14,2.954430332927900000e+13,-3.363677798399000000e+12,2.834209996790000000e+11,-1.697439743900000000e+10,6.820070390000000000e+08,-1.683967900000000000e+07,2.207990000000000000e+05,-1.151000000000000000e+03,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			5.629499534213120000e+14,-7.036874417766400000e+15,4.134163720437760000e+16,-1.517326046330880000e+17,3.900517499535360000e+17,-7.462990149110988800e+17,1.102487181118668800e+18,-1.287455960675123200e+18,1.206989963132928000e+18,-9.175086305116160000e+17,5.688553509172019200e+17,-2.884056849055744000e+17,1.195365667700735840e+17,-4.038397526015998400e+16,1.105751703551999200e+16,-2.432653747814392000e+15,4.248200478719960000e+14,-5.793000652799900000e+13,6.034375680000000000e+12,-4.661529600000000000e+11,2.563841280000000000e+10,-9.472320000000000000e+08,2.152800000000000000e+07,-2.600000000000000000e+05,1.250000000000000000e+03,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,
			0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00
			};	
	
//def log_int_beta_xn(n,mutation_rate):
//  return (np.log(gamma(mutation_rate))+np.log(gamma(n+mutation_rate)) - (np.log(beta(mutation_rate,mutation_rate)) + np.log(gamma(n+2*mutation_rate))))


      private double log_int_beta_xn(int n, double mutation_rate) {

                      double log_gamma = Gamma.logGamma(mutation_rate)+Gamma.logGamma(n+mutation_rate) - (Beta.logBeta( mutation_rate, mutation_rate) + Gamma.logGamma(n+2*mutation_rate));
                              
              
              return log_gamma;
      }
// def Tn_ito_xn_lookup_fast(n,mutation_rate,j):
//  int_sum = D[j,-1]
//  for i in range(0,n):
//      log = np.log(abs(D[j,i])) + (log_int_beta_xn(n-i,mutation_rate))
//      if D[j,i]<=0:
//          int_sum += -np.e**log    
//      if D[j,i]>0:
//          int_sum += np.e**log  
//  return int_sum
      private double[] Tn_ito_xn_lookup_fast() {
    	  if (N != 33) {
    		  throw new IllegalArgumentException("N must be 33");
    	  }
              freqs = new ChebyshevPolynomial(N);
              Double [] thetas = m_substitutionmodel.thetaInput.get().getValues();
              double theta = thetas[thetas.length - 1]; // theta for the root  
              //double theta = 1;
              for (int j = 1; j < N; j++) {
              double int_sum = D[32 + (j-1)*33];
              double llog = 0;
              //System.out.println("here");
              for (int i = 0; i < j; i++) {   
                      llog = Math.log(Math.abs(D[i + (j-1)*33])) + log_int_beta_xn((j)-i,theta);
                      //System.out.println((j));                      
                      //System.out.println(i);
                      //System.out.println(log_int_beta_xn((j)-i,theta));
                      //System.out.println(D[i + (j-1)*33]);
                      //System.out.println(llog);
                      if (D[i + (j-1)*33]<=0) {
                              int_sum = int_sum -Math.exp(llog); }    
                      if (D[i + (j-1)*33]>0) {
                              int_sum = int_sum + Math.exp(llog); } 

                              }
              //System.out.println(int_sum);
              freqs.f[j] = int_sum;
              }
              return freqs.f;
      }

    @Override
    public List<Input<?>> listInputs() {
    	List<Input<?>> list =  super.listInputs();
    	if (!Beauti.isInBeauti() && System.getProperty("beast.is.junit.testing") == null) {
    		// do not expose internal likelihoods to BEAUti or junit tests
    		list.add(likelihoodsInput);
    	}
    	return list;
    }

} // class SnapperTreeLikelihood