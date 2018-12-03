package snap;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.DataType;
import snap.Data;
import snap.SNPSequence;



@Description("Alignment based on a filter operation on another alignment")
public class FilteredAlignment extends Data {
    final public Input<String> filterInput = new Input<>("filter", "specifies which of the sites in the input alignment should be selected " +
            "First site is 1." +
            "Filter specs are comma separated, either a singleton, a range [from]-[to] or iteration [from]:[to]:[step]; " +
            "1-100 defines a range, " +
            "1-100\3 or 1:100:3 defines every third in range 1-100, " +
            "1::3,2::3 removes every third site. " +
            "Default for range [1]-[last site], default for iterator [1]:[last site]:[1]", Validate.REQUIRED);
    final public Input<Alignment> alignmentInput = new Input<>("data", "alignment to be filtered", Validate.REQUIRED);
    final public Input<IntegerParameter> constantSiteWeightsInput = new Input<>("constantSiteWeights", "if specified, constant " +
    		"sites will be added with weights specified by the input. The dimension and order of weights must match the datatype. " +
    		"For example for nucleotide data, a 4 dimensional " +
    		"parameter with weights for A, C, G and T respectively need to be specified.");

    // these triples specify a range for(i=From; i <= To; i += Step)
    int[] from;
    int[] to;
    int[] step;
    /**
     * list of indices filtered from input alignment *
     */
    int[] filter;
    
    boolean convertDataType = false;

    public FilteredAlignment() {
        sequenceInput.setRule(Validate.OPTIONAL);
        // it does not make sense to set weights on sites, since they can be scrambled by the filter
        siteWeightsInput.setRule(Validate.FORBIDDEN);
    }

    @Override
    public void initAndValidate() {
        parseFilterSpec();
        calcFilter();
        Alignment data = alignmentInput.get();
        m_dataType = data.getDataType();
        // see if this filter changes data type
        if (userDataTypeInput.get() != null) {
            m_dataType = userDataTypeInput.get();
            convertDataType = true;
        }

        if (constantSiteWeightsInput.get() != null) {
        	if (constantSiteWeightsInput.get().getDimension() != m_dataType.getStateCount()) {
        		throw new IllegalArgumentException("constantSiteWeights should be of the same dimension as the datatype " +
        				"(" + constantSiteWeightsInput.get().getDimension() + "!="+ m_dataType.getStateCount() +")");
        	}
    	}
        
        taxaNames = data.getTaxaNames();
        stateCounts = data.getStateCounts();
        if (convertDataType && m_dataType.getStateCount() > 0) {
        	for (int i = 0; i < stateCounts.size(); i++) {
                stateCounts.set(i, m_dataType.getStateCount());
        	}
        }

        if (alignmentInput.get().siteWeightsInput.get() != null) {
    		String str = alignmentInput.get().siteWeightsInput.get().trim();
    		String [] strs = str.split(",");
    		siteWeights = new int[strs.length];
    		for (int i = 0; i< strs.length; i++) {
    			siteWeights[i] = Integer.parseInt(strs[i].trim());
    		}    		
        }

        calcPatterns();
        //setupAscertainment();
    }

    private void parseFilterSpec() {
        // parse filter specification
        String filterString = filterInput.get();
        String[] filters = filterString.split(",");
        from = new int[filters.length];
        to = new int[filters.length];
        step = new int[filters.length];
        for (int i = 0; i < filters.length; i++) {
            filterString = " " + filters[i] + " ";
            if (filterString.matches(".*-.*")) {
                // range, e.g. 1-100/3
                if (filterString.indexOf('\\') >= 0) {
                	String str2 = filterString.substring(filterString.indexOf('\\') + 1); 
                	step[i] = parseInt(str2, 1);
                	filterString = filterString.substring(0, filterString.indexOf('\\'));
                } else {
                	step[i] = 1;
                }
                String[] strs = filterString.split("-");
                from[i] = parseInt(strs[0], 1) - 1;
                to[i] = parseInt(strs[1], alignmentInput.get().getSiteCount()) - 1;
            } else if (filterString.matches(".*:.*:.+")) {
                // iterator, e.g. 1:100:3
                String[] strs = filterString.split(":");
                from[i] = parseInt(strs[0], 1) - 1;
                to[i] = parseInt(strs[1], alignmentInput.get().getSiteCount()) - 1;
                step[i] = parseInt(strs[2], 1);
            } else if (filterString.trim().matches("[0-9]*")) {
                from[i] = parseInt(filterString.trim(), 1) - 1;
                to[i] = from[i];
            	step[i] = 1;
            } else {
                throw new IllegalArgumentException("Don't know how to parse filter " + filterString);
            }
        }
    }

    int parseInt(String str, int defaultValue) {
        str = str.replaceAll("\\s+", "");
        try {
            return Integer.parseInt(str);
        } catch (Exception e) {
            return defaultValue;
        }
    }

    private void calcFilter() {
        boolean[] isUsed = new boolean[alignmentInput.get().getSiteCount()];
        for (int i = 0; i < to.length; i++) {
            for (int k = from[i]; k <= to[i]; k += step[i]) {
                isUsed[k] = true;
            }
        }
        // count
        int k = 0;
        for (int i = 0; i < isUsed.length; i++) {
            if (isUsed[i]) {
                k++;
            }
        }
        // set up index set
        filter = new int[k];
        k = 0;
        for (int i = 0; i < isUsed.length; i++) {
            if (isUsed[i]) {
                filter[k++] = i;
            }
        }
    }

    
    @Override
    public List<List<Integer>> getCounts() {
    	if (counts == null) {
			counts = new ArrayList<>();
			for (int j = 0; j < getTaxonCount(); j++) {
				counts.add(new ArrayList<>());
			}
			
	        int nrOfSites = filter.length;
	        for (int i = 0; i < getTaxonCount(); i++) {
	        	List<Integer> seq = alignmentInput.get().getCounts().get(i);
		        for (int j = 0; j < nrOfSites; j++) {	            
	    			counts.get(i).add(seq.get(filter[j]));
		        }
	        }
//	                if (convertDataType) {
//	                	try {
//	                		String code = baseType.getCharacter(data[j][i]);
//							data[j][i] = m_dataType.stringToEncoding(code).get(0);
//	                	} catch (Exception e) {
//	                		e.printStackTrace();
//	                	}
//	                }
//	            }

    	}
    	return counts;
    }
    
    
    @Override
    protected void calcPatterns() {
        // determine # lineages for each site for each taxon
		nrOfLineages = new ArrayList<List<Integer>>();
        counts = null;
   		counts = getCounts();
		//sortByTaxonName(seqs);
		
        for (Sequence seq : alignmentInput.get().sequenceInput.get()) {
        	if (seq instanceof SNPSequence) {
        		// it can vary over sites
        		try {
        			nrOfLineages.add(((SNPSequence) seq).getLineageCounts(m_dataType));
        		} catch (Exception e) {
					e.printStackTrace();
				}
        	} else {
        		// it is constant over all sites 
        		List<Integer> statecounts = new ArrayList<Integer>();
        		int lineageCount = stateCounts.get(nrOfLineages.size()) - 1;
        		int nSites = counts.get(0).size();
        		for (int i = 0; i < nSites; i++) {
        			statecounts.add(lineageCount);
        		}
        		nrOfLineages.add(statecounts);
        	}
        }

		// remove constant sites
		int nTaxa = counts.size();
		int nZeroSitesCount = 0;
		int nAllSitesCount = 0;
		for (int i = 0; i < counts.get(0).size(); i++) {
			// TODO: isConstant should return ALL_REDS ALL_GREEN or MIXED
			if (isConstant(i)) {
				if (counts.get(0).get(i) == 0) {
					nZeroSitesCount++;
				} else {
					nAllSitesCount++;
				}
				for (int j = 0; j < nTaxa; j++) {
					counts.get(j).remove(i);
					nrOfLineages.get(j).remove(i);
				}
				if (siteWeights != null) {
					int [] tmp = new int[siteWeights.length-1];
					System.arraycopy(siteWeights, 0, tmp, 0, i - 1);
					System.arraycopy(siteWeights, i + 1, tmp, i, siteWeights.length - i);
					siteWeights = tmp;
				}
				i--;
			}
		}

		// remove sites that have no data in some branches
		int removed = 0;
		int weightRemoved = 0;
		for (int i = 0; i < counts.get(0).size(); i++) {
			boolean hasZeroCount = false;
			for (int j = 0; j < nTaxa; j++) {
				if (nrOfLineages.get(j).get(i) == 0) {
					hasZeroCount = true;
				}					
			}
			if (hasZeroCount) {
				for (int j = 0; j < nTaxa; j++) {
					counts.get(j).remove(i);
					nrOfLineages.get(j).remove(i);
				}
				removed++;
				if (siteWeights != null) {
				    weightRemoved += siteWeights[i];
					int [] tmp = new int[siteWeights.length-1];
					System.arraycopy(siteWeights, 0, tmp, 0, i - 1);
					System.arraycopy(siteWeights, i + 1, tmp, i, siteWeights.length - i);
					siteWeights = tmp;
				} else {
				    weightRemoved ++;
				}
				i--;
			}
		}
		if (removed > 0) {
			System.out.println("WARNING: removed " + removed + " patterns (" + weightRemoved + " sites) because they have one or more branches without data.");
		}
		
		// find unique patterns
		int nSites = counts.get(0).size();
		int [] weights = new int[nSites];
		for (int i = 0; i < weights.length; i++) {
			int j = 0;
			for (j = 0; j < i; j++) {
				if (isEqual(i,j)) {
					break;
				}
			}
			weights[j]++;
		}
		// count nr of patterns
		int nPatterns = 0;
		for (int i = 0; i < weights.length; i++) {
			if (weights[i]>0) {
				nPatterns++;
			}
		}		
		patternWeight = new int[nPatterns+2];
		sitePatterns = new int[nPatterns+2][nTaxa];
		m_nPatternLineageCounts = new int[nPatterns+2][nTaxa];
//		m_nPatterns = new int[nPatterns][nTaxa];
		patternIndex = new int[nSites];

		nPatterns = 0;
		int iSite = 0;
		// instantiate patterns
		for (int i = 0; i < nSites; i++) {
			if (weights[i]>0) {
				patternWeight[nPatterns] = weights[i];
				for (int j = 0; j < nTaxa; j++) {
					sitePatterns[nPatterns][j] = counts.get(j).get(i);
					m_nPatternLineageCounts[nPatterns][j] = nrOfLineages.get(j).get(i);
				}
				for (int k = 0; k < weights[i]; k++) {
					patternIndex[iSite++] = nPatterns;
				}
				nPatterns++;
			}
		}
        
        Arrays.fill(patternIndex, -1);
        for (int i = 0; i < nSites; i++) {
            int[] sites = new int[nTaxa];
            for (int j = 0; j < nTaxa; j++) {
                sites[j] = counts.get(j).get(i);
            }
            for (int j = 0; j < nPatterns; j++) {
            	boolean found = true;
            	for (int k = 0; k < nTaxa; k++) {
            		if (sites[k] != counts.get(k).get(j)) {
            			found = false;
            			break;
            		}
            	}
            	if (found) {
            		patternIndex[i] = j;
            		j = nPatterns;
            	}
            }
        }
		
        if (siteWeights != null) {
        	Arrays.fill(patternWeight, 0);
            for (int i = 0; i < nSites; i++) {
            	patternWeight[patternIndex[i]] += siteWeights[i];
            }        	
        }
		
		maxStateCount = 0;
		for (int i = 0; i < stateCounts.size(); i++) {
			maxStateCount = Math.max(maxStateCount, stateCounts.get(i));
		}
		// add one for the zero state
		maxStateCount++;
		// report some statistics
		for (int i = 0; i < taxaNames.size(); i++) {
			System.err.println(taxaNames.get(i) + ": " + counts.get(i).size() + " " + stateCounts.get(i));
		}
		
		
		// add dummy patterns
//		TODO: set up weights for dummy patterns
		patternWeight[nPatterns] = nZeroSitesCount;
		patternWeight[nPatterns+1] =nAllSitesCount;
		for (int i = 0; i < nTaxa; i++) {
			int lineageCount = stateCounts.get(i) - 1;
			sitePatterns[nPatterns + 1][i]             = lineageCount;
			m_nPatternLineageCounts[nPatterns][i]     = lineageCount;
			m_nPatternLineageCounts[nPatterns + 1][i] = lineageCount;
		}
//		for (int i = 0; i < nTaxa; i++) {
//			for (int j = 0; j < nPatterns + 2; j++) {
//				System.out.print(m_nPatterns[j][i]);
//			}
//			System.out.println();
//			for (int j = 0; j < nPatterns + 2; j++) {
//				System.out.print(m_nPatternLineageCounts[j][i]);
//			}
//			System.out.println();
//		}

        int totalWeight = 0;
        for (int weight : patternWeight) {
        	totalWeight += weight;
        }
        
        System.out.println(getNrTaxa() + " taxa");
        System.out.println(getSiteCount() + " sites" + (totalWeight == getSiteCount() ? "" : " with weight " + totalWeight));
        System.out.println(getPatternCount() + " patterns");
    }

    /** return indices of the sites that the filter uses **/
    public int [] indices() {
    	return filter.clone();
    }
}
