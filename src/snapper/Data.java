/*
 * File Data.java
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.datatype.Binary;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.datatype.Nucleotide;
import beast.evolution.datatype.StandardData;

import beast.app.BeastMCMC;

import snap.SNPSequence;

@Description("Represents sequence data for SnAP analysis. "+
 "The difference with standard sequence data is that constant sites "+
 "are removed + 2 patterns are added at the end representing these "+
 "constant sites, but with zero weight. The likelihood calculator "+
 "deals with these different sites.")
//use snap.Data instead
@Deprecated
public class Data extends snap.Data {
// extends beast.evolution.alignment.Alignment {
	public Input<beast.evolution.alignment.Alignment> m_rawData = new Input<beast.evolution.alignment.Alignment>("rawdata","raw binary sequences");
	public Input<List<TaxonSet>> m_taxonsets = new Input<List<TaxonSet>>("taxonset","set of taxons that group a number of SNP sequences into a single sequence", 
			new ArrayList<TaxonSet>());
	
	public List<List<Integer>> nrOfLineages;
	public int [][]  m_nPatternLineageCounts;
	private int threadCount;

	public Data() {
		sequenceInput.setRule(Validate.OPTIONAL);
	}
	
	@Override
	public void initAndValidate() {
		
		// guess taxon set if no sequences and no taxonsets are known
		if (/*m_pSequences.get().size() == 0 && */m_taxonsets.get().size() == 0 && m_rawData.get() != null) {
			while (sequenceInput.get().size() > 0) {
				sequenceInput.get().remove(0);
			}
			// by last separator
			int nIgnored = guessTaxonSets("^(.+)[-_\\. ](.*)$", 1);
			if (nIgnored > 0) {
				// by first separator
				nIgnored = guessTaxonSets("^([^-_\\. ]+)[-_\\. ](.*)$", 1);
			}
			if (nIgnored > 0) {
				// by taxon name 
				nIgnored = guessTaxonSets("^(.*)$", 0);
			}
		}
		
		// amalgamate binary sequences into count sequences by taxon sets
		if (m_taxonsets.get().size() > 0) {
			List<Sequence> sequences = m_rawData.get().sequenceInput.get();
			List<Sequence> SNPsequences = sequenceInput.get();
			
			DataType rawDataType = m_rawData.get().getDataType();
			if (rawDataType instanceof Binary || rawDataType instanceof IntegerData || rawDataType instanceof StandardData) {
		
				// sequences are defined through taxon sets, so construct
				// SNPSequences from binary sequences as defined by the taxon sets
				for(TaxonSet set : m_taxonsets.get()) {
					SNPSequence SNPSequence = new SNPSequence();
					SNPSequence.setID(set.getID());
					SNPSequence.taxonInput.setValue(set.getID(), SNPSequence);
					for (Taxon taxon : set.taxonsetInput.get()) {
						boolean bFound = false;
						for (int i = 0; i < sequences.size() && !bFound; i++) {
							if (sequences.get(i).taxonInput.get().equals(taxon.getID())) {
								SNPSequence.m_sequences.setValue(sequences.get(i), SNPSequence);
								bFound = true;
							}
						}
						if (!bFound) {
							String seq = m_rawData.get().defaultInput.get().get(taxon.getID());
							if (seq != null) {
								Sequence sequence = new Sequence(taxon.getID(), seq);
								sequence.totalCountInput.setValue(2, sequence);
								SNPSequence.m_sequences.setValue(sequence, SNPSequence);
							} else {
								throw new IllegalArgumentException("Could not find taxon " + taxon.getID() + " in alignment");
							}
						}
					}
					SNPSequence.initAndValidate();
					SNPsequences.add(SNPSequence);
				}
			} else if (rawDataType instanceof Nucleotide) {
				// call SNPs by setting first sequence to all zero
				// any character in subsequent sequences that is not the same will be set to one
				String refferenceSeq = sequences.get(0).dataInput.get().replaceAll("\\s", "");
				
				for(TaxonSet set : m_taxonsets.get()) {
					SNPSequence SNPSequence = new SNPSequence();
					SNPSequence.setID(set.getID());
					SNPSequence.taxonInput.setValue(set.getID(), SNPSequence);
					for (Taxon taxon : set.taxonsetInput.get()) {
						boolean bFound = false;
						for (int i = 0; i < sequences.size() && !bFound; i++) {
							if (sequences.get(i).taxonInput.get().equals(taxon.getID())) {
								Sequence binarySequence = toBinarySequence(taxon.getID(), refferenceSeq, sequences.get(i).dataInput.get());
								SNPSequence.m_sequences.setValue(binarySequence, SNPSequence);
								bFound = true;
							}
						}
						if (!bFound) {
							String seq = m_rawData.get().defaultInput.get().get(taxon.getID());
							if (seq != null) {
								Sequence binarySequence = toBinarySequence(taxon.getID(), refferenceSeq, seq);
								SNPSequence.m_sequences.setValue(binarySequence, SNPSequence);
							} else {
								throw new IllegalArgumentException("Could not find taxon " + taxon.getID() + " in alignment");
							}
						}
					}
					SNPSequence.initAndValidate();
					SNPsequences.add(SNPSequence);
				}
			} else {
				throw new RuntimeException("Cannot handle data of type " + rawDataType.getTypeDescription() +
						". Use binary or nucleotide data instead.");
			}
		}
		String oldSiteWeights = siteWeightsInput.get(); 
		if (oldSiteWeights == null && m_rawData.get() != null) {
			siteWeightsInput.setValue(m_rawData.get().siteWeightsInput.get(), this);
		}
		
		//System.out.println("Sequence input");
		//System.out.println(sequenceInput.get());
		
		findDataTypes();
		super.initAndValidate();

		threadCount = BeastMCMC.m_nThreads;

        if (threadCount<=1){
			calcPatterns();
			//if (m_rawData.get() != null) {
			//	sequenceInput.get().clear();
			//	siteWeightsInput.setValue(oldSiteWeights, this);
			//}
			
		}
		
		
	} // initAndValidate
	
	private Sequence toBinarySequence(String id, String refferenceSeq,
			String seqStr) {
		Sequence binarySequence = new Sequence();
		binarySequence.setID(id);
		seqStr = seqStr.replaceAll("\\s", "");
		StringBuilder binStr = new StringBuilder();
		for (int i = 0; i < refferenceSeq.length(); i++) {
			if (refferenceSeq.charAt(i) == seqStr.charAt(i)) {
				binStr.append('0');
			} else {
				binStr.append('1');
			}
		}
		binarySequence.initByName("value", binStr.toString(),
				"taxon", id,
				"totalcount", 2);
		return binarySequence;
	}

	/** guesses taxon sets based on pattern in sRegExp based
	 * on the taxa in m_rawData 
	 */
	public int guessTaxonSets(String sRegexp, int nMinSize) {
		m_taxonsets.get().clear();
		List<Taxon> taxa = new ArrayList<Taxon>();
		for (String name : m_rawData.get().getTaxaNames()) { //Sequence sequence : m_rawData.get().sequenceInput.get()) {
			Taxon taxon = new Taxon();
			// ensure sequence and taxon do not get same ID
			//if (sequence.getID() == null || sequence.getID().equals(sequence.taxonInput.get())) {
			//	sequence.setID("_"+sequence.getID());
			//}
			taxon.setID(name); //sequence.taxonInput.get());
			taxa.add(taxon);
		}
		HashMap<String, TaxonSet> map = new HashMap<String, TaxonSet>();
    	Pattern m_pattern = Pattern.compile(sRegexp);
    	int nIgnored = 0;
    	for (Taxon taxon : taxa) {
    		Matcher matcher = m_pattern.matcher(taxon.getID());
			if (matcher.find()) {
				String sMatch = matcher.group(1);
				try {
					if (map.containsKey(sMatch)) {
						TaxonSet set = map.get(sMatch);
						set.taxonsetInput.setValue(taxon, set);
					} else {
						TaxonSet set = new TaxonSet();
						if (sMatch.equals(taxon.getID())) {
							set.setID(sMatch + "_");
						} else {
							set.setID(sMatch);
						}
						set.taxonsetInput.setValue(taxon, set);
						map.put(sMatch, set);
					}
				} catch (Exception ex) {
					ex.printStackTrace();
				}
			} else {
		    	nIgnored++;
			}
    	}
    	// add taxon sets
    	for (TaxonSet set : map.values()) {
    		if (set.taxonsetInput.get().size() > nMinSize) {
                m_taxonsets.setValue(set, this);
    		} else {
    			nIgnored += set.taxonsetInput.get().size();
    		}
    	}
    	return nIgnored;
	}
	
	public int getPatternWeight(int id) {
		if (id < patternWeight.length) {
			return patternWeight[id];
		}
		return 0;
	}

	/** check whether a pattern is all red or all green **/
	// TODO: needs to return whether patterns are all red or all green
	// to deal with zero nrOfLineages when all sites are missing for a 
	// species 
	boolean isConstant(int iSite) {
		int nTaxa = counts.size();
		boolean bAllZero = true;
		boolean bAllMax = true;
		for (int i = 0; i < nTaxa; i++) {
			int iValue = counts.get(i).get(iSite);
			if (iValue > 0) {
				bAllZero = false;
				if (bAllMax == false) {
					return false;
				}
			}
//			if (iValue < m_nStateCounts.get(i)) {
			if (iValue < nrOfLineages.get(i).get(iSite)) {
				bAllMax = false;
				if (bAllZero == false) {
					return false;
				}
			}
		}
		return true;
	} // isConstant
	
	/** calculate patterns from sequence data
	 * The difference with standard sequence data is that constant sites
	 * are removed + 2 patterns are added at the end representing these
	 * constant sites, but with zero weight. The likelihood calculator
	 * deals with these different sites.
	 * 
	 * **/
	@Override
	public void calcPatterns() {
        // determine # lineages for each site for each taxon
		nrOfLineages = new ArrayList<List<Integer>>();
		List<Sequence> seqs = sequenceInput.get();
		//sortByTaxonName(seqs);
		
        for (Sequence seq : seqs) {
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
		//for (int i = 0; i < nTaxa; i++) {
		//	for (int j = 0; j < nPatterns + 2; j++) {
		//		System.out.print(m_nPatterns[j][i]);
		//	}
		//	System.out.println();
		//	for (int j = 0; j < nPatterns + 2; j++) {
		//		System.out.print(m_nPatternLineageCounts[j][i]);
		//	}
		//	System.out.println();
		//}	

        int totalWeight = 0;
        for (int weight : patternWeight) {
        	totalWeight += weight;
        }
        
        System.out.println(getNrTaxa() + " taxa");
        System.out.println(getSiteCount() + " sites" + (totalWeight == getSiteCount() ? "" : " with weight " + totalWeight));
        System.out.println(getPatternCount() + " patterns");
	} // calc


	/** test whether two columns (e.g. sites) contain equal values **/
	protected boolean isEqual(int iSite1, int iSite2) {
		for (int i = 0; i < counts.size(); i++) {
			if (counts.get(i).get(iSite1)
					!= counts.get(i).get(iSite2)) {
				return false;
			}
			if (nrOfLineages.get(i).get(iSite1)
					!= nrOfLineages.get(i).get(iSite2)) {
				return false;
			}
		}
		return true;
	} // isEqual

    public int[] getPatternLineagCounts(int id) {
        return m_nPatternLineageCounts[id];
    }
    
    public double getProportionZeros() {
    	int zeroCount = 0;
    	int oneCount  = 0;
    	for (int i = 0; i < sitePatterns.length; i++) {
    		int [] p = getPattern(i);
    		int [] n = getPatternLineagCounts(i);
    		int w = getPatternWeight(i);
    		for (int j = 0; j < p.length; j++) {
    			oneCount += w * p[j];
    			zeroCount += w * (n[j] - p[j]);
    		}
    	}
    	return (double) zeroCount / ((double) zeroCount + oneCount);
    }

}
