package beast.app.draw;



import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.regex.PatternSyntaxException;

import javax.swing.*;
import javax.swing.event.CellEditorListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;

import beast.app.beauti.BeautiDoc;
import beast.app.beauti.GuessPatternDialog;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;



public class DataInputEditor extends InputEditor.Base {
	public DataInputEditor(BeautiDoc doc) {
		super(doc);
	}

	private static final long serialVersionUID = 1L;
	List<TaxonSet> m_taxonset;
	List<Taxon> m_lineageset;
	Map<String,String> m_taxonMap;
	JTable m_table;
	DefaultTableModel m_model = new DefaultTableModel();
	
	JTextField filterEntry;
	String m_sFilter = ".*";
	int m_sortByColumn = 0;
	boolean m_bIsAscending = true;
	
	@Override
	public Class<?> type() {
		return snapper.Data.class;
	}

	
	@Override
	public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpand, boolean bAddButtons) {
		m_input = input;
		m_beastObject = plugin;
		this.itemNr = itemNr;
		List<TaxonSet> taxonset = ((snapper.Data)input.get()).m_taxonsets.get();
		add(getContent(taxonset));
	}
	
	
	private Component getContent(List<TaxonSet> taxonset) {
		m_taxonset = taxonset;
		m_taxonMap = new HashMap<String, String>();
		m_lineageset = new ArrayList<Taxon>();
		for (Taxon taxonset2 : m_taxonset) {
			for (Taxon taxon : ((TaxonSet)taxonset2).taxonsetInput.get()) {
				m_lineageset.add(taxon);
				m_taxonMap.put(taxon.getID(), taxonset2.getID());
			}
		}
		
		// set up table.
		// special features: background shading of rows
		// custom editor allowing only Date column to be edited.		
		m_model = new DefaultTableModel();
		m_model.addColumn("Taxon");
		m_model.addColumn("Species/Population");
		taxonSetToModel();

		m_table = new JTable(m_model) {
			private static final long serialVersionUID = 1L;

			// method that induces table row shading 
			@Override
			public Component prepareRenderer (TableCellRenderer renderer,int Index_row, int Index_col) {
				Component comp = super.prepareRenderer(renderer, Index_row, Index_col);
				//even index, selected or not selected
				if (isCellSelected(Index_row, Index_col)) {
					comp.setBackground(Color.gray);
				} else 	if (Index_row % 2 == 0) {
					comp.setBackground(new Color(237,243,255));
				} else {
					comp.setBackground(Color.white);
				}
				return comp;
			}
		};
		
		// set up editor that makes sure only doubles are accepted as entry
		// and only the Date column is editable.
		m_table.setDefaultEditor(Object.class, new TableCellEditor() {
			JTextField m_textField = new JTextField();
			int m_iRow, m_iCol;
			@Override
			public boolean stopCellEditing() {
				m_table.removeEditor();
				String sText = m_textField.getText();
				System.err.println(sText);
				m_model.setValueAt(sText, m_iRow, m_iCol);
//				try {
//					Double.parseDouble(sText);
//				} catch (Exception e) {
//					return false;
//				}
				modelToTaxonset();
				return true;
			}
		
			@Override
			public boolean isCellEditable(EventObject anEvent) {
				return m_table.getSelectedColumn() == 1;
			}
			
			
			@Override
			public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int iRow, int iCol) {
				if (!isSelected) {
					return null;
				}
				m_iRow = iRow;
				m_iCol = iCol;
				m_textField.setText((String)value);
				return m_textField; 			
			}

			@Override
			public boolean shouldSelectCell(EventObject anEvent) {return false;}
			@Override
			public void removeCellEditorListener(CellEditorListener l) {}
			@Override
			public Object getCellEditorValue() {return null;}
			@Override
			public void cancelCellEditing() {}
			@Override
			public void addCellEditorListener(CellEditorListener l) {}
		
		});				
		m_table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		m_table.getColumnModel().getColumn(0).setPreferredWidth(250);
		m_table.getColumnModel().getColumn(1).setPreferredWidth(250);
		
		JTableHeader header = m_table.getTableHeader();
		header.addMouseListener(new ColumnHeaderListener());

		JScrollPane pane = new JScrollPane(m_table);
		Box tableBox = Box.createHorizontalBox();
		tableBox.add(Box.createHorizontalGlue());
		tableBox.add(pane);
		tableBox.add(Box.createHorizontalGlue());

		Box box = Box.createVerticalBox();
		box.add(createFilterBox());
		box.add(tableBox);
		box.add(createButtonBox());
		return box;
	}
	
	private Component createButtonBox() {
		Box buttonBox = Box.createHorizontalBox();
		
		JButton fillDownButton = new JButton("Fill down");
		fillDownButton.setToolTipText("replaces all taxons in selection with the one that is selected at the top");
		fillDownButton.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				int [] rows = m_table.getSelectedRows();
				if (rows.length < 2) {
					return;
				}
				String sTaxon = (String) ((Vector<?>) m_model.getDataVector().elementAt(rows[0])).elementAt(1);
				for (int i = 1; i < rows.length; i++) {
					m_model.setValueAt(sTaxon, rows[i], 1);
				}
				modelToTaxonset();
			}
		});
		
		JButton guessButton = new JButton("Guess");
		guessButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				guess();
			}
		});
		
		buttonBox.add(Box.createHorizontalGlue());
		buttonBox.add(fillDownButton);
		buttonBox.add(Box.createHorizontalGlue());
		buttonBox.add(guessButton);
		buttonBox.add(Box.createHorizontalGlue());
		return buttonBox;
	}
	
	

	public class ColumnHeaderListener extends MouseAdapter {
	    public void mouseClicked(MouseEvent evt) {
	        // The index of the column whose header was clicked
	        int vColIndex = m_table.getColumnModel().getColumnIndexAtX(evt.getX());
	        if (vColIndex == -1) {
	            return;
	        }
	        if (vColIndex != m_sortByColumn) {
	        	m_sortByColumn = vColIndex;
	        	m_bIsAscending = true; 
	        } else { 
	        	m_bIsAscending = !m_bIsAscending;
	        }
	        taxonSetToModel();
	    }
	}

	void parseTrait(String trait) {
    	Map<String,String> traitmap = new HashMap<String, String>();
    	for (String line : trait.split(",")) {
    		String [] strs = line.split("=");
    		if (strs.length == 2) {
    			traitmap.put(strs[0].trim(), strs[1].trim());
    		}
    	}

    	m_taxonset.clear();
		List<Taxon> taxa = new ArrayList<Taxon>();
        for (Alignment alignment : getDoc().alignments) {
            for (Sequence sequence : alignment.sequenceInput.get()) {
				Taxon taxon = new Taxon();
				// ensure sequence and taxon do not get same ID
				if (sequence.getID() == null || sequence.getID().equals(sequence.taxonInput.get())) {
					sequence.setID("_"+sequence.getID());
				}
				taxon.setID(sequence.taxonInput.get());
				taxa.add(taxon);
            }
		}
		HashMap<String, TaxonSet> map = new HashMap<String, TaxonSet>();
        for (Taxon taxon : taxa) {
            if (!(taxon instanceof TaxonSet)) {
                String sMatch = traitmap.get(taxon.getID());
                if (sMatch != null) {
                    try {
                        if (map.containsKey(sMatch)) {
                            TaxonSet set = map.get(sMatch);
                            set.taxonsetInput.setValue(taxon, set);
                        } else {
                            TaxonSet set = new TaxonSet();
                            set.setID(sMatch);
                            set.taxonsetInput.setValue(taxon, set);
                            map.put(sMatch, set);
                        }
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }
            }
        }
    	// add taxon sets
    	for (TaxonSet set : map.values()) {
            m_taxonset.add(set);
    	}
	}


	String m_sPattern = "^(.+)[-_\\. ](.*)$";
	private void guess() {
        GuessPatternDialog dlg = new GuessPatternDialog(this, m_sPattern);
		try {
	        switch (dlg.showDialog("Guess taxon sets")) {
	        case canceled: return;
	        case trait: parseTrait(dlg.getTrait());
	        	break;
	        case pattern:
	            m_sPattern = dlg.getPattern(); 
				((snapper.Data)m_input.get()).guessTaxonSets(m_sPattern, 0);
				
				break;
	        }
			for (Taxon taxonset2 : m_taxonset) {
				for (Taxon taxon : ((TaxonSet)taxonset2).taxonsetInput.get()) {
					m_lineageset.add(taxon);
					m_taxonMap.put(taxon.getID(), taxonset2.getID());
				}
			}
			taxonSetToModel();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private Component createFilterBox() {
		Box filterBox = Box.createHorizontalBox();
		filterBox.add(new JLabel("filter: "));
		//Dimension size = new Dimension(100,20);
		filterEntry = new JTextField();
		filterEntry.setColumns(20);
//		filterEntry.setMinimumSize(size);
//		filterEntry.setPreferredSize(size);
//		filterEntry.setSize(size);
		filterEntry.setToolTipText("Enter regular expression to match taxa");
		filterEntry.setMaximumSize(new Dimension(1024, 20));
		filterBox.add(filterEntry);
		filterBox.add(Box.createHorizontalGlue());
		filterEntry.getDocument().addDocumentListener(new DocumentListener() {
			@Override
			public void removeUpdate(DocumentEvent e) {
				processFilter();
			}
			@Override
			public void insertUpdate(DocumentEvent e) {
				processFilter();
			}
			@Override
			public void changedUpdate(DocumentEvent e) {
				processFilter();
			}
			private void processFilter() {
				String sFilter = ".*" + filterEntry.getText() + ".*";
				try {
					// sanity check: make sure the filter is legit
					sFilter.matches(sFilter);
					m_sFilter = sFilter;
					taxonSetToModel();
					m_table.repaint();
				} catch (PatternSyntaxException e) {
					// ignore
				}
			}
		});
		return filterBox;
	}




	/** for convert taxon sets to table model **/
	@SuppressWarnings("unchecked")
	private void taxonSetToModel() {
		// count number if lineages that match the filter
		int i = 0;
		for (String sLineageID : m_taxonMap.keySet()) {
			if (sLineageID.matches(m_sFilter)) {
				i++;
			}
		}
		
		// clear table model
		while (m_model.getRowCount() > 0) {
			m_model.removeRow(0);
		}	
		
		// fill table model with lineages matching the filter
		for (String sLineageID : m_taxonMap.keySet()) {
			if (sLineageID.matches(m_sFilter)) {
				Object [] rowData = new Object[2];
				rowData[0] = sLineageID;
				rowData[1] = m_taxonMap.get(sLineageID);
				m_model.addRow(rowData);
			}
		}
		
	    @SuppressWarnings("rawtypes")
		Vector data = m_model.getDataVector();
	    Collections.sort(data, new Comparator<Vector<?>>() {
			@Override
			public int compare(Vector<?> v1, Vector<?> v2) {
		        String o1 = (String) v1.get(m_sortByColumn);
		        String o2 = (String) v2.get(m_sortByColumn);
		        if (o1.equals(o2)) {
			        o1 = (String) v1.get(1 - m_sortByColumn);
			        o2 = (String) v2.get(1 - m_sortByColumn);
		        }
		        if (m_bIsAscending) {
		        	return o1.compareTo(o2);
		        } else {
		        	return o2.compareTo(o1);
		        }
			}
	    	
		});
	    m_model.fireTableRowsInserted(0, m_model.getRowCount());
	}

	/** for convert table model to taxon sets **/
	private void modelToTaxonset() {
		// update map
		for (int i = 0; i < m_model.getRowCount();i++) {
			String sLineageID = (String) ((Vector<?>)m_model.getDataVector().elementAt(i)).elementAt(0);
			String sTaxonSetID = (String) ((Vector<?>)m_model.getDataVector().elementAt(i)).elementAt(1);
			
			// new taxon set?
			if (!m_taxonMap.containsValue(sTaxonSetID)) {
				// create new taxon set
				TaxonSet taxonset = new TaxonSet();
				taxonset.setID(sTaxonSetID);
				m_taxonset.add(taxonset);
			}
			m_taxonMap.put(sLineageID, sTaxonSetID);
		}	
		
		// clear old taxon sets
		for (TaxonSet set : m_taxonset) {
			set.taxonsetInput.get().clear();
		}
		
		// group lineages with their taxon sets
		for (String sLineageID : m_taxonMap.keySet()) {
			for (Taxon taxon : m_lineageset) {
				if (taxon.getID().equals(sLineageID)) {
					String sTaxonSet = m_taxonMap.get(sLineageID);
					for (TaxonSet set : m_taxonset) {
						if (set.getID().equals(sTaxonSet)) {
							try {
								set.taxonsetInput.setValue(taxon, set);
							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}
				}
			}
		}
		
		// remove unused taxon sets
		for (int i = m_taxonset.size()-1; i >= 0; i--) {
			if (m_taxonset.get(i).taxonsetInput.get().size() == 0) {
				m_taxonset.remove(i);
			}
		}
	}

}
