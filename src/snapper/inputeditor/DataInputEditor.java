package snapper.inputeditor;




import java.util.*;
import java.util.regex.PatternSyntaxException;


import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.GuessPatternDialog;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.Alert;
import beastfx.app.util.FXUtils;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.EventHandler;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.SelectionMode;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.control.TableColumn.CellEditEvent;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.VBox;
import snapper.Data;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.datatype.Nucleotide;



public class DataInputEditor extends InputEditor.Base {
	public DataInputEditor(BeautiDoc doc) {
		super(doc);
	}

	List<TaxonSet> m_taxonset;
	List<Taxon> m_lineageset;
	Map<String,String> m_taxonMap;
	TableView<TaxonMap> m_table;
    ObservableList<TaxonMap> taxonMapping;
	//DefaultTableModel m_model = new DefaultTableModel();
	
    public class TaxonMap {
		String taxon;
    	String taxon2;

    	TaxonMap(String taxon, String taxon2) {
    		this.taxon = taxon;
    		this.taxon2 = taxon2;
    	}
    	
    	public String getTaxon() {
			return taxon;
		}
		public void setTaxon(String taxon) {
			this.taxon = taxon;
		}
		public String getTaxon2() {
			return taxon2;
		}
		public void setTaxon2(String taxon2) {
			this.taxon2 = taxon2;
		}
    }
	
	TextField filterEntry;
	String m_sFilter = ".*";
	int m_sortByColumn = 0;
	boolean m_bIsAscending = true;
	
	@Override
	public Class<?> type() {
		return snapper.Data.class;
	}

	
	private static boolean firstNucleotideOccurance = false;
	private final static String nucleotideMessage = 
			"A nucleotide alignment was loaded and SNPs will be called by setting the first "
			+ "sequence as all zeros and other sequences as 1 when they differ from the first sequence.\n\n"
			+ "This could lead to problems if you have sites with more than two values and the first sequence "
			+ "contains a low frequency state.\n\n"
			+ "Consider using phrynomics or other tools to preprocess your alignment instead of using the"
			+ "raw nucleotide alignment.";
	
	@Override
	public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpand, boolean bAddButtons) {
		m_input = input;
		
		if (!firstNucleotideOccurance) {
			Data data = (Data) m_input.get();
			Alignment rawdata = data.m_rawData.get();
			if (rawdata.getDataType() instanceof Nucleotide) {
				Platform.runLater(() -> Alert.showMessageDialog(this, nucleotideMessage));
				firstNucleotideOccurance = true;	
			}
		}
		
		m_beastObject = plugin;
		this.itemNr = itemNr;
		List<TaxonSet> taxonset = ((snapper.Data)input.get()).m_taxonsets.get(); 
		getChildren().add(getContent(taxonset));
	}
	
	
	private Pane getContent(List<TaxonSet> taxonset) {
		m_taxonset = taxonset;
		m_taxonMap = new HashMap<String, String>();
		m_lineageset = new ArrayList<Taxon>();
		for (Taxon taxonset2 : m_taxonset) {
			for (Taxon taxon : ((TaxonSet)taxonset2).taxonsetInput.get()) {
				m_lineageset.add(taxon);
				m_taxonMap.put(taxon.getID(), taxonset2.getID());
			}
		}
		
		taxonMapping = FXCollections.observableArrayList();		
		// set up table.
		// special features: background shading of rows
		// custom editor allowing only Date column to be edited.		
		taxonSetToModel();

		
        // set up table.
        // special features: background shading of rows
        // custom editor allowing only Date column to be edited.
        m_table = new TableView<>();        
        m_table.setPrefWidth(1024);
        m_table.setEditable(true);
        m_table.setItems(taxonMapping);

        TableColumn<TaxonMap, String> col1 = new TableColumn<>("Taxon");
        col1.setPrefWidth(500);
        col1.setEditable(false);
        col1.setCellValueFactory(
        	    new PropertyValueFactory<TaxonMap,String>("Taxon")
        	);
        m_table.getColumns().add(col1);
//        col1.getSortNode().setOnMouseClicked(e -> {
//            // The index of the column whose header was clicked
//			int vColIndex = 0;
//            if (vColIndex != m_sortByColumn) {
//                m_sortByColumn = vColIndex;
//                m_bIsAscending = true;
//            } else {
//                m_bIsAscending = !m_bIsAscending;
//            }
//            taxonSetToModel();
//        });

        TableColumn<TaxonMap, String> col2 = new TableColumn<>("Species/Population");
        col2.setPrefWidth(500);
        col2.setEditable(true);
        col2.setCellValueFactory(
        	    new PropertyValueFactory<TaxonMap,String>("Taxon2")
        	);
        col2.setCellFactory(TextFieldTableCell.forTableColumn());
        col2.setOnEditCommit(
                new EventHandler<CellEditEvent<TaxonMap, String>>() {
  					@Override
  					public void handle(CellEditEvent<TaxonMap, String> event) {
  						String newValue = event.getNewValue();
  						TaxonMap location = event.getRowValue();
  						location.taxon2 = newValue;
  						modelToTaxonset();
  						m_table.refresh();
    					validateInput();
  					}
  				}                
            );
//        col2.getSortNode().setOnMouseClicked(e -> {
//                    // The index of the column whose header was clicked
//        			int vColIndex = 1;
//                    if (vColIndex != m_sortByColumn) {
//                        m_sortByColumn = vColIndex;
//                        m_bIsAscending = true;
//                    } else {
//                        m_bIsAscending = !m_bIsAscending;
//                    }
//                    taxonSetToModel();
//            });
        m_table.getColumns().add(col2);
        m_table.getSelectionModel().setSelectionMode(SelectionMode.MULTIPLE);
		
        
        col2.setOnEditCommit(
                new EventHandler<CellEditEvent<TaxonMap, String>>() {
					@Override
					public void handle(CellEditEvent<TaxonMap, String> event) {
						String newValue = event.getNewValue();
						TaxonMap tipDate = event.getRowValue();
						tipDate.setTaxon2(newValue);
						modelToTaxonset();
					}
				}                
            );
//
//		m_table = new JTable(m_model) {
//			private static final long serialVersionUID = 1L;
//
//			// method that induces table row shading 
//			@Override
//			public Component prepareRenderer (TableCellRenderer renderer,int Index_row, int Index_col) {
//				Component comp = super.prepareRenderer(renderer, Index_row, Index_col);
//				//even index, selected or not selected
//				if (isCellSelected(Index_row, Index_col)) {
//					comp.setBackground(Color.gray);
//				} else 	if (Index_row % 2 == 0) {
//					comp.setBackground(new Color(237,243,255));
//				} else {
//					comp.setBackground(Color.white);
//				}
//				return comp;
//			}
//		};
//		
//		// set up editor that makes sure only doubles are accepted as entry
//		// and only the Date column is editable.
//		m_table.setDefaultEditor(Object.class, new TableCellEditor() {
//			JTextField m_textField = new JTextField();
//			int m_iRow, m_iCol;
//			@Override
//			public boolean stopCellEditing() {
//				m_table.removeEditor();
//				String sText = m_textField.getText();
//				System.err.println(sText);
//				m_model.setValueAt(sText, m_iRow, m_iCol);
////				try {
////					Double.parseDouble(sText);
////				} catch (Exception e) {
////					return false;
////				}
//				modelToTaxonset();
//				return true;
//			}
//		
//			@Override
//			public boolean isCellEditable(EventObject anEvent) {
//				return m_table.getSelectedColumn() == 1;
//			}
//			
//			
//			@Override
//			public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int iRow, int iCol) {
//				if (!isSelected) {
//					return null;
//				}
//				m_iRow = iRow;
//				m_iCol = iCol;
//				m_textField.setText((String)value);
//				return m_textField; 			
//			}
//
//			@Override
//			public boolean shouldSelectCell(EventObject anEvent) {return false;}
//			@Override
//			public void removeCellEditorListener(CellEditorListener l) {}
//			@Override
//			public Object getCellEditorValue() {return null;}
//			@Override
//			public void cancelCellEditing() {}
//			@Override
//			public void addCellEditorListener(CellEditorListener l) {}
//		
//		});				
//		m_table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
//		m_table.getColumnModel().getColumn(0).setPreferredWidth(250);
//		m_table.getColumnModel().getColumn(1).setPreferredWidth(250);
//		
//		JTableHeader header = m_table.getTableHeader();
//		header.addMouseListener(new ColumnHeaderListener());

//		JScrollPane pane = new JScrollPane(m_table);
//		Box tableBox = Box.createHorizontalBox();
//		tableBox.add(Box.createHorizontalGlue());
//		tableBox.add(pane);
//		tableBox.add(Box.createHorizontalGlue());

		VBox box = FXUtils.newVBox();
		box.getChildren().add(createFilterBox());
		box.getChildren().add(m_table);
		box.getChildren().add(createButtonBox());
		return box;
	}
	
	private HBox createButtonBox() {
		HBox buttonBox = FXUtils.newHBox();
		
		Button fillDownButton = new Button("Fill down");
		fillDownButton.setTooltip(new Tooltip("replaces all taxons in selection with the one that is selected at the top"));
		fillDownButton.setOnAction(e-> {
            List<Integer> rows = m_table.getSelectionModel().getSelectedIndices();
			if (rows.size() < 2) {
				return;
			}
            String taxon = taxonMapping.get(rows.get(0)).getTaxon2();
            for (int i = 1; i < rows.size(); i++) {
                taxonMapping.get(rows.get(i)).setTaxon2(taxon);
            }
			modelToTaxonset();
			taxonSetToModel();
		});
		
		Button guessButton = new Button("Guess");
		guessButton.setOnAction(e->guess());
		
		buttonBox.getChildren().add(fillDownButton);
		buttonBox.getChildren().add(guessButton);
		return buttonBox;
	}
	
	

//	public class ColumnHeaderListener extends MouseAdapter {
//	    public void mouseClicked(MouseEvent evt) {
//	        // The index of the column whose header was clicked
//	        int vColIndex = m_table.getColumnModel().getColumnIndexAtX(evt.getX());
//	        if (vColIndex == -1) {
//	            return;
//	        }
//	        if (vColIndex != m_sortByColumn) {
//	        	m_sortByColumn = vColIndex;
//	        	m_bIsAscending = true; 
//	        } else { 
//	        	m_bIsAscending = !m_bIsAscending;
//	        }
//	        taxonSetToModel();
//	    }
//	}

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
	
	private HBox createFilterBox() {
		HBox filterBox = FXUtils.newHBox();
		filterBox.getChildren().add(new Label("filter: "));
		//Dimension size = new Dimension(100,20);
		filterEntry = new TextField();
		filterEntry.setMinWidth(20*12);
		//filterEntry.setColumns(20);
//		filterEntry.setMinimumSize(size);
//		filterEntry.setPreferredSize(size);
//		filterEntry.setSize(size);
		filterEntry.setTooltip(new Tooltip("Enter regular expression to match taxa"));
		filterEntry.setMaxSize(1024, 20);
		filterBox.getChildren().add(filterEntry);
		filterBox.setOnKeyReleased(e-> {
				String sFilter = ".*" + filterEntry.getText() + ".*";
				try {
					// sanity check: make sure the filter is legit
					sFilter.matches(sFilter);
					m_sFilter = sFilter;
					taxonSetToModel();
					m_table.refresh();
				} catch (PatternSyntaxException ex) {
					// ignore
				}
		});
		return filterBox;
	}


	/** for convert taxon sets to table model **/
    private void taxonSetToModel() {
    	
		// count number if lineages that match the filter
		int i = 0;
		for (String sLineageID : m_taxonMap.keySet()) {
			if (sLineageID.matches(m_sFilter)) {
				i++;
			}
		}
		
		// clear table model
		taxonMapping.clear();

		// fill table model with lineages matching the filter
		for (String sLineageID : m_taxonMap.keySet()) {
			if (sLineageID.matches(m_sFilter)) {
				taxonMapping.add(new TaxonMap(sLineageID, m_taxonMap.get(sLineageID)));
			}
		}

//	    Collections.sort(taxonMapping, new Comparator<TaxonMap>() {
//			@Override
//			public int compare(TaxonMap v1, TaxonMap v2) {
//		        String o1 = (String) v1.get(m_sortByColumn);
//		        String o2 = (String) v2.get(m_sortByColumn);
//		        if (o1.equals(o2)) {
//			        o1 = (String) v1.get(1 - m_sortByColumn);
//			        o2 = (String) v2.get(1 - m_sortByColumn);
//		        }
//		        if (m_bIsAscending) {
//		        	return o1.compareTo(o2);
//		        } else {
//		        	return o2.compareTo(o1);
//		        }
//			}
//		});
		if (m_table != null) {
			m_table.refresh();
		}
    }



	/** for convert table model to taxon sets **/
	private void modelToTaxonset() {
		// update map
		for (int i = 0; i < taxonMapping.size();i++) {
			String sLineageID = taxonMapping.get(i).taxon;
			String sTaxonSetID = taxonMapping.get(i).taxon2;
			
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
