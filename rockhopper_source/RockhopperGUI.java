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

import java.util.ArrayList;
import java.util.Scanner;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.DocumentEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.net.URL;
import java.net.URI;
import java.io.IOException;

public class RockhopperGUI extends JApplet implements ActionListener,DocumentListener,MouseListener {

    /*****************************************
     **********   CLASS VARIABLES   **********
     *****************************************/

    public static final String directory = Rockhopper.output_DIR + "intermediary/";
    public static final String fileName = "experiments.txt";

    private static final Color white = Color.white;
    public static final Color lightBlue = new Color(150,200,255);
    public static final Color gold = new Color(255,215,0);
    public static final Color lightPurple = new Color(240,240,255);
    public static ClassLoader cl;
    public static URL penguin;
    private static JFrame frame;



    /********************************************
     **********   INSTANCE VARIABLES   **********
     ********************************************/

    private NCBI ncbi;
    private String genome_DIR;
    private JMenuBar menuBar;
    private JMenuItem parameterMenuItem;
    private JMenuItem cacheMenuItem;
    private JMenuItem exitMenuItem;
    private JMenuItem userGuideMenuItem;
    private JMenuItem aboutMenuItem;

    private JPanel mainPanel;
    private JScrollPane outerRepliconPanel;
    private JPanel innerRepliconPanel;
    private ArrayList<JTextField> genomeTextFields;
    private JPopupMenu popupMenu;
    private int currentReplicon;  // Replicon associated with popup menu
    private ArrayList<IndexedMenuItem> popupReplicons;
    private ArrayList<JLabel> genomeLabels;
    private ArrayList<JLabel> customGenomeLabels;
    private ArrayList<String> customGenomeDirs;
    private JLabel genomeAddLabel;
    private JPanel denovoPanel;
    private JLabel denovoLabel;
    private boolean isDenovo;  // Performing de novo assembly (true) or reference based assembly (false)
    private JPanel denovoInfo;
    private JLabel submitLabel;
    private ImageIcon submitButton, submitButtonPressed, genomeAddButton, genomeAddButtonPressed, genomeRemoveButton, genomeRemoveButtonPressed, customGenomeButton, customGenomeButtonPressed, denovoNoButton, denovoYesButton, denovoNoButtonPressed, denovoYesButtonPressed;
    private JPanel spacer;
    private boolean readingFromFile;
    private JFileChooser fc;

    // RNAseq files
    private ExperimentsGUI exps;
    private ParametersGUI params;
    private ProgressGUI progress;
    private Rockhopper rh;



    /**************************************
     **********   CONSTRUCTORS   **********
     **************************************/

    public RockhopperGUI() {
	cl = this.getClass().getClassLoader();
	penguin = cl.getResource("images/Rockhopper.jpg");
	fc = new JFileChooser();
	fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
	customGenomeDirs = new ArrayList<String>();
	ncbi = new NCBI();
	exps = new ExperimentsGUI(this);
	params = new ParametersGUI(this);
	progress = new ProgressGUI(this);
	makeMenuBar();
	makeMainPanel();
	readExperimentsFromFile();
    }



    /*************************************************
     **********   PUBLIC INSTANCE METHODS   **********
     *************************************************/

    public void init() {
	this.setJMenuBar(menuBar);
	add(mainPanel);
    }

    public void actionPerformed(ActionEvent event) {
	Object source = event.getSource();

	if (source.equals(aboutMenuItem)) newDialogAbout();
	else if (source.equals(parameterMenuItem)) params.displayParameters();
	else if (source.equals(cacheMenuItem)) clearCache();
	else if (source.equals(exitMenuItem)) System.exit(0);
	else if (source.equals(userGuideMenuItem)) openUserGuide();
	else if (params.isSource(source)) params.handleEvent(source);
	else {  // Pop up menu selection
	    // Check if event was one of the replicon text fields and, if so, get index.
	    int z = -1;
	    for (int i=0; i<genomeTextFields.size(); i++) {
		if (source.equals(genomeTextFields.get(i))) z = i;
	    }
	    if (z >= 0) handleGenomeTextField(z);  // Replicon text field
	    else {  // Replicon popup menu
		for (IndexedMenuItem item : popupReplicons) {
		    if (source.equals(item)) {
			if (item.getText().indexOf("*** MORE ***") >= 0) {
			    updateRepliconMenu(currentReplicon, item.index);
			} else {
			    genomeTextFields.get(currentReplicon).setText(item.getText());
			    //removePopupMenu();
			    customGenomeDirs.set(currentReplicon, null);
			    handleGenomeTextField(currentReplicon);
			}
			break;
		    }
		}
	    }
	}
    }

    public void insertUpdate(DocumentEvent event) {
	for (int i=0; i<genomeTextFields.size(); i++) {
	    if (event.getDocument().equals(genomeTextFields.get(i).getDocument()) && !readingFromFile) updateRepliconMenu(i, 0);
	}
    }

    public void removeUpdate(DocumentEvent event) {
	for (int i=0; i<genomeTextFields.size(); i++) {
	    if (event.getDocument().equals(genomeTextFields.get(i).getDocument())) updateRepliconMenu(i, 0);
	}
    }

    public void changedUpdate(DocumentEvent event) {
	// Plain text components do not fire these events.
    }

    public void mouseClicked(MouseEvent event) {
	Object source = event.getSource();
	if (exps.isSource(source)) {
	    mainPanel.remove(exps.getScrollPane());
	    mainPanel.remove(submitLabel);
	    exps.handleEvent(source);
	    updateSeqFilePanel();
	} else if (params.isSource(source)) {
	    params.handleEvent(source);
	} else if (source.equals(genomeAddLabel) && !isDenovo) {  // Genome add button is pressed
	    addReplicon();
	    if (exps.getNumExperiments() > 0) {
		mainPanel.remove(exps.getScrollPane());
		mainPanel.remove(submitLabel);
		updateSeqFilePanel();
	    } else {
		validate();
		repaint();
	    }
	} else if (source.equals(denovoLabel)) {  // Toggle de novo assembly vs. reference based assembly
	    if (isDenovo) {  // Changing from de novo to reference based assembly
		toggleRepliconDeNovoPanels(denovoInfo, innerRepliconPanel);
	    } else {  // Changing from reference based to de novo assembly
		toggleRepliconDeNovoPanels(innerRepliconPanel, denovoInfo);
	    }
	    isDenovo = !isDenovo;
	    mainPanel.remove(exps.getScrollPane());
	    mainPanel.remove(submitLabel);
	    if (exps.getNumExperiments() == 0) exps.addExperiment();
	    updateSeqFilePanel();
	} else if (source.equals(submitLabel)) {  // Submit button is pressed
	    if (!isDenovo) {  // Check replicons if performing reference based assembly
		if (!prepareRepliconFiles()) return;
	    }
	    if (!outputExperimentsToFile()) return;
	    progress.createProgressWindow();
	    Peregrine.output = progress.progressArea;
	    Assembler.output = progress.progressArea;
	    Rockhopper.output = progress.progressArea;
	    Rockhopper.commandLineArguments(generateCommandLineArguments());
	    Peregrine.wigFiles = null;

	    class Worker extends SwingWorker<Void, Void> {
		public Void doInBackground() {
		    Rockhopper.worker = this;
		    rh = new Rockhopper();
		    return null;
		}
		public void done() {
		    setProgress(100);
		    if (!isDenovo) progress.addIGV();  // Add IGV button only if NOT de novo assembly
		    validate();
		    repaint();
		}
	    }
	    Worker w = new Worker();
	    w.addPropertyChangeListener(progress);
	    w.execute();
	} else if (progress.isSource(source)) {
	    if (isDenovo) return;  // No need to launch IGV if completed de novo assembly
	    // Launch IGV!
	    try {
		boolean success = true;
		ArrayList<Genome> genomes = rh.getGenomes();
		File genomeRegistryFile = new File(directory + "genomes.txt");
		success &= IGV_Ops.createGenomeRegistry(genomeRegistryFile, genomes);
		for (int i=0; i<genomes.size(); i++) {
		    Genome genome = genomes.get(i);
		    String baseFileName = genome.getBaseFileName();
		    success &= IGV_Ops.createGenomeArchive(baseFileName + ".genome", genome.getName(), baseFileName + ".fna", baseFileName + ".gff");
		    success &= IGV_Ops.createBatchFileToLoadData(directory + "batch_" + genome.getName().replace(' ', '_') + ".txt", genome.getName(), formatWigFiles(Peregrine.wigFiles.get(i)), getRockhopperWigFiles(genome));
		}
		final String firstBatchFile = directory + "batch_" + genomes.get(0).getName().replace(' ', '_') + ".txt";
		rh = null;
		System.gc();  // Clean up memory before launching IGV

		class Worker extends SwingWorker<Void, Void> {
		    public Void doInBackground() {
			String[] args = {"-u", "file:" + directory + "genomes.txt", "-b", firstBatchFile};
			org.broad.igv.ui.Main.main(args);
			return null;
		    }
		    public void done() { }
		}
		(new Worker()).execute();
		if (!success) error("IGV error", "There was an error launching and initializing IGV for viewing results.");
	    } catch (Exception ex) {
		error("IGV error", "Could not open IGV for viewing results.");
	    }
	} else {  // Check if replicon remove button was clicked
	    for (int i=0; i<genomeLabels.size(); i++) {
		if (source.equals(genomeLabels.get(i)) && !isDenovo) {
		    removeReplicon(i);
		    if (exps.getNumExperiments() > 0) {
			mainPanel.remove(exps.getScrollPane());
			mainPanel.remove(submitLabel);
			updateSeqFilePanel();
		    } else {
			validate();
			repaint();
		    }
		}
	    }
	    for (int i=0; i<customGenomeLabels.size(); i++) {
		if (source.equals(customGenomeLabels.get(i)) && !isDenovo) {
		    addCustomGenome(i);
		}
	    }
	}
    }

    public void mousePressed(MouseEvent event) { 
	Object source = event.getSource();
	exps.mousePressed(source);
	params.mousePressed(source);
	progress.mousePressed(source);
	if (source.equals(genomeAddLabel) && !isDenovo) genomeAddLabel.setIcon(genomeAddButtonPressed);
	if (source.equals(submitLabel)) submitLabel.setIcon(submitButtonPressed);
	for (JLabel genomeLabel : genomeLabels) {
	    if (source.equals(genomeLabel) && !isDenovo) genomeLabel.setIcon(genomeRemoveButtonPressed);
	}
	for (JLabel customGenomeLabel : customGenomeLabels) {
	    if (source.equals(customGenomeLabel) && !isDenovo) customGenomeLabel.setIcon(customGenomeButtonPressed);
	}
	if (source.equals(denovoLabel)) {
	    if (isDenovo) denovoLabel.setIcon(denovoYesButtonPressed);
	    else denovoLabel.setIcon(denovoNoButtonPressed);
	}

    }

    public void mouseReleased(MouseEvent event) { 
	Object source = event.getSource();
	exps.mouseReleased(source);
	params.mouseReleased(source);
	progress.mouseReleased(source);
	if (source.equals(genomeAddLabel) && !isDenovo) genomeAddLabel.setIcon(genomeAddButton);
	if (source.equals(submitLabel)) submitLabel.setIcon(submitButton);
	for (JLabel genomeLabel : genomeLabels) {
	    if (source.equals(genomeLabel) && !isDenovo) genomeLabel.setIcon(genomeRemoveButton);
	}
	for (JLabel customGenomeLabel : customGenomeLabels) {
	    if (source.equals(customGenomeLabel) && !isDenovo) customGenomeLabel.setIcon(customGenomeButton);
	}
	if (source.equals(denovoLabel)) {
	    if (isDenovo) denovoLabel.setIcon(denovoNoButton);
	    else denovoLabel.setIcon(denovoYesButton);
	}
    }

    public void mouseEntered(MouseEvent event) { }

    public void mouseExited(MouseEvent event) { }



    /**************************************************
     **********   PRIVATE INSTANCE METHODS   **********
     **************************************************/

    private void makeMenuBar() {
	menuBar = new JMenuBar();
	menuBar.setOpaque(true);
	menuBar.setBackground(gold);
	menuBar.setPreferredSize(new Dimension(200,20));

	JMenu fileMenu = new JMenu("Options");
	menuBar.add(fileMenu);
	fileMenu.setBackground(gold);
	fileMenu.addSeparator();
	parameterMenuItem = new JMenuItem("Parameter settings");
	parameterMenuItem.setBackground(gold);
	fileMenu.add(parameterMenuItem);
	parameterMenuItem.addActionListener(this);
	fileMenu.addSeparator();
	cacheMenuItem = new JMenuItem("Clear cache");
	cacheMenuItem.setBackground(gold);
	fileMenu.add(cacheMenuItem);
	cacheMenuItem.addActionListener(this);
	fileMenu.addSeparator();
	exitMenuItem = new JMenuItem("Exit");
	exitMenuItem.setBackground(gold);
	fileMenu.add(exitMenuItem);
	exitMenuItem.addActionListener(this);

	JMenu helpMenu = new JMenu("Help");
	menuBar.add(helpMenu);
	helpMenu.setBackground(gold);
	userGuideMenuItem = new JMenuItem("User guide");
	userGuideMenuItem.setBackground(gold);
	helpMenu.add(userGuideMenuItem);
	userGuideMenuItem.addActionListener(this);
	helpMenu.addSeparator();
	aboutMenuItem = new JMenuItem("About Rockhopper");
	aboutMenuItem.setBackground(gold);
	helpMenu.add(aboutMenuItem);
	aboutMenuItem.addActionListener(this);
    }

    private void newDialogAbout() {
	JPanel panel = new JPanel();
	panel.setBackground(lightBlue);
	panel.setLayout(new GridLayout(1,2));
	JLabel picLabel = new JLabel();
	picLabel.setIcon(new ImageIcon(penguin));
	JPanel textPanel = new JPanel();
	textPanel.setBackground(lightBlue);
	textPanel.setLayout(new GridLayout(7,1));
	textPanel.add(new JLabel());
	JLabel textLabel1 = new JLabel("Rockhopper", SwingConstants.CENTER);
	textLabel1.setFont(new Font("Times", Font.ITALIC, 30));
	textPanel.add(textLabel1);
	JLabel textLabel2 = new JLabel("version " + Rockhopper.version, SwingConstants.CENTER);
	textLabel2.setFont(new Font("Times", Font.PLAIN, 16));
	textPanel.add(textLabel2);
	textPanel.add(new JLabel());
	JLabel textLabel3 = new JLabel("A tool for analyzing", SwingConstants.CENTER);
	textLabel3.setFont(new Font("Times", Font.PLAIN, 12));
	textPanel.add(textLabel3);
	JLabel textLabel4 = new JLabel("bacterial RNAseq data", SwingConstants.CENTER);
	textLabel4.setFont(new Font("Times", Font.PLAIN, 12));
	textPanel.add(textLabel4);
	panel.add(picLabel);
	panel.add(textPanel);
	JDialog about = new JDialog(frame, "About Rockhopper", true);
	about.setSize(500,200);
	about.setResizable(false);
	about.setLocationRelativeTo(this);
	about.add(panel);
	about.setVisible(true);
    }

    private void makeMainPanel() {
	mainPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
	mainPanel.setBackground(lightBlue);
	mainPanel.add(getSpacerPanel(790,20));
	genomeAddButton = RockhopperGUI.getScaledImage("images/addGenome.gif", 20, 20);
	genomeAddButtonPressed = RockhopperGUI.getScaledImage("images/addGenomePressed.gif", 20, 20);
	genomeRemoveButton = RockhopperGUI.getScaledImage("images/removeGenome.gif", 20, 20);
	genomeRemoveButtonPressed = RockhopperGUI.getScaledImage("images/removeGenomePressed.gif", 20, 20);
	customGenomeButton = RockhopperGUI.getScaledImage("images/customGenome.gif", 20, 20);
	customGenomeButtonPressed = RockhopperGUI.getScaledImage("images/customGenomePressed.gif", 20, 20);
	denovoNoButton = RockhopperGUI.getScaledImage("images/denovoNoButton.png", 20, 20);
	denovoNoButtonPressed = RockhopperGUI.getScaledImage("images/denovoNoButtonPressed.png", 20, 20);
	denovoYesButton = RockhopperGUI.getScaledImage("images/denovoYesButton.png", 20, 20);
	denovoYesButtonPressed = RockhopperGUI.getScaledImage("images/denovoYesButtonPressed.png", 20, 20);
	isDenovo = false;
	denovoInfo = new JPanel(new GridBagLayout());
	denovoInfo.setBackground(lightBlue);
	denovoInfo.setBorder(BorderFactory.createEmptyBorder());
	JLabel denovoInfoLabel = new JLabel("Perform de novo transcript assembly without reference replicons");
	denovoInfoLabel.setFont(new Font("Times", Font.BOLD, 14));
	denovoInfoLabel.setForeground(gold);
	denovoInfo.add(denovoInfoLabel);
	innerRepliconPanel = new JPanel(new FlowLayout(FlowLayout.LEFT,2,2));
	outerRepliconPanel = new JScrollPane(innerRepliconPanel);
	outerRepliconPanel.setPreferredSize(new Dimension(500, 100));
	outerRepliconPanel.setBackground(lightBlue);
	outerRepliconPanel.setBorder(BorderFactory.createEmptyBorder());
	genomeAddLabel = new JLabel(genomeAddButton);
	genomeAddLabel.setToolTipText("Add replicon");
	genomeAddLabel.addMouseListener(this);
	genomeTextFields = new ArrayList<JTextField>();
	genomeLabels = new ArrayList<JLabel>();
	customGenomeLabels = new ArrayList<JLabel>();
	denovoLabel = new JLabel(denovoNoButton);
	denovoLabel.setToolTipText("Press here to toggle between de novo transcript assembly and reference based transcript assembly");
	denovoLabel.addMouseListener(this);
	denovoPanel = new JPanel();
	denovoPanel.setPreferredSize(new Dimension(20, 100));
	denovoPanel.setBackground(lightBlue);
	denovoPanel.setBorder(BorderFactory.createEmptyBorder());
	denovoPanel.add(denovoLabel);
	addReplicon();
	popupReplicons = new ArrayList<IndexedMenuItem>();
	mainPanel.add(outerRepliconPanel);
	mainPanel.add(denovoPanel);

	submitButton = getScaledImage("images/Submit.gif", 117, 41);
	submitButtonPressed = getScaledImage("images/SubmitPressed.gif", 117, 41);
	submitLabel = new JLabel(submitButton);
	submitLabel.setToolTipText("Submit RNAseq files for analysis");
	submitLabel.addMouseListener(this);
	spacer = getSpacerPanel(790,5);
    }

    public void addReplicon() {
	outerRepliconPanel.remove(innerRepliconPanel);
	genomeTextFields.add(new JTextField("", 25));
	genomeTextFields.get(genomeTextFields.size()-1).setBackground(gold);
	genomeTextFields.get(genomeTextFields.size()-1).addActionListener(this);
	genomeTextFields.get(genomeTextFields.size()-1).getDocument().addDocumentListener(this);
	genomeLabels.add(new JLabel(genomeRemoveButton));
	genomeLabels.get(genomeLabels.size()-1).setToolTipText("Remove replicon");
	genomeLabels.get(genomeLabels.size()-1).addMouseListener(this);
	customGenomeLabels.add(new JLabel(customGenomeButton));
	customGenomeLabels.get(customGenomeLabels.size()-1).setToolTipText("Use a replicon whose sequence is not available from GenBank");
	customGenomeLabels.get(customGenomeLabels.size()-1).addMouseListener(this);
	customGenomeDirs.add(null);
	currentReplicon = genomeTextFields.size()-1;
	populateOuterRepliconPanel();
    }

    public void removeReplicon(int i) {
	outerRepliconPanel.remove(innerRepliconPanel);
	genomeTextFields.remove(i);
	genomeLabels.remove(i);
	customGenomeLabels.remove(i);
	customGenomeDirs.remove(i);
	populateOuterRepliconPanel();
    }

    private JPanel generateRepliconPanel(int i) {
	JPanel newPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
	newPanel.setPreferredSize(new Dimension(480, 27));
	newPanel.setBackground(lightBlue);
	newPanel.add(new JLabel("Replicon name: "));
	newPanel.add(genomeTextFields.get(i));
	newPanel.add(genomeLabels.get(i));
	newPanel.add(customGenomeLabels.get(i));
	JLabel spacer = new JLabel("");
	spacer.setPreferredSize(new Dimension(70, 0));
	newPanel.add(spacer);
	return newPanel;
    }

    private void populateOuterRepliconPanel() {
	innerRepliconPanel = new JPanel(new FlowLayout(FlowLayout.LEFT,2,2));
	innerRepliconPanel.setPreferredSize(new Dimension(480, genomeTextFields.size()*29 + 40));
	innerRepliconPanel.setBackground(lightBlue);
	innerRepliconPanel.setBorder(BorderFactory.createEmptyBorder());
	for (int i=0; i<genomeTextFields.size(); i++) innerRepliconPanel.add(generateRepliconPanel(i));
	innerRepliconPanel.add(getSpacerPanel(25, 20));
	innerRepliconPanel.add(genomeAddLabel);
	toggleRepliconDeNovoPanels(innerRepliconPanel, innerRepliconPanel);

	// Force scroll to bottom
	SwingUtilities.invokeLater( new Runnable() {
		public void run()
		{
		    JScrollBar v = outerRepliconPanel.getVerticalScrollBar();
		    v.setValue(v.getMaximum());
		}
	    } );

	validate();
	repaint();
    }

    /**
     * Update panel containing the tree of RNAseq file names.
     */
    private void updateSeqFilePanel() {
	mainPanel.add(exps.getScrollPane());
	mainPanel.add(spacer);
	mainPanel.add(submitLabel);
	validate();
	repaint();
    }

    /**
     * Action performed on genome text field.
     */
    private void handleGenomeTextField(int index) {
	if (genomeTextFields.get(index).getText().length() > 0) genomeTextFields.get(index).setBackground(white);
	removePopupMenu();
	if (exps.getNumExperiments() == 0) {
	    exps.addExperiment();
	    updateSeqFilePanel();
	}
    }

    private void addCustomGenome(int i) {
	// First show pop-up information message about custom replicon
	String title = "Use replicon whose sequence is not available from GenBank?";
	String message = "If the sequence of your replicon is avilable in GenBank, do not use this option, instead select the Cancel button \nand type the name of your replicon in the Replicon name field. \n\nIf the sequence of your replicion is NOT available in GenBank, select the OK button and indicate a folder that contains \na single FASTA file containing the replicon genomic sequence as well as a gene annotation file (in PTT format).";
	Object[] options = {"OK", "Cancel"};
	int n = JOptionPane.showOptionDialog(SwingUtilities.getAncestorOfClass(JFrame.class,this), message, title, JOptionPane.YES_NO_OPTION, JOptionPane.INFORMATION_MESSAGE, customGenomeButtonPressed, options, options[0]);

	// Select replicon folder
	if (n == 0) {
	    fc.setSelectedFile(null);
	    int returnVal = fc.showOpenDialog(SwingUtilities.getAncestorOfClass(JFrame.class,this));
	    if (returnVal == JFileChooser.APPROVE_OPTION) {
		try {
		    File dir = fc.getSelectedFile();
		    String genomeFileName = null;
		    String geneFileName = null;
		    String rnaFileName = null;
		    String[] files = dir.list();
		    for (String f : files) {
			if (f.toLowerCase().endsWith(".fna") || f.toLowerCase().endsWith(".fa") || f.toLowerCase().endsWith(".fasta")) genomeFileName = f;
			if (f.toLowerCase().endsWith(".ptt")) geneFileName = f;
			if (f.toLowerCase().endsWith(".rnt")) rnaFileName = f;
		    }
		    String genomeName = validateFastaFile(dir, genomeFileName);
		    if (genomeName != null) genomeFileName = getBaseFileName(genomeFileName) + ".fna";
		    if (genomeFileName == null) {
			error("Problem reading in FASTA file", "Could not identify in the specified folder a FASTA file with appropriate extension (fna, fa, or fasta) and in the expected format");
			return;
		    }
		    if (geneFileName == null) {
			error("Problem reading in gene annotation file", "Could not identify in the specified folder a gene annotation file with appropriate extension (ptt) and in the expected format");
			return;
		    }

		    genomeTextFields.get(i).setText(genomeName);  // Fill in custom replicon name in text field
		    customGenomeDirs.set(i, dir.getAbsolutePath());
		    handleGenomeTextField(i);
		    
		} catch (Exception e) {
		    error("Problem reading replicon files", "Unable to read replicon files in supplied folder");
		}
	    }
	}
    }

    /**
     * Switch between replicon panel and de novo panel.
     */
    private void toggleRepliconDeNovoPanels(JPanel p1, JPanel p2) {
	outerRepliconPanel.remove(p1);
	mainPanel.remove(outerRepliconPanel);
	mainPanel.remove(denovoPanel);
	outerRepliconPanel = new JScrollPane(p2);
	outerRepliconPanel.setPreferredSize(new Dimension(500, 100));
	outerRepliconPanel.setBackground(lightBlue);
	outerRepliconPanel.setBorder(BorderFactory.createEmptyBorder());	
	mainPanel.add(outerRepliconPanel);
	mainPanel.add(denovoPanel);
    }

    /**
     * Check FASTA file is in correct format. Rename the file's extension
     * to ".fna". If not in the correct format, return null. If in the
     * correct format, return the name of the genome.
     */
    private String validateFastaFile(File dir, String genomeFileName) {
	if (genomeFileName == null) return null;
	String baseFileName = getBaseFileName(genomeFileName);
	try {
	    File f = new File(dir, genomeFileName);
	    Scanner reader = new Scanner(f);
	    String line = reader.nextLine();
	    String name = baseFileName;
	    if (line.length() > 0) {
		if (line.charAt(0) == '>') {
		    name = line.substring(1);
		    if (name.split("\\|").length >= 5) {
			reader.close();
			boolean result = f.renameTo(new File(dir, baseFileName + ".fna"));
			String[] parse_line = name.split("\\|");
			if (result) return parse_line[4].trim();
			return null;
		    } else {  // Modify header line of FASTA file
			StringBuilder sequence = new StringBuilder();
			while (reader.hasNext()) sequence.append(reader.nextLine());
			reader.close();
			PrintWriter writer = new PrintWriter(new File(dir, baseFileName + ".fna"));
			writer.println(">gi|99999999|ref|NC_999999.9|" + name);
			writer.println(sequence.toString());
			writer.close();
			return name.trim();
		    }
		} else {  // Read in file
		    StringBuilder sequence = new StringBuilder(line);
		    while (reader.hasNext()) sequence.append(reader.nextLine());
		    reader.close();
		    PrintWriter writer = new PrintWriter(new File(dir, baseFileName + ".fna"));
		    writer.println(">gi|99999999|ref|NC_999999.9|" + baseFileName);
		    writer.println(sequence.toString());
		    writer.close();
		    return baseFileName;
		}
	    }
	    reader.close();
	} catch (FileNotFoundException e) {
	    return null;
	}
	return null;
    }

    private String getBaseFileName(String fileName) {
	if (fileName == null) return null;
	int slashIndex = fileName.lastIndexOf(File.separatorChar);
	int periodIndex = fileName.lastIndexOf('.');
	if (periodIndex == -1) return fileName.substring(slashIndex+1);
	else return fileName.substring(slashIndex+1, periodIndex); 
    }

    /**
     * Update the popup menu for the replicon name.
     */
    private void updateRepliconMenu(int z, int index) {
	String s = genomeTextFields.get(z).getText();
	if (s.length() >= NCBI.HASH_SIZE) {
	    popupMenu = new JPopupMenu("Replicons");
	    mainPanel.add(popupMenu);
	    currentReplicon = z;
	    for (IndexedMenuItem item : popupReplicons) item.removeActionListener(this);
	    popupReplicons.clear();
	    ArrayList<String> reps = ncbi.getRepliconNames(s);
	    for (int i=index; i<Math.min(reps.size(), index+25); i++) {
		IndexedMenuItem item = new IndexedMenuItem(reps.get(i));
		item.setBackground(white);
		popupReplicons.add(item);
		popupMenu.add(item);
		item.addActionListener(this);
	    }
	    if (reps.size() > index+25) {
		IndexedMenuItem item = new IndexedMenuItem("     *** MORE ***");
		item.setBackground(white);
		item.index = index+25;
		popupReplicons.add(item);
		popupMenu.add(item);
		item.addActionListener(this);
	    }
	    int popupHeight = Math.min(genomeTextFields.get(z).getHeight()+37+z*27, 95);
	    popupMenu.show(mainPanel, genomeTextFields.get(z).getX()+147, popupHeight);
	    genomeTextFields.get(z).requestFocusInWindow();
	} else removePopupMenu();
    }

    /**
     * Delete the popup menu.
     */
    private void removePopupMenu() {
	if (popupMenu != null) popupMenu.setVisible(false);
	popupMenu = null;
	currentReplicon = -1;
	for (IndexedMenuItem item : popupReplicons) item.removeActionListener(this);
	popupReplicons.clear();
	validate();
	repaint();
    }

    private boolean outputExperimentsToFile() {
	try {
	    File f = new File(directory);
	    if (!f.exists()) f.mkdir();
	    PrintWriter writer = new PrintWriter(new File(directory + fileName));
	    writer.println(Rockhopper.version);
	    writer.println("" + isDenovo);
	    if (!isDenovo) {  // Only write replicon info if we are performing reference based assembly
		writer.println(genomeTextFields.size());
		for (int i=0; i<genomeTextFields.size(); i++) {
		    if (customGenomeDirs.get(i) == null) writer.println(genomeTextFields.get(i).getText());
		    else writer.println(genomeTextFields.get(i).getText() + "\t" + customGenomeDirs.get(i));
		}
	    }
	    String s = exps.experimentsToString();
	    writer.print(s);
	    writer.close();
	    if (s.trim().length() == 0) {
		error("Invalid experiment files", "Please choose one or more valid experiment files");
		return false;
	    }
	} catch (FileNotFoundException ex) {
	    // Do nothing if experiment file could not be written.
	}
	return true;
    }

    private void readExperimentsFromFile() {
	try {
	    Scanner reader = new Scanner(new File(directory + fileName));
	    String version = reader.nextLine();
	    if (!version.equals(Rockhopper.version)) return;
	    isDenovo = Boolean.parseBoolean(reader.nextLine());
	    if (isDenovo) {  // De novo assembly
		denovoLabel.setIcon(denovoYesButton);
		toggleRepliconDeNovoPanels(innerRepliconPanel, denovoInfo);
	    } else {  // Reference based assembly
		denovoLabel.setIcon(denovoNoButton);
		int numReplicons = Integer.parseInt(reader.nextLine().trim());
		readingFromFile = true;
		if (numReplicons > 0) {
		    String repliconLine = reader.nextLine();
		    String[] parse_line = repliconLine.split("\t");
		    genomeTextFields.get(0).setText(parse_line[0]);
		    genomeTextFields.get(0).setBackground(white);
		    if (parse_line.length > 1) customGenomeDirs.set(0, parse_line[1]);
		}
		for (int i=1; i<numReplicons; i++) {
		    addReplicon();
		    String repliconLine = reader.nextLine();
		    String[] parse_line = repliconLine.split("\t");
		    genomeTextFields.get(genomeTextFields.size()-1).setText(parse_line[0]);
		    genomeTextFields.get(genomeTextFields.size()-1).setBackground(white);
		    if (parse_line.length > 1) customGenomeDirs.set(customGenomeDirs.size()-1, parse_line[1]);
		}
	    }
	    readingFromFile = false;
	    StringBuilder sb = new StringBuilder();
	    while (reader.hasNext()) sb.append(reader.nextLine() + "\n");
	    exps.experimentsFromString(sb.toString());
	    updateSeqFilePanel();
	    reader.close();
	} catch (Exception ex) {
	    // Do nothing if experiment file could not be read.
	}
    }

    private boolean prepareRepliconFiles() {
	if (genomeTextFields.size() == 0) {
	    error("No replicons selected", "Please designate one or more replicons");
	    return false;
	}
	for (int i=0; i<genomeTextFields.size(); i++) {
	    String repliconName = genomeTextFields.get(i).getText();
	    if ((customGenomeDirs.get(i) == null) && !ncbi.prepGenomeFiles(repliconName)) {
		error("Invalid replicon name", "Please enter a valid replicon name");
		return false;
	    }
	    if (customGenomeDirs.get(i) == null) {  // Use NCBI replicon
		if (i == 0) genome_DIR = NCBI.genome_DIR + NCBI.getRepliconDir(repliconName);
		else genome_DIR += "," + NCBI.genome_DIR + NCBI.getRepliconDir(repliconName);
	    } else {  // Use custom replicon
		if (i == 0) genome_DIR = customGenomeDirs.get(0);
		else genome_DIR += "," + customGenomeDirs.get(i);
	    }
	}
	return true;
    }

    /**
     * Returns a comma-delimited String of the WIG files generated
     * by Rockhopper that can be viewed in IGV.
     */
    private String getRockhopperWigFiles(Genome genome) {
	String s = "";
	String dir = Rockhopper.output_DIR + Rockhopper.browser_DIR;
	if (params.testDifferentialExp) s += "," + dir + genome.getID() + "_diffExpressedGenes.wig";
	if (params.transcriptBoundaries) s += "," + dir + genome.getID() + "_UTRs.wig";
	if (params.transcriptBoundaries) s += "," + dir + genome.getID() + "_ncRNAs.wig";
	if (params.predictOperons) s += "," + dir + genome.getID() + "_operons.wig";
	if (s.length() == 0) return s;
	else return s.substring(1);
    }

    /**
     * Given a comma-delimited list of WIG files representing
     * reads mapped to a genome by Peregrine, we clean/format
     * the comma-delimited list and return the 
     * cleaned/formatted version.
     */
    private String formatWigFiles(String wigFiles) {
	if (wigFiles.length() == 0) return "";
	ArrayList<String> wigs_plus = new ArrayList<String>();
	ArrayList<String> wigs_minus = new ArrayList<String>();
	String[] wigList = wigFiles.split(",");
	for (int i=0; i<wigList.length; i++) {
	    if (wigList[i].indexOf(".minus.") >= 0) wigs_minus.add(wigList[i]);
	    else wigs_plus.add(wigList[i]);
	}
	String s = "";
	for (int i=0; i<wigs_plus.size(); i++) {
	    if (s.length() > 0) s += ",";
	    s += wigs_plus.get(i);
	}
	for (int i=0; i<wigs_minus.size(); i++) {
	    if (s.length() > 0) s += ",";
	    s += wigs_minus.get(i);
	}
	return s;
    }

    private void clearCache() {
	removeAllFiles(new File(Rockhopper.output_DIR));
	File f = new File(Rockhopper.output_DIR);  // Create "Rockhopper_Results" directory
	f.mkdir();
	f = new File(directory);  // Create "intermediary" directory
	f.mkdir();
	mainPanel.removeAll();
	mainPanel.add(getSpacerPanel(790,20));
	ncbi = new NCBI();
	exps = new ExperimentsGUI(this);
	params = new ParametersGUI(this);
	progress = new ProgressGUI(this);
	removePopupMenu();
	genomeTextFields.clear();
	genomeLabels.clear();
	customGenomeLabels.clear();
	customGenomeDirs.clear();
	addReplicon();
	if (isDenovo) {  // Clearing cache in de novo assembly mode
	    toggleRepliconDeNovoPanels(innerRepliconPanel, denovoInfo);
	    exps.addExperiment();
	    updateSeqFilePanel();
	}
    }

    private void openUserGuide() {
	String address = "http://cs.wellesley.edu/~btjaden/Rockhopper/user_guide.html";
	if (Desktop.isDesktopSupported()) {
	    try {
		Desktop.getDesktop().browse(new URI(address));
	    } catch (Exception e) {
		JOptionPane.showMessageDialog(frame, "The user guide is available at the following URL:<BR><u>" + address + "</u>", "User Guide", JOptionPane.INFORMATION_MESSAGE);
	    }
	} else {
	    JOptionPane.showMessageDialog(frame, "The user guide is available at the following URL:<BR><u>" + address + "</u>", "User Guide", JOptionPane.INFORMATION_MESSAGE);
	}
    }

    /**
     * Recursively deletes all files in the specified directory.
     */
    private void removeAllFiles(File f) {
	if (!f.exists()) return;
	if (f.isDirectory()) {  // We have a directory. Look inside directory.
	    for (File x : f.listFiles()) removeAllFiles(x);
	}
	f.delete();  // Remove file or empty directory.
    }

    private String[] generateCommandLineArguments() {
	ArrayList<String> args = new ArrayList<String>();
	args.add("-m");
	args.add("" + params.mismatches);
	args.add("-l");
	args.add("" + params.minReadLength);
	args.add("-c");
	args.add("" + params.singleEndOrientation);
	args.add("-" + params.matePairOrientation);
	if (params.numProcessors > 0) {
	    args.add("-p");
	    args.add("" + params.numProcessors);
	}
	args.add("-d");
	args.add("" + params.maxPairedEndLength);
	args.add("-v");
	args.add("" + params.verboseOutput);
	args.add("-s");
	args.add("" + params.isStrandSpecific);
	args.add("-e");
	args.add("" + params.testDifferentialExp);
	args.add("-y");
	args.add("" + params.predictOperons);
	args.add("-t");
	args.add("" + params.transcriptBoundaries);
	args.add("-z");
	args.add("" + params.expressionParameter);
	args.add("-o");
	args.add("" + params.outputDIR);
	args.add("-L");
	args.add(exps.getExperimentsList());
	if (!isDenovo) {
	    args.add("-g");
	    args.add(genome_DIR);
	}
	args.add("-b");
	args.add("" + params.minReadsMapping);
	args.add("-u");
	args.add("" + params.minTranscriptLength);
	args.add("-w");
	args.add("" + params.minSeedExpression);
	args.add("-x");
	args.add("" + params.minExtendExpression);
	if (params.outputSAM) args.add("-SAM");
	String[] reps = exps.getReplicatesList();
	for (int i=0; i<reps.length; i++) args.add(reps[i]);
	return args.toArray(new String[args.size()]);
    }

    public static JPanel getSpacerPanel(int x, int y) {
        JPanel spacer = new JPanel(new FlowLayout(FlowLayout.CENTER));
        spacer.setPreferredSize(new Dimension(x, y));
        spacer.setBackground(RockhopperGUI.lightBlue);
        return spacer;
    }

    public static ImageIcon getScaledImage(String fileName, int x, int y) {
        return new ImageIcon(new ImageIcon(cl.getResource(fileName)).getImage().getScaledInstance(x,y,java.awt.Image.SCALE_SMOOTH));
    }

    private void error(String title, String message) {
        JOptionPane.showMessageDialog(frame, message, title, JOptionPane.ERROR_MESSAGE);
    }



    /*************************************
     **********   MAIN METHOD   **********
     *************************************/

    private static void createAndShowGUI() {    
	RockhopperGUI applet = new RockhopperGUI();
	//JFrame.setDefaultLookAndFeelDecorated(true);
	frame = new JFrame("Rockhopper");
	frame.setSize(800, 750);
	frame.setResizable(false);
	frame.setLocationRelativeTo(null);
	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	frame.setIconImage(new ImageIcon(penguin).getImage());
	applet.init();
	frame.add(applet);
	frame.setVisible(true);
    }

    public static void main(String[] args) {
	javax.swing.SwingUtilities.invokeLater(new Runnable() {
		public void run() { createAndShowGUI();}
	    });
    }

}



class IndexedMenuItem extends JMenuItem {
    public int index;
    public IndexedMenuItem(String s) { super(s); }
}

