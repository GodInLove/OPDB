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

import java.util.Hashtable;
import java.util.Scanner;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.MouseListener;
import javax.swing.*;
import java.io.PrintWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.NoSuchElementException;

public class ParametersGUI {

    private final String fileName = "parameters.txt";
    JApplet applet;

    // Default values
    private final double mismatches_DEFAULT = Peregrine.percentMismatches;
    private final double minReadLength_DEFAULT = Peregrine.percentSeedLength;
    private final boolean singleEndOrientation_DEFAULT = Peregrine.singleEndOrientationReverseComplement;
    private final String matePairOrientation_DEFAULT = Peregrine.pairedEndOrientation;
    private final int numProcessors_DEFAULT = 0;
    private final int maxPairedEndLength_DEFAULT = Peregrine.maxPairedEndLength;
    private final boolean verboseOutput_DEFAULT = Rockhopper.verbose;
    private final boolean outputSAM_DEFAULT = false;
    private final boolean isStrandSpecific_DEFAULT = !Rockhopper.unstranded;
    private final boolean testDifferentialExp_DEFAULT = Rockhopper.computeExpression;
    private final boolean predictOperons_DEFAULT = Rockhopper.computeOperons;
    private final boolean transcriptBoundaries_DEFAULT = Rockhopper.computeTranscripts;
    private final double expressionParameter_DEFAULT = Rockhopper.transcriptSensitivity;
    private final int minReadsMapping_DEFAULT = Assembler.MIN_READS_MAPPING;
    private final int minTranscriptLength_DEFAULT = 2 * Assembler.k;
    private final int minSeedExpression_DEFAULT = Assembler.minSeedExpression;
    private final int minExtendExpression_DEFAULT = Assembler.minExpression;
    private final String outputDIR_DEFAULT = Rockhopper.output_DIR;

    // Editable values
    public double mismatches;
    public double minReadLength;
    public boolean singleEndOrientation;
    public String matePairOrientation;
    public int numProcessors;
    public int maxPairedEndLength;
    public boolean verboseOutput;
    public boolean outputSAM;
    public boolean isStrandSpecific;
    public boolean testDifferentialExp;
    public boolean predictOperons;
    public boolean transcriptBoundaries;
    public double expressionParameter;
    public int minReadsMapping;
    public int minTranscriptLength;
    public int minSeedExpression;
    public int minExtendExpression;
    public String outputDIR;

    // GUI components
    private JFrame parameterFrame;
    private JTextField qualityTextField;
    private JTextField mismatchTextField;
    private JTextField readLengthTextField;
    private JCheckBox orientationCheckBox;
    private JComboBox orientationComboBox;
    private JTextField processorsTextField;
    private JTextField pairedEndLengthTextField;
    private JCheckBox verbose;
    private JCheckBox sam;
    private JCheckBox strand;
    private JCheckBox diffExpression;
    private JCheckBox operons;
    private JCheckBox boundaries;
    private JSlider minExpression;
    private JTextField minReadsMappingTextField;
    private JTextField minTranscriptLengthTextField;
    private JTextField minSeedExpressionTextField;
    private JTextField minExtendExpressionTextField;
    private JTextField dirTextField;
    private JLabel browseLabel;
    private JFileChooser fc;
    private File f;
    private ImageIcon browseButton, browseButtonPressed;
    private JButton defaultButton;
    private JButton saveButton;



    /**
     * Constructor.
     */
    public ParametersGUI(JApplet applet) {
	this.applet = applet;
	browseButton = RockhopperGUI.getScaledImage("images/browse.gif", 77, 29);
	browseButtonPressed = RockhopperGUI.getScaledImage("images/browsePressed.gif", 77, 29);

	// Initialize to default values
	mismatches = mismatches_DEFAULT;
	minReadLength = minReadLength_DEFAULT;
	singleEndOrientation = singleEndOrientation_DEFAULT;
	matePairOrientation = matePairOrientation_DEFAULT;
	numProcessors = numProcessors_DEFAULT;
	maxPairedEndLength = maxPairedEndLength_DEFAULT;
	verboseOutput = verboseOutput_DEFAULT;
	outputSAM = outputSAM_DEFAULT;
	isStrandSpecific = isStrandSpecific_DEFAULT;
	testDifferentialExp = testDifferentialExp_DEFAULT;
	predictOperons = predictOperons_DEFAULT;
	transcriptBoundaries = transcriptBoundaries_DEFAULT;
	expressionParameter = expressionParameter_DEFAULT;
	minReadsMapping = minReadsMapping_DEFAULT;
	minTranscriptLength = minTranscriptLength_DEFAULT;
	minSeedExpression = minSeedExpression_DEFAULT;
	minExtendExpression = minExtendExpression_DEFAULT;
	outputDIR = outputDIR_DEFAULT;

	// Read in from file previous parameter values.
	Scanner reader = null;
	fc = new JFileChooser();
	fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
	try {
	    File f = new File(RockhopperGUI.directory + fileName);
	    if (!f.exists()) return;
	    reader = new Scanner(f);
	} catch (FileNotFoundException ex) {
	    // Do nothing. Use default parameter values.
	} catch (NoSuchElementException ex) {
	    // Do nothing. Use default parameter values.
	}
	if (reader != null) {
	    try {
		String version = reader.nextLine();
		if (!version.equals(Rockhopper.version)) throw new Exception();
		double mismatches_TEMP = Double.parseDouble(reader.nextLine());
		double minReadLength_TEMP = Double.parseDouble(reader.nextLine());
		boolean singleEndOrientation_TEMP = Boolean.valueOf(reader.nextLine()).booleanValue();
		String matePairOrientation_TEMP = reader.nextLine();
		int numProcessors_TEMP = Integer.parseInt(reader.nextLine());
		int maxPairedEndLength_TEMP = Integer.parseInt(reader.nextLine());
		boolean verboseOutput_TEMP = Boolean.valueOf(reader.nextLine()).booleanValue();
		boolean outputSAM_TEMP = Boolean.valueOf(reader.nextLine()).booleanValue();
		boolean isStrandSpecific_TEMP = Boolean.valueOf(reader.nextLine()).booleanValue();
		boolean testDifferentialExp_TEMP = Boolean.valueOf(reader.nextLine()).booleanValue();
		boolean predictOperons_TEMP = Boolean.valueOf(reader.nextLine()).booleanValue();
		boolean transcriptBoundaries_TEMP = Boolean.valueOf(reader.nextLine()).booleanValue();
		double expressionParameter_TEMP = Double.parseDouble(reader.nextLine());
		int minReadsMapping_TEMP = Integer.parseInt(reader.nextLine());
		int minTranscriptLength_TEMP = Integer.parseInt(reader.nextLine());
		int minSeedExpression_TEMP = Integer.parseInt(reader.nextLine());
		int minExtendExpression_TEMP = Integer.parseInt(reader.nextLine());
		String outputDIR_TEMP = reader.nextLine();
		reader.close();

		mismatches = mismatches_TEMP;
		minReadLength = minReadLength_TEMP;
		singleEndOrientation = singleEndOrientation_TEMP;
		matePairOrientation = matePairOrientation_TEMP;
		numProcessors = numProcessors_TEMP;
		maxPairedEndLength = maxPairedEndLength_TEMP;
		verboseOutput = verboseOutput_TEMP;
		outputSAM = outputSAM_TEMP;
		isStrandSpecific = isStrandSpecific_TEMP;
		testDifferentialExp = testDifferentialExp_TEMP;
		predictOperons = predictOperons_TEMP;
		transcriptBoundaries = transcriptBoundaries_TEMP;
		expressionParameter = expressionParameter_TEMP;
		minReadsMapping = minReadsMapping_TEMP;
		minTranscriptLength = minTranscriptLength_TEMP;
		minSeedExpression = minSeedExpression_TEMP;
		minExtendExpression = minExtendExpression_TEMP;
		outputDIR = outputDIR_TEMP;
	    } catch (Exception e) { }
	}
    }

    public boolean isSource(Object source) {
	if (source.equals(defaultButton)) return true;
	if (source.equals(saveButton)) return true;
	if (source.equals(browseLabel)) return true;
	return false;
    }

    public void handleEvent(Object source) {
	if (source.equals(defaultButton)) {
	    mismatches = mismatches_DEFAULT;
	    minReadLength = minReadLength_DEFAULT;
	    singleEndOrientation = singleEndOrientation_DEFAULT;
	    matePairOrientation = matePairOrientation_DEFAULT;
	    numProcessors = numProcessors_DEFAULT;
	    maxPairedEndLength = maxPairedEndLength_DEFAULT;
	    verboseOutput = verboseOutput_DEFAULT;
	    outputSAM = outputSAM_DEFAULT;
	    isStrandSpecific = isStrandSpecific_DEFAULT;
	    testDifferentialExp = testDifferentialExp_DEFAULT;
	    predictOperons = predictOperons_DEFAULT;
	    transcriptBoundaries = transcriptBoundaries_DEFAULT;
	    expressionParameter = expressionParameter_DEFAULT;
	    minReadsMapping = minReadsMapping_DEFAULT;
	    minTranscriptLength = minTranscriptLength_DEFAULT;
	    minSeedExpression = minSeedExpression_DEFAULT;
	    minExtendExpression = minExtendExpression_DEFAULT;
	    outputDIR = outputDIR_DEFAULT;

	    mismatchTextField.setText("" + mismatches);
	    readLengthTextField.setText("" + minReadLength);
	    orientationCheckBox.setSelected(singleEndOrientation);
	    orientationComboBox.setSelectedItem(matePairOrientation);
	    processorsTextField.setText("" + numProcessors);
	    pairedEndLengthTextField.setText("" + maxPairedEndLength);
	    verbose.setSelected(verboseOutput);
	    sam.setSelected(outputSAM);
	    strand.setSelected(isStrandSpecific);
	    diffExpression.setSelected(testDifferentialExp);
	    operons.setSelected(predictOperons);
	    boundaries.setSelected(transcriptBoundaries);
	    minExpression.setValue((int)(expressionParameter*10));
	    minReadsMappingTextField.setText("" + minReadsMapping);
	    minTranscriptLengthTextField.setText("" + minTranscriptLength);
	    minSeedExpressionTextField.setText("" + minSeedExpression);
	    minExtendExpressionTextField.setText("" + minExtendExpression);
	    dirTextField.setText(outputDIR);
	    outputParametersToFile();
	} else if (source.equals(saveButton)) {
	    try { mismatches = Double.parseDouble(mismatchTextField.getText()); }
	    catch (NumberFormatException ex) {
		error("Parameter input error", "An invalid number was entered for the allowed mismatches.");
		return;
	    }
	    try { minReadLength = Double.parseDouble(readLengthTextField.getText()); }
	    catch (NumberFormatException ex) {
		error("Parameter input error", "An invalid number was entered for the minimum seed.");
		return;
	    }
	    try { numProcessors = Integer.parseInt(processorsTextField.getText()); }
	    catch (NumberFormatException ex) {
		error("Parameter input error", "An invalid integer was entered for the number of processors.");
		return;
	    }
	    try { maxPairedEndLength = Integer.parseInt(pairedEndLengthTextField.getText()); }
	    catch (NumberFormatException ex) {
		error("Parameter input error", "An invalid integer was entered for the maximum paired end length.");
		return;
	    }
	    try { minReadsMapping = Integer.parseInt(minReadsMappingTextField.getText()); }
	    catch (NumberFormatException ex) {
		error("Parameter input error", "An invalid integer was entered for the minimum number of mapping reads.");
		return;
	    }
	    try { minTranscriptLength = Integer.parseInt(minTranscriptLengthTextField.getText()); }
	    catch (NumberFormatException ex) {
		error("Parameter input error", "An invalid integer was entered for the minimum transcript length.");
		return;
	    }
	    try { minSeedExpression = Integer.parseInt(minSeedExpressionTextField.getText()); }
	    catch (NumberFormatException ex) {
		error("Parameter input error", "An invalid integer was entered for the minimum seed count.");
		return;
	    }
	    try { minExtendExpression = Integer.parseInt(minExtendExpressionTextField.getText()); }
	    catch (NumberFormatException ex) {
		error("Parameter input error", "An invalid integer was entered for the minimum extension count.");
		return;
	    }
	    singleEndOrientation = orientationCheckBox.isSelected();
	    matePairOrientation = (String)orientationComboBox.getSelectedItem();
	    verboseOutput = verbose.isSelected();
	    outputSAM = sam.isSelected();
	    isStrandSpecific = strand.isSelected();
	    testDifferentialExp = diffExpression.isSelected();
	    predictOperons = operons.isSelected();
	    transcriptBoundaries = boundaries.isSelected();
	    expressionParameter = minExpression.getValue()/10.0;
	    outputDIR = dirTextField.getText();
	    outputParametersToFile();
	} else if (source.equals(browseLabel)) {
	    fileChooser();
	}
    }

    public void mousePressed(Object source) {
	if (source.equals(browseLabel)) {
	    browseLabel.setIcon(browseButtonPressed);
	}
    }

    public void mouseReleased(Object source) {
	if (source.equals(browseLabel)) {
	    browseLabel.setIcon(browseButton);
	}
    }

    /** 
     * Open frame to display parameter options.
     */
    public void displayParameters() {
	int width = 700;
	parameterFrame = new JFrame("Parameter Settings");
	parameterFrame.setSize(width, 620);
	parameterFrame.setResizable(false);
	parameterFrame.setLocationRelativeTo(null);
	parameterFrame.setIconImage(new ImageIcon(RockhopperGUI.penguin).getImage());
	JPanel p1 = new JPanel(new FlowLayout(FlowLayout.CENTER));
	p1.setPreferredSize(new Dimension(width, 620));
	p1.setBackground(RockhopperGUI.lightBlue);
	JPanel spacerPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
	spacerPanel.setPreferredSize(new Dimension(width, 5));
	spacerPanel.setBackground(RockhopperGUI.lightBlue);
	p1.add(spacerPanel);

	JLabel generalLabel = new JLabel("General parameters");
	generalLabel.setForeground(RockhopperGUI.gold);
	generalLabel.setFont(new Font("Times", Font.ITALIC, 18));
	p1.add(generalLabel);
	JPanel gen = new JPanel(new GridLayout(4,2));
	gen.setPreferredSize(new Dimension(width-150, 120));
	gen.setBackground(RockhopperGUI.lightBlue);

	// Strand specific
	JPanel gen_a = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	gen_a.setBackground(RockhopperGUI.lightBlue);
	JLabel strandLabel = new JLabel("Strand specific");
	strandLabel.setToolTipText("Indicates whether RNA-seq experiments are strand specific (checked) or strand ambiguous (unchecked)");
	gen_a.add(strandLabel);
	strand = new JCheckBox();
	strand.setBackground(RockhopperGUI.lightBlue);
	strand.setSelected(isStrandSpecific);
	strand.setToolTipText("Indicates whether RNA-seq experiments are strand specific (checked) or strand ambiguous (unchecked)");
	gen_a.add(strand);
	gen.add(gen_a);

	// Test for differential expression
	JPanel gen_b = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	gen_b.setBackground(RockhopperGUI.lightBlue);
	JLabel diffExpressionLabel = new JLabel("Test for differential expression");
	diffExpressionLabel.setToolTipText("Test for differential expression of each transcripts in pairs of experiments");
	gen_b.add(diffExpressionLabel);
	diffExpression = new JCheckBox();
	diffExpression.setBackground(RockhopperGUI.lightBlue);
	diffExpression.setSelected(testDifferentialExp);
	diffExpression.setToolTipText("Test for differential expression of each transcripts in pairs of experiments");
	gen_b.add(diffExpression);
	gen.add(gen_b);

	// Single-end orientation
	JPanel gen_c = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	gen_c.setBackground(RockhopperGUI.lightBlue);
	JLabel singleEndReverseComplementLabel = new JLabel("Reverse complement reads");
	singleEndReverseComplementLabel.setToolTipText("Reverse complement all single-end reads");
	gen_c.add(singleEndReverseComplementLabel);
	orientationCheckBox = new JCheckBox();
	orientationCheckBox.setBackground(RockhopperGUI.lightBlue);
	orientationCheckBox.setSelected(singleEndOrientation);
	orientationCheckBox.setToolTipText("Reverse complement all single-end reads");
	gen_c.add(orientationCheckBox);
	gen.add(gen_c);

	// Paired-end orientation
	JPanel gen_d = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	gen_d.setBackground(RockhopperGUI.lightBlue);
	JLabel orientationLabel = new JLabel("Orientation of mate-pair reads");
	orientationLabel.setToolTipText("Orientation of two mate reads for paired-end reads, f=forward and r=reverse_complement");
	gen_d.add(orientationLabel);
	orientationComboBox = new JComboBox();
	orientationComboBox.setToolTipText("Orientation of two mate reads for paired-end reads, f=forward and r=reverse_complement");
	orientationComboBox.addItem("ff");
	orientationComboBox.addItem("fr");
	orientationComboBox.addItem("rf");
	orientationComboBox.addItem("rr");
	orientationComboBox.setSelectedItem(matePairOrientation);
	gen_d.add(orientationComboBox);
	gen.add(gen_d);

	// Verbose
	JPanel gen_e = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	gen_e.setBackground(RockhopperGUI.lightBlue);
	JLabel verboseLabel = new JLabel("Verbose output");
	verboseLabel.setToolTipText("Verbose output to transcript.txt file, including raw counts, normalized counts, and expression variances");
	gen_e.add(verboseLabel);
	verbose = new JCheckBox();
	verbose.setBackground(RockhopperGUI.lightBlue);
	verbose.setSelected(verboseOutput);
	verbose.setToolTipText("Verbose output to transcripts.txt file, including raw counts, normalized counts, and expression variances");
	gen_e.add(verbose);
	gen.add(gen_e);

	// Output SAM
	JPanel gen_f = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	gen_f.setBackground(RockhopperGUI.lightBlue);
	JLabel samLabel = new JLabel("Output SAM file");
	samLabel.setToolTipText("Output sequencing reads to a SAM file");
	gen_f.add(samLabel);
	sam = new JCheckBox();
	sam.setBackground(RockhopperGUI.lightBlue);
	sam.setSelected(outputSAM);
	sam.setToolTipText("Output sequencing reads to a SAM file");
	gen_f.add(sam);
	gen.add(gen_f);

	// Maximum paired-end length
	JPanel gen_g = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	gen_g.setBackground(RockhopperGUI.lightBlue);
	JLabel pairedEndLengthLabel = new JLabel("Max bases between paired reads");
	pairedEndLengthLabel.setToolTipText("Maximum number of bases between mate pairs for paired-end reads");
	gen_g.add(pairedEndLengthLabel);
	pairedEndLengthTextField = new JTextField("" + maxPairedEndLength, 4);
	pairedEndLengthTextField.setToolTipText("Maximum number of bases between mate pairs for paired-end reads");
	gen_g.add(pairedEndLengthTextField);
	gen.add(gen_g);

	// Number of processors
	JPanel gen_h = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	gen_h.setBackground(RockhopperGUI.lightBlue);
	JLabel processorsLabel = new JLabel("Number of processors");
	processorsLabel.setToolTipText("Number of parallel processors to use; a value of zero indicates that Rockhopper should determine the number (recommended)");
	gen_h.add(processorsLabel);
	processorsTextField = new JTextField("" + numProcessors, 3);
	processorsTextField.setToolTipText("Number of parallel processors to use; a value of zero indicates that Rockhopper should determine the number (recommended)");
	gen_h.add(processorsTextField);
	gen.add(gen_h);
	p1.add(gen);

	// Output file location
	JPanel output = new JPanel(new FlowLayout(FlowLayout.CENTER));
	output.setPreferredSize(new Dimension(width-50, 40));
	output.setBackground(RockhopperGUI.lightBlue);
	JLabel dirLabel = new JLabel("Output file location:");
	dirLabel.setBackground(RockhopperGUI.lightBlue);
	dirLabel.setToolTipText("Location where output files should be written");
	dirTextField = new JTextField(outputDIR, 25);
	dirTextField.setToolTipText("Location where output files should be written");
	browseLabel = new JLabel(browseButton);
	browseLabel.setToolTipText("Search for location");
	browseLabel.addMouseListener((MouseListener)applet);
	output.add(dirLabel);
	output.add(dirTextField);
	output.add(browseLabel);

	// Spacer
	p1.add(output);
	spacerPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
	spacerPanel.setPreferredSize(new Dimension(width, 10));
	spacerPanel.setBackground(RockhopperGUI.lightBlue);
	p1.add(spacerPanel);

	// Parameters for reference based assembly
	JLabel referenceBasedLabel = new JLabel("Parameters for reference based assembly");
	referenceBasedLabel.setForeground(RockhopperGUI.gold);
	referenceBasedLabel.setFont(new Font("Times", Font.ITALIC, 18));
	p1.add(referenceBasedLabel);
	JPanel ref = new JPanel(new GridLayout(2,2));
	ref.setPreferredSize(new Dimension(width-150, 70));
	ref.setBackground(RockhopperGUI.lightBlue);

	// Allowed mismatches
	JPanel ref_a = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	ref_a.setBackground(RockhopperGUI.lightBlue);
	JLabel mismatchLabel = new JLabel("Allowed mismatches");
	mismatchLabel.setToolTipText("Maximum number of mismatches/gaps/insertions allowed when aligning read to genome (as percentage of read length)");
	ref_a.add(mismatchLabel);
	mismatchTextField = new JTextField("" + mismatches, 4);
	mismatchTextField.setToolTipText("Maximum number of mismatches/gaps/insertions allowed when aligning read to genome (as percentage of read length)");
	ref_a.add(mismatchTextField);
	ref.add(ref_a);

	// Minimum seed length
	JPanel ref_b = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	ref_b.setBackground(RockhopperGUI.lightBlue);
	JLabel readLengthLabel = new JLabel("Minimum seed length");
	readLengthLabel.setToolTipText("Minimum number of consecutive nucleotides aligning to genome exactly (as percentage of read length)");
	ref_b.add(readLengthLabel);
	readLengthTextField = new JTextField("" + minReadLength, 4);
	readLengthTextField.setToolTipText("Minimum number of consecutive nucleotides aligning to genome exactly (as percentage of read length)");
	ref_b.add(readLengthTextField);
	ref.add(ref_b);

	// Transcript boundaries
	JPanel ref_c = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	ref_c.setBackground(RockhopperGUI.lightBlue);
	JLabel boundariesLabel = new JLabel("Identify transcript boundaries");
	boundariesLabel.setToolTipText("Determine UTRs and novel noncoding RNAs");
	ref_c.add(boundariesLabel);
	boundaries = new JCheckBox();
	boundaries.setBackground(RockhopperGUI.lightBlue);
	boundaries.setSelected(transcriptBoundaries);
	boundaries.setToolTipText("Determine UTRs and novel noncoding RNAs");
	ref_c.add(boundaries);
	ref.add(ref_c);

	// Operons
	JPanel ref_d = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	ref_d.setBackground(RockhopperGUI.lightBlue);
	JLabel operonsLabel = new JLabel("Predict operons");
	operonsLabel.setToolTipText("Identify polycistronic messages");
	ref_d.add(operonsLabel);
	operons = new JCheckBox();
	operons.setBackground(RockhopperGUI.lightBlue);
	operons.setSelected(predictOperons);
	operons.setToolTipText("Identify polycistronic messages");
	ref_d.add(operons);
	ref.add(ref_d);

	// Panel for min expression of UTRs and ncRNAs
	p1.add(ref);
	JPanel ref_ef = new JPanel(new FlowLayout(FlowLayout.CENTER));
	//ref_ef.setPreferredSize(new Dimension(width-50, 75));
	ref_ef.setBackground(RockhopperGUI.lightBlue);
	JLabel minExpressionLabel = new JLabel("Minimum expression of UTRs and ncRNAs");
	minExpressionLabel.setToolTipText("Parameter indicating minimum expression level used to determine UTRs and novel noncoding RNAs");
	ref_ef.add(minExpressionLabel);
	minExpression = new JSlider(JSlider.HORIZONTAL, 0, 10, (int)(expressionParameter*10));
	minExpression.setBackground(RockhopperGUI.lightBlue);
	minExpression.setForeground(Color.BLACK);
	minExpression.setMajorTickSpacing(5);
	minExpression.setMinorTickSpacing(1);
	minExpression.setPaintTicks(true);
	Hashtable<Integer, JLabel> labelTable = new Hashtable<Integer, JLabel>();
	labelTable.put(new Integer(0), new JLabel("0.0"));
	labelTable.put(new Integer(5), new JLabel("0.5"));
	labelTable.put(new Integer(10), new JLabel("1.0"));
	minExpression.setLabelTable(labelTable);
	minExpression.setPaintLabels(true);
	minExpression.setToolTipText("Parameter indicating minimum expression level used to determine UTRs and novel noncoding RNAs");
	ref_ef.add(minExpression);

	// Spacer
	p1.add(ref_ef);
	spacerPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
	spacerPanel.setPreferredSize(new Dimension(width, 10));
	spacerPanel.setBackground(RockhopperGUI.lightBlue);
	p1.add(spacerPanel);

	// Parameters for de novo assembly
	JLabel denovoLabel = new JLabel("Parameters for de novo assembly");
	denovoLabel.setForeground(RockhopperGUI.gold);
	denovoLabel.setFont(new Font("Times", Font.ITALIC, 18));
	p1.add(denovoLabel);
	JPanel denovo = new JPanel(new GridLayout(2,2));
	denovo.setPreferredSize(new Dimension(width-150, 70));
	denovo.setBackground(RockhopperGUI.lightBlue);

	// Min reads mapping
	JPanel denovo_a = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	denovo_a.setBackground(RockhopperGUI.lightBlue);
	JLabel minReadsMappingLabel = new JLabel("Min reads mapping to a transcript");
	minReadsMappingLabel.setToolTipText("Minimum number of full length reads required to map to a de novo assembled transcript");
	denovo_a.add(minReadsMappingLabel);
	minReadsMappingTextField = new JTextField("" + minReadsMapping, 3);
	minReadsMappingTextField.setToolTipText("Minimum number of full length reads required to map to a de novo assembled transcript");
	denovo_a.add(minReadsMappingTextField);
	denovo.add(denovo_a);

	// Min transcript length
	JPanel denovo_b = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	denovo_b.setBackground(RockhopperGUI.lightBlue);
	JLabel minTranscriptLengthLabel = new JLabel("Minimum transcript length");
	minTranscriptLengthLabel.setToolTipText("Minimum length of de novo assembled transcripts (in nts)");
	denovo_b.add(minTranscriptLengthLabel);
	minTranscriptLengthTextField = new JTextField("" + minTranscriptLength, 3);
	minTranscriptLengthTextField.setToolTipText("Minimum length of de novo assembled transcripts (in nts)");
	denovo_b.add(minTranscriptLengthTextField);
	denovo.add(denovo_b);

	// Min seed expression
	JPanel denovo_c = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	denovo_c.setBackground(RockhopperGUI.lightBlue);
	JLabel minSeedExpressionLabel = new JLabel("Min count to seed a transcript");
	minSeedExpressionLabel.setToolTipText("Minimum count to seed a new de novo assembled transcript");
	denovo_c.add(minSeedExpressionLabel);
	minSeedExpressionTextField = new JTextField("" + minSeedExpression, 3);
	minSeedExpressionTextField.setToolTipText("Minimum count to seed a new de novo assembled transcript");
	denovo_c.add(minSeedExpressionTextField);
	denovo.add(denovo_c);

	// Min extend expression
	JPanel denovo_d = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	denovo_d.setBackground(RockhopperGUI.lightBlue);
	JLabel minExtendExpressionLabel = new JLabel("Min count to extend a transcript");
	minExtendExpressionLabel.setToolTipText("Minimum count to extend an existing de novo assembled transcript");
	denovo_d.add(minExtendExpressionLabel);
	minExtendExpressionTextField = new JTextField("" + minExtendExpression, 3);
	minExtendExpressionTextField.setToolTipText("Minimum count to extend an existing de novo assembled transcript");
	denovo_d.add(minExtendExpressionTextField);
	denovo.add(denovo_d);

	// Spacer
	p1.add(denovo);
	spacerPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
	spacerPanel.setPreferredSize(new Dimension(width, 25));
	spacerPanel.setBackground(RockhopperGUI.lightBlue);
	p1.add(spacerPanel);

	// Defaults and Save buttons
	JPanel buttons = new JPanel(new GridLayout(1,2,30,0));
	buttons.setPreferredSize(new Dimension(width-400, 25));
	buttons.setBackground(RockhopperGUI.lightBlue);
	defaultButton = new JButton("DEFAULTS");
	defaultButton.setToolTipText("Restore default parameters values");
	defaultButton.setForeground(Color.red);
	defaultButton.setBackground(RockhopperGUI.gold);
	defaultButton.addActionListener((ActionListener)applet);
	saveButton = new JButton("SAVE");
	saveButton.setToolTipText("Save parameter values");
	saveButton.setForeground(RockhopperGUI.gold);
	saveButton.setBackground(Color.red);
	saveButton.addActionListener((ActionListener)applet);
	buttons.add(defaultButton);
	buttons.add(saveButton);
	p1.add(buttons);
	parameterFrame.add(p1);
	parameterFrame.setVisible(true);
	saveButton.requestFocus();
    }

    private void fileChooser() {
	fc.setSelectedFile(null);
	int returnVal = fc.showOpenDialog(SwingUtilities.getAncestorOfClass(JFrame.class,applet));
	f = fc.getSelectedFile();
	if (f != null) dirTextField.setText(f.getAbsolutePath());
    }

    private void outputParametersToFile() {
	try {
	    File f = new File(RockhopperGUI.directory);
	    if (!f.exists()) f.mkdir();
	    PrintWriter writer = new PrintWriter(new File(RockhopperGUI.directory + fileName));
	    writer.println(Rockhopper.version);
	    writer.println(mismatches);
	    writer.println(minReadLength);
	    writer.println(singleEndOrientation);
	    writer.println(matePairOrientation);
	    writer.println(numProcessors);
	    writer.println(maxPairedEndLength);
	    writer.println(verboseOutput);
	    writer.println(outputSAM);
	    writer.println(isStrandSpecific);
	    writer.println(testDifferentialExp);
	    writer.println(predictOperons);
	    writer.println(transcriptBoundaries);
	    writer.println(expressionParameter);
	    writer.println(minReadsMapping);
	    writer.println(minTranscriptLength);
	    writer.println(minSeedExpression);
	    writer.println(minExtendExpression);
	    writer.println(outputDIR);
	    writer.close();
	} catch (FileNotFoundException ex) {
	    // Do nothing if parameter file could not be written.
	}
    }

    private void error(String title, String message) {
        JOptionPane.showMessageDialog(parameterFrame, message, title, JOptionPane.ERROR_MESSAGE);
    }

}

