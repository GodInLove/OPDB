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
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.MouseListener;
import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.io.IOException;
import java.io.File;

public class ExperimentsGUI {

    private JApplet applet;
    private JFileChooser fc;
    private int numExperiments;
    private int currentExperimentNumber;
    private JPanel panel;
    private JScrollPane scrollPane;
    private ArrayList<ArrayList<JTextField>> fileNames;
    private ArrayList<ArrayList<File>> files;
    private ArrayList<ArrayList<File>> pairedEndFiles;
    private ArrayList<ArrayList<JLabel>> replicateBrowseLabels;
    private ArrayList<JLabel> replicateAddLabels;
    private ArrayList<ArrayList<JLabel>> replicateRemoveLabels;
    private ArrayList<ArrayList<JLabel>> pairedEndLabels;
    private ArrayList<JTextField> experimentNames;
    private JLabel experimentAddLabel;
    private ImageIcon experimentAddButton, experimentAddButtonPressed, replicateAddButton, replicateAddButtonPressed, experimentRemoveButton, experimentRemoveButtonPressed, replicateRemoveButton, replicateRemoveButtonPressed, browseButton, browseButtonPressed, pairedEndNoButton, pairedEndYesButton, pairedEndNoButtonPressed, pairedEndYesButtonPressed;
    private ArrayList<JLabel> experimentRemoveLabels;
    private ArrayList<JPanel> experimentPanels;

    public ExperimentsGUI(JApplet applet) {
	this.applet = applet;
	fc = new JFileChooser();
	numExperiments = 0;
	currentExperimentNumber = 0;
	panel = new JPanel(new FlowLayout(FlowLayout.CENTER,5,5));
	scrollPane = new JScrollPane(panel);
	scrollPane.setPreferredSize(new Dimension(720, 500));
	scrollPane.setBorder(BorderFactory.createEmptyBorder());
	fileNames = new ArrayList<ArrayList<JTextField>>();
	files = new ArrayList<ArrayList<File>>();
	pairedEndFiles = new ArrayList<ArrayList<File>>();
	replicateBrowseLabels = new ArrayList<ArrayList<JLabel>>();
	replicateAddLabels = new ArrayList<JLabel>();
	replicateRemoveLabels = new ArrayList<ArrayList<JLabel>>();
	pairedEndLabels = new ArrayList<ArrayList<JLabel>>();
	experimentNames = new ArrayList<JTextField>();
	experimentAddButton = RockhopperGUI.getScaledImage("images/addExp.gif", 20, 20);
	experimentAddButtonPressed = RockhopperGUI.getScaledImage("images/addExpPressed.gif", 20, 20);
	replicateAddButton = RockhopperGUI.getScaledImage("images/addRep.gif", 20, 20);
	replicateAddButtonPressed = RockhopperGUI.getScaledImage("images/addRepPressed.gif", 20, 20);
	experimentRemoveButton = RockhopperGUI.getScaledImage("images/removeExp.gif", 20, 20);
	experimentRemoveButtonPressed = RockhopperGUI.getScaledImage("images/removeExpPressed.gif", 20, 20);
	replicateRemoveButton = RockhopperGUI.getScaledImage("images/removeRep.gif", 20, 20);
	replicateRemoveButtonPressed = RockhopperGUI.getScaledImage("images/removeRepPressed.gif", 20, 20);
	browseButton = RockhopperGUI.getScaledImage("images/browse.gif", 77, 29);
	browseButtonPressed = RockhopperGUI.getScaledImage("images/browsePressed.gif", 77, 29);
	pairedEndNoButton = RockhopperGUI.getScaledImage("images/pairedEndRed.gif", 20, 20);
	pairedEndYesButton = RockhopperGUI.getScaledImage("images/pairedEndGreen.gif", 20, 20);
	pairedEndNoButtonPressed = RockhopperGUI.getScaledImage("images/pairedEndRedPressed.gif", 20, 20);
	pairedEndYesButtonPressed = RockhopperGUI.getScaledImage("images/pairedEndGreenPressed.gif", 20, 20);
	experimentAddLabel = new JLabel(experimentAddButton);
	experimentAddLabel.setToolTipText("Add experiment");
	experimentAddLabel.addMouseListener((MouseListener)applet);
	experimentRemoveLabels = new ArrayList<JLabel>();
	experimentPanels = new ArrayList<JPanel>();
    }

    public int getNumExperiments() {
	return numExperiments;
    }

    public JScrollPane getScrollPane() {
	return scrollPane;
    }

    public boolean isSource(Object source) {
	if (source.equals(experimentAddLabel)) return true;
	for (int i=0; i<numExperiments; i++) {
	    if (source.equals(experimentRemoveLabels.get(i))) return true;
	    if (source.equals(replicateAddLabels.get(i))) return true;
	    for (int j=0; j<fileNames.get(i).size(); j++) {
		if (source.equals(replicateBrowseLabels.get(i).get(j))) return true;
		if (source.equals(replicateRemoveLabels.get(i).get(j))) return true;
		if (source.equals(pairedEndLabels.get(i).get(j))) return true;
	    }
	}
	return false;
    }

    public void handleEvent(Object source) {
	if (source.equals(experimentAddLabel)) {
	    addExperiment();
	    return;
	}
	for (int i=0; i<numExperiments; i++) {
	    if (source.equals(experimentRemoveLabels.get(i))) {
		removeExperiment(i);
		return;
	    }
	    if (source.equals(replicateAddLabels.get(i))) {
		addReplicate(i);
		return;
	    }
	    for (int j=0; j<fileNames.get(i).size(); j++) {
		if (source.equals(replicateBrowseLabels.get(i).get(j))) {
		    fileChooser(i,j);
		    return;
		}
		if (source.equals(replicateRemoveLabels.get(i).get(j))) {
		    removeReplicate(i,j);
		    return;
		}
		if (source.equals(pairedEndLabels.get(i).get(j))) {
		    addMatePairFile(i, j);
		}
	    }
	}
    }

    public void mousePressed(Object source) {
	if (source.equals(experimentAddLabel)) {
	    experimentAddLabel.setIcon(experimentAddButtonPressed);
	    return;
	}
	for (int i=0; i<numExperiments; i++) {
	    if (source.equals(experimentRemoveLabels.get(i))) {
		experimentRemoveLabels.get(i).setIcon(experimentRemoveButtonPressed);
		return;
	    }
	    if (source.equals(replicateAddLabels.get(i))) {
		replicateAddLabels.get(i).setIcon(replicateAddButtonPressed);
		return;
	    }
	    for (int j=0; j<fileNames.get(i).size(); j++) {
		if (source.equals(replicateBrowseLabels.get(i).get(j))) {
		    replicateBrowseLabels.get(i).get(j).setIcon(browseButtonPressed);
		    return;
		}
		if (source.equals(replicateRemoveLabels.get(i).get(j))) {
		    replicateRemoveLabels.get(i).get(j).setIcon(replicateRemoveButtonPressed);
		    return;
		}
		if (source.equals(pairedEndLabels.get(i).get(j))) {
		    if (pairedEndLabels.get(i).get(j).getIcon().equals(pairedEndNoButton))
			pairedEndLabels.get(i).get(j).setIcon(pairedEndNoButtonPressed);
		    if (pairedEndLabels.get(i).get(j).getIcon().equals(pairedEndYesButton))
			pairedEndLabels.get(i).get(j).setIcon(pairedEndYesButtonPressed);
		    return;
		}
	    }
	}
    }

    public void mouseReleased(Object source) {
	if (source.equals(experimentAddLabel)) {
	    experimentAddLabel.setIcon(experimentAddButton);
	    return;
	}
	for (int i=0; i<numExperiments; i++) {
	    if (source.equals(experimentRemoveLabels.get(i))) {
		experimentRemoveLabels.get(i).setIcon(experimentRemoveButton);
		return;
	    }
	    if (source.equals(replicateAddLabels.get(i))) {
		replicateAddLabels.get(i).setIcon(replicateAddButton);
		return;
	    }
	    for (int j=0; j<fileNames.get(i).size(); j++) {
		if (source.equals(replicateBrowseLabels.get(i).get(j))) {
		    replicateBrowseLabels.get(i).get(j).setIcon(browseButton);
		    return;
		}
		if (source.equals(replicateRemoveLabels.get(i).get(j))) {
		    replicateRemoveLabels.get(i).get(j).setIcon(replicateRemoveButton);
		    return;
		}
		if (source.equals(pairedEndLabels.get(i).get(j))) {
		    if (pairedEndLabels.get(i).get(j).getIcon().equals(pairedEndNoButtonPressed))
			pairedEndLabels.get(i).get(j).setIcon(pairedEndNoButton);
		    if (pairedEndLabels.get(i).get(j).getIcon().equals(pairedEndYesButtonPressed))
			pairedEndLabels.get(i).get(j).setIcon(pairedEndYesButton);
		    return;
		}
	    }
	}
    }

    public void addExperiment() {
	scrollPane.remove(panel);
	numExperiments++;
	currentExperimentNumber++;
	fileNames.add(new ArrayList<JTextField>());
	fileNames.get(fileNames.size()-1).add(new JTextField(30));
	fileNames.get(fileNames.size()-1).get(fileNames.get(fileNames.size()-1).size()-1).setBackground(RockhopperGUI.gold);
	files.add(new ArrayList<File>());
	files.get(files.size()-1).add(null);
	pairedEndFiles.add(new ArrayList<File>());
	pairedEndFiles.get(pairedEndFiles.size()-1).add(null);
	replicateBrowseLabels.add(new ArrayList<JLabel>());
	JLabel repBrowse = new JLabel(browseButton);
	repBrowse.setToolTipText("Search for file");
	repBrowse.addMouseListener((MouseListener)applet);
	replicateBrowseLabels.get(replicateBrowseLabels.size()-1).add(repBrowse);
	JLabel repAdd = new JLabel(replicateAddButton);
	repAdd.setToolTipText("Add replicate");
	repAdd.addMouseListener((MouseListener)applet);
	replicateAddLabels.add(repAdd);
	replicateRemoveLabels.add(new ArrayList<JLabel>());
	JLabel repRemove = new JLabel(replicateRemoveButton);
	repRemove.setToolTipText("Remove replicate");
	repRemove.addMouseListener((MouseListener)applet);
	replicateRemoveLabels.get(replicateRemoveLabels.size()-1).add(repRemove);
	pairedEndLabels.add(new ArrayList<JLabel>());
	JLabel pairedEndAdd = new JLabel(pairedEndNoButton);
	pairedEndAdd.setToolTipText("Using paired-end reads; add file of mate-pair reads");
	pairedEndAdd.addMouseListener((MouseListener)applet);
	pairedEndLabels.get(pairedEndLabels.size()-1).add(pairedEndAdd);
	JTextField expName = new JTextField("" + currentExperimentNumber, 15);
	expName.addActionListener((ActionListener)applet);
	experimentNames.add(expName);
	JLabel expRemove = new JLabel(experimentRemoveButton);
	expRemove.setToolTipText("Remove experiment");
	expRemove.addMouseListener((MouseListener)applet);
	experimentRemoveLabels.add(expRemove);
	experimentPanels.add(generateExperimentPanel(experimentNames.size()-1));
	populateScrollPane();
    }

    public void removeExperiment(int i) {
	scrollPane.remove(panel);
	numExperiments--;
	fileNames.remove(i);
	files.remove(i);
	pairedEndFiles.remove(i);
	replicateBrowseLabels.remove(i);
	replicateAddLabels.remove(i);
	replicateRemoveLabels.remove(i);
	pairedEndLabels.remove(i);
	experimentNames.remove(i);
	experimentRemoveLabels.remove(i);
	experimentPanels.remove(i);
	populateScrollPane();
    }

    public void fileChooser(int i, int j) {
	fc.setSelectedFile(null);
	int returnVal = fc.showOpenDialog(SwingUtilities.getAncestorOfClass(JFrame.class,applet));
	if (returnVal == JFileChooser.APPROVE_OPTION) {
	    files.get(i).set(j, fc.getSelectedFile());
	    fileNames.get(i).get(j).setText(files.get(i).get(j).getName());
	    fileNames.get(i).get(j).setBackground(Color.WHITE);
	}
    }

    public void addMatePairFile(int i, int j) {
	// First show pop-up information message about paired-end reads
	String title = "Use mate-pair file for paired-end reads?";
	String message = "If you are using paired-end reads, use this option if you have two mate-pair files for a particular replicate experiment.";
	Object[] options = {"OK", "Cancel"};
	int n = JOptionPane.showOptionDialog(SwingUtilities.getAncestorOfClass(JFrame.class,applet), message, title, JOptionPane.YES_NO_OPTION, JOptionPane.INFORMATION_MESSAGE, pairedEndYesButton, options, options[0]);

	// Select mate-pair file
	if (n == 0) {
	    fc.setSelectedFile(null);
	    int returnVal = fc.showOpenDialog(SwingUtilities.getAncestorOfClass(JFrame.class,applet));
	    if (returnVal == JFileChooser.APPROVE_OPTION) {
		pairedEndFiles.get(i).set(j, fc.getSelectedFile());
		pairedEndLabels.get(i).get(j).setIcon(pairedEndYesButton);
	    }
	}
    }

    public void addReplicate(int i) {
	scrollPane.remove(panel);
	fileNames.get(i).add(new JTextField(30));
	fileNames.get(i).get(fileNames.get(i).size()-1).setBackground(RockhopperGUI.gold);
	files.get(i).add(null);
	pairedEndFiles.get(i).add(null);
	JLabel repBrowse = new JLabel(browseButton);
	repBrowse.setToolTipText("Search for file");
	repBrowse.addMouseListener((MouseListener)applet);
	replicateBrowseLabels.get(i).add(repBrowse);
	JLabel repRemove = new JLabel(replicateRemoveButton);
	repRemove.setToolTipText("Remove replicate");
	repRemove.addMouseListener((MouseListener)applet);
	replicateRemoveLabels.get(i).add(repRemove);
	JLabel pairedEndAdd = new JLabel(pairedEndNoButton);
	pairedEndAdd.setToolTipText("Using paired-end reads; add file of mate-pair reads");
	pairedEndAdd.addMouseListener((MouseListener)applet);
	pairedEndLabels.get(i).add(pairedEndAdd);
	experimentPanels.remove(i);
	experimentPanels.add(i, generateExperimentPanel(i));
	populateScrollPane();
    }

    public void removeReplicate(int i, int j) {
	scrollPane.remove(panel);
	fileNames.get(i).remove(j);
	files.get(i).remove(j);
	pairedEndFiles.get(i).remove(j);
	replicateBrowseLabels.get(i).remove(j);
	replicateRemoveLabels.get(i).remove(j);
	pairedEndLabels.get(i).remove(j);
	experimentPanels.remove(i);
	experimentPanels.add(i, generateExperimentPanel(i));
	populateScrollPane();
    }

    /**
     * Returns a String representation of the experiments' names and
     * all replicate files (including paths to files) for each experiment.
     * This String is output to a file, which can be read in to populate
     * the set of experiments and replicate files.
     *
     * There is one line per experiment. Each line is tab-delimited. 
     * A line contains the name of an experiment
     * followed by the path to each replicate file for the experiment.
     * If mate-pair files are included for paired-end reads, each
     * replicate file is followed by a comma and the path to the mate-pair file.
     */
    public String experimentsToString() {
	StringBuilder sb = new StringBuilder();
	for (int i=0; i<files.size(); i++) {

	    boolean validExperiment = false;
	    int count = 0;  // Count replicate file names for an experiment
	    for (int j=0; j<files.get(i).size(); j++) {
		if (files.get(i).get(j) != null) count++;
	    }
	    if (count  > 0) {
		// Output name of experiment
		String experimentName = experimentNames.get(i).getText();
		if (experimentName.length() == 0) experimentName = "N/A";
		if (files.get(i).size() > 0) {
		    sb.append(experimentName);
		    validExperiment = true;
		}
	    }
	    for (int j=0; j<files.get(i).size(); j++) {
		if (files.get(i).get(j) != null)
		    sb.append("\t" + files.get(i).get(j).getAbsolutePath());
		if (pairedEndFiles.get(i).get(j) != null)
		    sb.append("%" + pairedEndFiles.get(i).get(j).getAbsolutePath());
	    }
	    if (validExperiment) sb.append("\n");
	}
	return sb.toString();
    }

    public void experimentsFromString(String s) {
	String[] experiments = s.split("\n");
	for (int i=0; i<experiments.length; i++) {
	    String[] replicates = experiments[i].split("\t");
	    addExperiment();
	    experimentNames.get(i).setText(replicates[0]);
	    removeReplicate(i, 0);
	    for (int j=1; j<replicates.length; j++) {
		addReplicate(i);
		String[] parse_replicates = replicates[j].split("%");
		File f = new File(parse_replicates[0]);
		fileNames.get(i).get(j-1).setText(f.getName());
		fileNames.get(i).get(j-1).setBackground(Color.WHITE);
		files.get(i).set(j-1, f);
		if (parse_replicates.length > 1) {  // Using paired-end reads
		    pairedEndFiles.get(i).set(j-1, new File(parse_replicates[1]));
		    pairedEndLabels.get(i).get(j-1).setIcon(pairedEndYesButton);
		}
	    }
	}
    }

    /**
     * Returns a comma-separated list of experiment names.
     */
    public String getExperimentsList() {
	StringBuilder sb = new StringBuilder();
	if (experimentNames.size() > 0) sb.append(experimentNames.get(0).getText());
	for (int i=1; i<experimentNames.size(); i++)
	    sb.append("," + experimentNames.get(i).getText());
	return sb.toString();
    }

    /**
     * Returns an array where each entry corresponds to one experiment.
     * Each entry is a comma-separated list of replicate files for the experiment.
     * If paired-end reads are used, each comma-separated entry is a percent-separated
     * list of mate pair files.
     */
    public String[] getReplicatesList() {
	ArrayList<String> repsList = new ArrayList<String>();
	for (int i=0; i<files.size(); i++) {
	    ArrayList<String> reps = new ArrayList<String>();
	    for (int j=0; j<files.get(i).size(); j++) {
		if (files.get(i).get(j) != null) {
		    if (files.get(i).get(j).getAbsolutePath().length() > 0) {
			StringBuilder sb = new StringBuilder(files.get(i).get(j).getAbsolutePath());
			if (pairedEndFiles.get(i).get(j) != null) {
			    if (pairedEndFiles.get(i).get(j).getAbsolutePath().length() > 0)
				sb.append("%" + pairedEndFiles.get(i).get(j).getAbsolutePath());
			}
			reps.add(sb.toString());
		    }
		}
	    }
	    if (reps.size() > 0) {
		StringBuilder sb = new StringBuilder();
		sb.append(reps.get(0));
		for (int j=1; j<reps.size(); j++) sb.append("," + reps.get(j));
		repsList.add(sb.toString());
	    }
	}
	return repsList.toArray(new String[repsList.size()]);
    }

    private JPanel generateExperimentPanel(int i) {
	Color color = RockhopperGUI.lightPurple;
	int numReplicates = fileNames.get(i).size();
	JPanel newPanel = new JPanel(new GridLayout(numReplicates+2,1,0,0));
	newPanel.setPreferredSize(new Dimension(680, 30*numReplicates+70));
	newPanel.setBackground(color);
	JPanel p1 = new JPanel(new GridLayout(1,2,0,0));
	p1.setBackground(color);
	JPanel leftPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
	leftPanel.setBackground(color);
	JLabel expLabel = new JLabel("Experiment Name:");
	expLabel.setForeground(Color.red);
	leftPanel.add(expLabel);
	leftPanel.add(experimentNames.get(i));
	JPanel rightPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	rightPanel.setBackground(color);
	rightPanel.add(experimentRemoveLabels.get(i));
	p1.add(leftPanel);
	p1.add(rightPanel);
	newPanel.add(p1);
	for (int j=0; j<numReplicates; j++) {
	    JPanel p2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
	    p2.setBackground(color);
	    JLabel spacer2 = new JLabel("");
	    spacer2.setPreferredSize(new Dimension(70,0));
	    p2.add(spacer2);
	    JLabel repLabel = new JLabel("Replicate #" + (j+1) + ":");
	    repLabel.setForeground(Color.blue);
	    p2.add(repLabel);
	    p2.add(fileNames.get(i).get(j));
	    p2.add(replicateBrowseLabels.get(i).get(j));
	    p2.add(replicateRemoveLabels.get(i).get(j));
	    p2.add(pairedEndLabels.get(i).get(j));
	    newPanel.add(p2);
	}
	JPanel p3 = new JPanel(new FlowLayout(FlowLayout.LEFT));
	p3.setBackground(color);
	JLabel spacer3 = new JLabel("");
	spacer3.setPreferredSize(new Dimension(70, 0));
	p3.add(spacer3);
	p3.add(replicateAddLabels.get(i));
	newPanel.add(p3);
	return newPanel;
    }

    private void populateScrollPane() {
	panel = new JPanel(new FlowLayout(FlowLayout.CENTER,2,2));
	int height = 0;
	for (int i=0; i<experimentPanels.size(); i++) height += experimentPanels.get(i).getPreferredSize().height;
	panel.setPreferredSize(new Dimension(700, height + (experimentPanels.size()*14) + 40));
	panel.setBackground(RockhopperGUI.lightBlue);
	for (int i=0; i<experimentPanels.size(); i+=2) {
	    if (i > 0) panel.add(RockhopperGUI.getSpacerPanel(700, 10));
	    JPanel twoPanel = new JPanel(new FlowLayout(FlowLayout.CENTER,3,3));
	    int twoHeight = experimentPanels.get(i).getPreferredSize().height + 6; 
	    if (i < experimentPanels.size()-1) twoHeight += experimentPanels.get(i+1).getPreferredSize().height + 6;
	    else twoHeight += 3;
	    twoPanel.setPreferredSize(new Dimension(690, twoHeight));
	    twoPanel.setBackground(RockhopperGUI.lightBlue);
	    twoPanel.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.RAISED));
	    twoPanel.add(experimentPanels.get(i));
	    if (i < experimentPanels.size()-1) twoPanel.add(experimentPanels.get(i+1));
	    panel.add(twoPanel);
	}

	JPanel p4 = new JPanel(new FlowLayout(FlowLayout.LEFT,5,5));
	p4.setPreferredSize(new Dimension(700, 30));
	p4.setBackground(RockhopperGUI.lightBlue);
	p4.add(experimentAddLabel);
	panel.add(p4);
	scrollPane = new JScrollPane(panel);
	scrollPane.setPreferredSize(new Dimension(720, 500));
	scrollPane.setBorder(BorderFactory.createEmptyBorder());

	// Force scroll to bottom
	SwingUtilities.invokeLater( new Runnable() {
		public void run()
		{
		    JScrollBar v = scrollPane.getVerticalScrollBar();
		    v.setValue(v.getMaximum());
		}
	    } );
    }

    private void error(String title, String message) {
	JOptionPane.showMessageDialog(SwingUtilities.getAncestorOfClass(JFrame.class,applet), message, title, JOptionPane.ERROR_MESSAGE);
    }

}

