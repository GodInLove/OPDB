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

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.MouseListener;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeEvent;
import javax.swing.*;

public class ProgressGUI implements PropertyChangeListener {

    private JFrame progressFrame;
    private JApplet applet;
    public JTextArea progressArea;
    public JProgressBar progressBar;
    private JLabel IGV_Label;
    private ImageIcon IGV_Button, IGV_ButtonPressed;
    private int width = 620;
    private int height = 650;



    /**
     * Constructor.
     */
    public ProgressGUI(JApplet applet) {
	this.applet = applet;
    }



    public void createProgressWindow() {

	progressFrame = new JFrame("Progress");
	progressFrame.setSize(width,height);
	progressFrame.setResizable(false);
	progressFrame.setLocationRelativeTo(null);
	progressFrame.setIconImage(new ImageIcon(RockhopperGUI.penguin).getImage());
	progressFrame.setVisible(true);
	JPanel progressPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
	progressPanel.setPreferredSize(new Dimension(width,height));
	progressPanel.setBackground(RockhopperGUI.lightBlue);
	
	progressBar = new JProgressBar(0, 100);
	progressBar.setValue(0);
	progressBar.setStringPainted(true);

	progressArea = new JTextArea();
	progressArea.setEditable(false);
	progressArea.setText("Initializing RNAseq analysis...\n");

	IGV_Button = RockhopperGUI.getScaledImage("images/IGV_button.gif", 130, 57);
	IGV_ButtonPressed = RockhopperGUI.getScaledImage("images/IGV_buttonPressed.gif", 130, 57);
	IGV_Label = new JLabel(IGV_Button);
	IGV_Label.setToolTipText("Launch the Integrative Genomics Viewer from the Broad Institute");
	IGV_Label.addMouseListener((MouseListener)applet);
	IGV_Label.setVisible(false);

	JScrollPane scrollPane = new JScrollPane(progressArea);
	scrollPane.setPreferredSize(new Dimension(width-20,height-130));
	scrollPane.setBorder(BorderFactory.createEmptyBorder());

	progressPanel.add(progressBar);
	progressPanel.add(scrollPane);
	progressPanel.add(IGV_Label);
	progressFrame.add(progressPanel);

	progressFrame.validate();
	progressFrame.repaint();
    }

    /**
     * Add button to launch IGV.
     */
    public void addIGV() {
	IGV_Label.setVisible(true);
    }

    public void propertyChange(PropertyChangeEvent event) {
	if ("progress" == event.getPropertyName()) {
	    int progress = (Integer) event.getNewValue();
	    progressBar.setValue(progress);
	}
    }

    public boolean isSource(Object source) {
	if (source.equals(IGV_Label)) return true;
	return false;
    }

    public void mousePressed(Object source) {
	if (source.equals(IGV_Label)) {
	    IGV_Label.setIcon(IGV_ButtonPressed);
	}
    }

    public void mouseReleased(Object source) {
	if (source.equals(IGV_Label)) {
	    IGV_Label.setIcon(IGV_Button);
	}
    }

}

