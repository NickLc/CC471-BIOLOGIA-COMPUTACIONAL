package edu.ucf.cs.hmm.visualizer;

import java.awt.Font;
import java.awt.Frame;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.border.Border;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

/**
 * 
 * Dialog for Pfam indexer
 *
 */

public class PfamPickDialog extends JDialog implements ActionListener, DocumentListener, MouseListener {
	Pfam pfam;
	PfamWorker worker=null;
	JProgressBar progressBar=null;
	JTextField filterField=null;
	JButton reindexButton=null;
	JButton stopButton=null;
	JButton selectButton=null;
	JButton cancelButton=null;
	JTable pfamIndexTable=null;
	PfamTableModel pfamTableModel=null;
	int selectedHmmIndex=-1;
	public PfamPickDialog(Pfam p)
	{
		super((Frame) null, "Pfam Indexer", true);
		pfam = p;
		pfamIndexTable = new JTable();
		pfamIndexTable.setFont(new Font("Courier New", Font.BOLD, 12));
		setLayout(null);
		setSize(600, 800);
		addWindowListener(new WindowAdapter()
		{
			public void windowClosed(WindowEvent e)
			{
				if (stopButton!=null && stopButton.isEnabled())
					actionPerformed(new ActionEvent(this, 0, "Stop"));
				dispose();
			}
		});
		
		JScrollPane pfamIndexScrollPane = new JScrollPane(pfamIndexTable,
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		pfamTableModel = new PfamTableModel(pfam);
		pfamIndexTable.setModel(pfamTableModel);
		pfamIndexTable.getColumnModel().getColumn(0).setPreferredWidth(
				70);
		pfamIndexTable.getColumnModel().getColumn(1).setPreferredWidth(
				100);
		pfamIndexTable.getColumnModel().getColumn(2).setPreferredWidth(
				430);
		pfamIndexTable.addMouseListener(this);
		pfamIndexScrollPane.setBounds(new Rectangle(10, 30, getWidth()-20,
				getHeight() - 150));
		add(pfamIndexScrollPane);
		
		progressBar = new JProgressBar();
		progressBar.setBounds(new Rectangle(10, 5, getWidth()-20, 20));
		Border pbBorder = progressBar.getBorder();
		Border spaceBelow = BorderFactory.createEmptyBorder(0, 0, 5, 0);
        progressBar.setBorder(BorderFactory.createCompoundBorder(
                                        spaceBelow,
                                        pbBorder));
		add(progressBar);
		
		JLabel filterLabel = new JLabel("Filter");
		filterLabel.setBounds(10, getHeight()-110, 60, 20);
		add(filterLabel);
		
		filterField = new JTextField();
		filterField.setBounds(100, getHeight()-110, 400, 20);
		filterField.setActionCommand("Filter");
		filterField.getDocument().addDocumentListener(this);
		add(filterField);
		
		reindexButton = new JButton("ReIndex");
		reindexButton.setBounds(new Rectangle(35, getHeight()-70, 80, 30));
		reindexButton.setActionCommand("ReIndex");
		reindexButton.addActionListener(this);
		add(reindexButton);
		
		stopButton = new JButton("Stop");
		stopButton.setBounds(new Rectangle(185, getHeight()-70, 80, 30));
		stopButton.setActionCommand("Stop");
		stopButton.addActionListener(this);
		add(stopButton);
		
		selectButton = new JButton("Select");
		selectButton.setBounds(new Rectangle(335, getHeight()-70, 80, 30));
		selectButton.setActionCommand("Select");
		selectButton.addActionListener(this);
		selectButton.setEnabled(false);
		add(selectButton);
		
		cancelButton = new JButton("Cancel");
		cancelButton.setBounds(new Rectangle(485, getHeight()-70, 80, 30));
		cancelButton.setActionCommand("Cancel");
		cancelButton.addActionListener(this);
		add(cancelButton);
		
		if (!pfam.complete)
		{
			reindexButton.setEnabled(false);
			stopButton.setEnabled(true);
			worker = new PfamWorker(pfam, pfamIndexTable, progressBar, reindexButton, stopButton, pfamTableModel);
			worker.execute();
		}
	}
	public void actionPerformed(ActionEvent e)
	{
		if (e.getActionCommand().equals("ReIndex"))
		{
			reindexButton.setEnabled(false);
			stopButton.setEnabled(true);
			selectButton.setEnabled(false);
			pfam.reset();
			pfamTableModel.reset();
			worker = new PfamWorker(pfam, pfamIndexTable, progressBar, reindexButton, stopButton, pfamTableModel);
			worker.execute();
			filter();
		}
		else if (e.getActionCommand().equals("Stop"))
		{
			if (worker!=null && !worker.isDone())
			{
				worker.cancel(true);
				reindexButton.setEnabled(true);
				stopButton.setEnabled(false);
			}
		}
		else if (e.getActionCommand().equals("Select"))
		{
			if (pfamIndexTable.getSelectedRow()!=-1)
			{
				if (pfamTableModel.filter==null)
					selectedHmmIndex=pfamIndexTable.getSelectedRow();
				else
					selectedHmmIndex=pfamTableModel.filtered.get(pfamIndexTable.getSelectedRow());
				if (stopButton!=null && stopButton.isEnabled())
					actionPerformed(new ActionEvent(this, 0, "Stop"));
				dispose();
			}
		}
		else if (e.getActionCommand().equals("Cancel"))
		{
			if (stopButton!=null && stopButton.isEnabled())
				actionPerformed(new ActionEvent(this, 0, "Stop"));
			dispose();
		}
	}
	public void changedUpdate(DocumentEvent e)
	{
		filter();
	}
	public void insertUpdate(DocumentEvent e)
	{
		filter();
	}
	public void removeUpdate(DocumentEvent e)
	{
		filter();
	}
	public void filter()
	{
		pfamTableModel.setFilter(filterField.getText());
		pfamIndexTable.revalidate();
		pfamIndexTable.repaint();
	}
	public void mouseClicked(MouseEvent e) {
		if (e.getSource() == pfamIndexTable) {
			if (e.getClickCount() == 2 && pfamIndexTable.getSelectedRow() != -1)
			{
				actionPerformed(new ActionEvent(this, 0, "Select"));
			}
			if (pfamIndexTable.getSelectedRow()!=-1)
				selectButton.setEnabled(true);
			else
				selectButton.setEnabled(false);
		}
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseExited(MouseEvent e) {
	}

	public void mousePressed(MouseEvent e) {
	}

	public void mouseReleased(MouseEvent e) {
	}
}
