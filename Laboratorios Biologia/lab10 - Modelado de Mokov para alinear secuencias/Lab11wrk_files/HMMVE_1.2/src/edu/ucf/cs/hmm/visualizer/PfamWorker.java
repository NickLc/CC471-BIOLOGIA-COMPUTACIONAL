package edu.ucf.cs.hmm.visualizer;

import java.util.List;

import javax.swing.JButton;
import javax.swing.JProgressBar;
import javax.swing.JTable;
import javax.swing.SwingWorker;

/**
 * 
 * Thread to read pfam index background
 */
public class PfamWorker extends SwingWorker<Pfam, Integer> implements IndexingListener {
	Pfam pfam;
	JTable table;
	JProgressBar progressBar;
	JButton enableButton;
	JButton disableButton;
	PfamTableModel model;
	public PfamWorker(Pfam p, JTable t, JProgressBar bar, JButton enable, JButton disable, PfamTableModel m)
	{
		pfam = p;
		table = t;
		progressBar = bar;
		enableButton = enable;
		disableButton = disable;
		progressBar.setMaximum(100);
		model = m;
	}

	@Override
	protected Pfam doInBackground() throws Exception {
		pfam.addIndexingListener(this);
		pfam.buildIndex();
		return pfam;
	}
	public void indexProgress(int i)
	{
		if (isCancelled())
			pfam.cancel();
		publish(i);
	}
	public void indexFinish()
	{
		publish(-1);
	}
	@Override
    protected void process(List<Integer> chunks) {
		int high = chunks.get(chunks.size()-1);
		if (high==-1)
		{
			enableButton.setEnabled(true);
			disableButton.setEnabled(false);
			
		}
		if (high==-1)
			high=pfam.indexEntries.size();
		progressBar.setValue((int)(((double)pfam.position)/pfam.fileSize*100));
		model.processNewEntries();
		table.revalidate();
    }
}
