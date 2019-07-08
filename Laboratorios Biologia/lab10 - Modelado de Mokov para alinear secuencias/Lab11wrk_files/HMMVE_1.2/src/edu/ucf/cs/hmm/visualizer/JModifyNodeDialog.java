package edu.ucf.cs.hmm.visualizer;

import java.awt.Dimension;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;

/**
 * Dialog to modify emission and transition probability
 */
public class JModifyNodeDialog extends JDialog implements ActionListener {
	int action = JOptionPane.CANCEL_OPTION;

	Probabilities tp;

	Probabilities ep;

	public JModifyNodeDialog(String node, Probabilities _tp, Probabilities _ep) {
		super((Frame) null, node, true);
		tp = _tp;
		ep = _ep;
		setResizable(false);
		JProbEditorPanel transitionPanel = new JProbEditorPanel(tp);
		JScrollPane transitionScrollPane = new JScrollPane(transitionPanel,
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		JProbEditorPanel emissionPanel = new JProbEditorPanel(ep);
		JScrollPane emissionScrollPane = new JScrollPane(emissionPanel,
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		JTabbedPane tabbedPane = new JTabbedPane();
		tabbedPane.addTab("Transition", transitionScrollPane);
		tabbedPane.addTab("Emission", emissionScrollPane);

		int width = Math.max(transitionPanel.getWidth(), emissionPanel
				.getWidth());
		int height = Math.max(transitionPanel.getHeight(), emissionPanel
				.getHeight());
		if (height < 120)
			height = 160;

		tabbedPane.setSize(new Dimension(width, height));
		getContentPane().setLayout(null);
		getContentPane().add(tabbedPane);

		setSize(new Dimension(width, height + 80));
		JButton saveButton = new JButton("Save");
		saveButton.setMnemonic(KeyEvent.VK_S);
		saveButton.setBounds(60, height + 10, 80, 20);
		saveButton.setActionCommand("Save");
		getContentPane().add(saveButton);
		saveButton.addActionListener(this);
		JButton cancelButton = new JButton("Cancel");
		cancelButton.setMnemonic(KeyEvent.VK_C);
		cancelButton.setBounds(200, height + 10, 80, 20);
		getContentPane().add(cancelButton);
		cancelButton.setActionCommand("Cancel");
		cancelButton.addActionListener(this);
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("Save")) {
			if (tp != null && Math.abs(tp.sum() - 1) > 0.0001 || ep != null
					&& Math.abs(ep.sum() - 1) > 0.0001) {
				JOptionPane.showMessageDialog(this,
						"Sum is not 1. Please normalize it.");
				return;
			}
			action = JOptionPane.OK_OPTION;
			setVisible(false);
		} else if (e.getActionCommand().equals("Cancel")) {
			action = JOptionPane.CANCEL_OPTION;
			setVisible(false);
		}
	}
}
