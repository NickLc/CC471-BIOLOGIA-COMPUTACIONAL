package edu.ucf.cs.hmm.visualizer;

import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * Panel inside JModifyNodedialog for probability editing
 */
public class JProbEditorPanel extends JPanel implements FocusListener,
		ActionListener {
	Probabilities probs;

	JLabel sumText;

	JTextField textFields[];

	JLabel modifiedFields[];

	StateChartLabel chartLabel;

	public JProbEditorPanel(Probabilities p) {
		if (p == null)
			return;
		probs = p;
		setLayout(null);
		setSize(new Dimension(350, Math.min(p.size, 20) * 22 + 50));
		textFields = new JTextField[p.size];
		modifiedFields = new JLabel[p.size];
		for (int i = 0; i < p.size; i++) {
			JLabel l = new JLabel(p.names[i]);
			l.setBounds(new Rectangle(10, 22 * i + 1, 30, 20));
			add(l);
			JTextField t = new JTextField("" + p.probs[i]);
			t.setBounds(new Rectangle(50, 22 * i + 1, 150, 20));
			add(t);
			textFields[i] = t;
			t.addFocusListener(this);
			JLabel modified = new JLabel("");
			modified.setBounds(new Rectangle(210, 22 * i + 1, 10, 20));
			add(modified);
			modifiedFields[i] = modified;
		}
		JLabel sumLabel = new JLabel("Sum");
		sumLabel.setBounds(new Rectangle(10, 22 * p.size + 1, 30, 20));
		add(sumLabel);
		sumText = new JLabel("" + probs.sum());
		sumText.setBounds(new Rectangle(50, 22 * p.size + 1, 150, 20));
		add(sumText);
		chartLabel = new StateChartLabel();
		chartLabel.setAlphabet(probs.names);
		chartLabel.setState("", p.probs);
		chartLabel.setBounds(230, 0, 100, 100);
		add(chartLabel);
		JButton adjButton = new JButton("Normalize");
		adjButton.setBounds(230, 100, 100, 30);
		add(adjButton);
		adjButton.setActionCommand("Normalize");
		adjButton.addActionListener(this);
		setPreferredSize(new Dimension(350, p.size * 22 + 50));
	}

	public void focusGained(FocusEvent e) {
	}

	public void focusLost(FocusEvent e) {
		for (int i = 0; i < probs.size; i++) {
			if (e.getComponent() == textFields[i]) {

				try {
					if (probs.probs[i] != Double.parseDouble(textFields[i]
							.getText())) {
						probs.probs[i] = Double.parseDouble(textFields[i]
								.getText());
						probs.modified[i] = true;
						modifiedFields[i].setText("*");

					} else
						modifiedFields[i].setText(" ");
					chartLabel.repaint();
				} catch (NumberFormatException e1) {
					modifiedFields[i].setText("?");
					textFields[i].requestFocus();
					return;
				}
			}
		}
		sumText.setText("" + probs.sum());
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("Normalize")) {
			probs.normalize();
			for (int i = 0; i < probs.size; i++) {
				if (probs.probs[i] != Double.parseDouble(textFields[i]
						.getText())) {
					textFields[i].setText("" + probs.probs[i]);
					probs.modified[i] = true;
					modifiedFields[i].setText("*");
				}
			}
			sumText.setText("" + probs.sum());
			chartLabel.repaint();
		}
	}
}
