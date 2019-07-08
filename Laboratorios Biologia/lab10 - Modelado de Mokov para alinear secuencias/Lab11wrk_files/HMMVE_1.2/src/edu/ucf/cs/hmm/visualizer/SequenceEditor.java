package edu.ucf.cs.hmm.visualizer;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
/**
 * A simple sequence editor with ruler
 */

class JRulerPanel extends JPanel
{
	JTextArea editor;
	String seq;
	int txtWidth;
	JRulerPanel(String _seq, boolean editable)
	{
		seq = _seq;
		editor = new JTextArea();
		Font font = new Font("Courier", Font.PLAIN, 14);
		editor.setFont(font);
		editor.setText(seq);
		editor.setEditable(editable);
		FontMetrics metrics = editor.getFontMetrics(font);
		txtWidth = metrics.stringWidth(" ");

		editor.setColumns(80);
		editor.setLineWrap(true);

		JScrollPane scrollPane = new JScrollPane(editor,
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

		setLayout(null);
		scrollPane.setBounds(10, 30, txtWidth * 80+10, 420);
		add(scrollPane);
		setBounds(10, 0, txtWidth * 80+30, 460);
	}
	@Override
	public void paintComponent(Graphics g)
	{
		Graphics2D g2d = (Graphics2D) g;
		g2d.setColor(new Color(0x00, 0x00, 0x00));
		g2d.drawLine(10, 2, txtWidth * 80 + 10, 2);
		for (int i = 0; i <= 80; i++) {
			if (i % 10 == 0) {
				g2d.drawLine(txtWidth * i + 10, 2, txtWidth * i + 10, 8);
				String s = "" + i;
				g2d.drawString(s,
						txtWidth * i - s.length() * txtWidth / 2 + 10, 18);
			} else
				g2d.drawLine(txtWidth * i + 10, 2, txtWidth * i + 10, 4);
		}
	}
}

public class SequenceEditor extends JDialog implements ActionListener {

	JRulerPanel rulerPanel;
	
	String text;

	String alphabet;

	public SequenceEditor(String seq, String _alphabet, boolean editable) {
		super((Frame) null, "Sequence Editor", true);
		alphabet = _alphabet;
		rulerPanel = new JRulerPanel(seq, editable);
		getContentPane().setLayout(null);
		setSize(rulerPanel.txtWidth * 80 + 60, 520);
		getContentPane().add(rulerPanel);

		JButton okButton = new JButton("Ok");
		okButton.setMnemonic(KeyEvent.VK_O);
		okButton.setBounds(60, 460, 80, 20);
		okButton.setActionCommand("Ok");
		okButton.addActionListener(this);
		getContentPane().add(okButton);

		JButton cancelButton = new JButton("Cancel");
		cancelButton.setMnemonic(KeyEvent.VK_C);
		cancelButton.setBounds(rulerPanel.txtWidth * 80 - 80 - 60, 460, 80, 20);
		cancelButton.setActionCommand("Cancel");
		cancelButton.addActionListener(this);
		getContentPane().add(cancelButton);
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("Ok")) {
			if (alphabet != null) {
				text = rulerPanel.editor.getText();
				for (int i = 0; i < text.length(); i++) {
					if (!alphabet.contains("" + text.charAt(i))) {
						JOptionPane.showMessageDialog(this, "Invalid alphabet");
						rulerPanel.editor.setCaretPosition(i);
						rulerPanel.editor.requestFocusInWindow();
						return;
					}
				}
			}
		} else
			text = null;
		setVisible(false);
	}
}
