package edu.ucf.cs.hmm.visualizer;

import java.io.IOException;
import java.io.StringWriter;

import javax.swing.JTextArea;

/**
 * Textarea of ascii version of hmm
 */
public class HmmText extends JTextArea {
	Hmm hmm;

	public HmmText() {
		super();
		setEditable(false);

	}

	public void setHmm(Hmm _hmm) {
		hmm = _hmm;
	}

	void updateHmmText() {
		if (hmm != null) {
			StringWriter w = new StringWriter();
			try {
				hmm.write(w);
				setText(w.toString());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
}
