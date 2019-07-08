package edu.ucf.cs.hmm.squid;

import java.io.IOException;

/**
 * Adapted from Hammer squid library
 *
 * 
 */
public class IGReader extends SequenceReader {

	@Override
	boolean[] endTest(String s) {
		boolean[] result = new boolean[2];
		result[1] = true; /* 1 or 2 occur in line w/ bases */
		if (s.contains("1")||s.contains("2"))
			result[0] = true;
		else
			result[0] = false;
		return result;
	}

	@Override
	void read(SQFILE V) throws IOException, InvalidFlagException
	{
		String nm[];
			/* position past ';' comments */
		do {
		V.SeqfileGetLine();
		} while (! (V.buf==null) || ((V.buf != null) && (V.buf.trim().equals("")||(V.buf.startsWith(";")))));
		
		if ((V.buf != null))
		{
			if ((nm = V.buf.split("[\n|\t| ]+")) != null)
				V.sqinfo.SetSeqinfoString(nm[0], SQINFO.SQINFO_NAME);
			readLoop(false, V);
		}
		
		while (V.buf != null && V.buf.startsWith(";"))
			V.SeqfileGetLine();
	}
}
