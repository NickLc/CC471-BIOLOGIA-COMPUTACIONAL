package edu.ucf.cs.hmm.squid;

import java.io.IOException;
/**
 * Adapted from Hammer squid library
 *
 * 
 */
public class PearsonReader extends SequenceReader {
	@Override
	boolean[] endTest(String s) {
		boolean[] result = new boolean[2];
		result[0] = (s==null || s.startsWith(">"));
		result[1] = false;
		return result;
	}
	@Override
	void read(SQFILE V) throws IOException, InvalidFlagException
	{
		  if (V.ssimode >= 0) V.r_off = V.ssioffset;

		  if (!V.buf.startsWith(">")) 
		    throw new InvalidSequenceFormat();

		  String[] tokens = V.buf.substring(1).trim().split("[\n|\t| ]+");
		  if (tokens.length>0)
			  V.sqinfo.SetSeqinfoString(tokens[0], SQINFO.SQINFO_NAME);
		  if (tokens.length>1)
			  V.sqinfo.SetSeqinfoString(tokens[1], SQINFO.SQINFO_DESC);

		  readLoop(false, V);

		  while (V.buf!=null && (!V.buf.startsWith(">"))) 
		  {
			  V.SeqfileGetLine();
		  }
	}
}
