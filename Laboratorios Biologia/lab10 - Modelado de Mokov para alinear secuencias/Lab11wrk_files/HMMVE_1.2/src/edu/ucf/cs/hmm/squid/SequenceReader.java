package edu.ucf.cs.hmm.squid;

import java.io.IOException;

/**
 * Adapted from Hammer squid library
 *
 * 
 */
abstract public class SequenceReader {
	void readLoop(boolean addfirst, SQFILE V) throws IOException
	{
	  boolean addend = false;
	  boolean done   = false;

	  V.seqlen = 0;
	  V.lastrpl = V.lastbpl = 0;
	  if (addfirst) {
	    if (V.ssimode >= 0) V.d_off = V.ssioffset;
	    V.addseq(V.buf);
	  } else if (V.ssimode >= 0)
		  throw new NotImplementedException();
	  do {
	    V.SeqfileGetLine();
	    if (V.buf == null)
	      done = true;
	    boolean[] result = endTest(V.buf);
	    done |= result[0];
	    addend = result[1];
	    if (addend || !done)
	      V.addseq(V.buf);
	  } while (!done);
	}
	abstract void read(SQFILE V) throws IOException, InvalidFlagException;
	abstract boolean[] endTest(String s);
}
