package edu.ucf.cs.hmm.squid;

import java.io.IOException;

abstract public class MSAReader {
	abstract public MSA Read(MSAFILE afp) throws IOException;
}
