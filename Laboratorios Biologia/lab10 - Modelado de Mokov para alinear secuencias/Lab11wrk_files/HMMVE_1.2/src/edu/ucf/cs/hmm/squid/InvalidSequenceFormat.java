package edu.ucf.cs.hmm.squid;

import java.io.IOException;

public class InvalidSequenceFormat extends IOException {
	String reason;
	public InvalidSequenceFormat()
	{
		reason = "";
	}
	public InvalidSequenceFormat(String format, Object... args)
	{
		reason = String.format(format, args);
	}
}
