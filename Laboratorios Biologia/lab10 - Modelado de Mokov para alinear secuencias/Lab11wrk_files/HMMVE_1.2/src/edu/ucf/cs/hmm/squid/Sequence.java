package edu.ucf.cs.hmm.squid;

public class Sequence {
	public SQINFO sqInfo;
	public String seqs;
	public double weight;
	public Sequence()
	{
		sqInfo = new SQINFO();
		seqs = "";
		weight = 0;
	}
}
