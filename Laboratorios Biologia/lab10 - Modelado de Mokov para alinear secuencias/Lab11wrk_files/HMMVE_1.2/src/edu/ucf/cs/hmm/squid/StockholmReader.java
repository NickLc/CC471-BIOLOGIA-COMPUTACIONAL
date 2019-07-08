package edu.ucf.cs.hmm.squid;

import java.io.IOException;
/**
 * Adapted from Hammer squid library
 *
 * 
 */
public class StockholmReader extends MSAReader {
	/* Format of a GF line:
	 *    #=GF <featurename> <text>
	 */
	int	parse_gf(MSA msa, String buf)
	{
	  String gf;
	  String featurename;
	  String text;

	  String[] tokens = buf.split("[ |\t|\n]+");
	  if (tokens.length<3)
	  {
		  return 0;
	  }
	  gf = tokens[0];
	  featurename = tokens[1].trim();
	  text = tokens[1].trim();
	  
	  if (featurename.equals("ID")) 
		  msa.name = text;
	  else if (featurename.equals("AC")) 
		  msa.acc = text;
	  else if (featurename.equals("DE")) 
		  msa.desc = text;
	  else if (featurename.equals("AU")) 
		  msa.au = text;
	  else if (featurename.equals("GA")) 
	  {				/* Pfam has GA1, GA2. Rfam just has GA1. */
	      String[] tokens1 = text.split(" |\t|\n");
		  if (tokens1.length < 1)
			  return 0;
		  msa.cutoff[MSA.MSA_CUTOFF_GA1]        = Double.parseDouble(tokens1[0]);
		  msa.cutoff_is_set[MSA.MSA_CUTOFF_GA1] = true;
	      if (tokens1.length > 1) {
	    	  msa.cutoff[MSA.MSA_CUTOFF_GA2]        = Double.parseDouble(tokens1[1]);
	    	  msa.cutoff_is_set[MSA.MSA_CUTOFF_GA2] = true;
	      }
	  }
	  else if (featurename.equals("NC")) 
	  {
		  String[] tokens1 = text.split("[ |\t|\n]+");
	      if (tokens1.length < 1)
	    	  return 0;
	      msa.cutoff[MSA.MSA_CUTOFF_NC1]        = Double.parseDouble(tokens1[0]);
	      msa.cutoff_is_set[MSA.MSA_CUTOFF_NC1] = true;
	      if (tokens1.length > 1) {
	    	  msa.cutoff[MSA.MSA_CUTOFF_NC2]        = Double.parseDouble(tokens1[1]);
	    	  msa.cutoff_is_set[MSA.MSA_CUTOFF_NC2] = true;
	      }
	  }
	  else if (featurename.equals("TC")) 
	  {
		  String[] tokens1 = text.split("[ |\t|\n]+");
	      if (tokens1.length < 1)
	    	  return 0;
	      msa.cutoff[MSA.MSA_CUTOFF_TC1]        = Double.parseDouble(tokens1[0]);
	      msa.cutoff_is_set[MSA.MSA_CUTOFF_TC1] = true;
	      if (tokens1.length > 1) {
	    	  msa.cutoff[MSA.MSA_CUTOFF_TC2]        = Double.parseDouble(tokens1[1]);
	    	  msa.cutoff_is_set[MSA.MSA_CUTOFF_TC2] = true;
	      }
	  }
	  else 
	    msa.MSAAddGF(featurename, text);

	  return 1;
	}
	
	int parse_sequence(MSA msa, String buf)
	{
	  String seqname;
	  String text;
	  int   seqidx;
	  int   len;

	  String[] tokens = buf.split("[ |\t|\n]+");
	  if (tokens.length<2)
		  return 0;
	  seqname = tokens[0];
	  text = tokens[1];
	  
	  /* seq usually follows another seq; guess msa->lastidx +1 */
	  seqidx = msa.MSAGetSeqidx(seqname, msa.lastidx+1);
	  msa.lastidx = seqidx;

	  msa.aseq.set(seqidx, new StringBuffer(msa.aseq.get(seqidx)+text));
	  return 1;
	}

	
	/** Parse the next alignment from an open Stockholm
	 *           format alignment file. Return the alignment, or
	 *           NULL if there are no more alignments in the file.
	 *
	 * @param afp open alignment file
	 *
	 */
	@Override
	public MSA Read(MSAFILE afp) throws IOException
	{
		  MSA   msa = new MSA();
		  String  s;
		  int    status = 1;

		  /* Initialize allocation of the MSA.
		   */

		  /* Check the magic Stockholm header line.
		   * We have to skip blank lines here, else we perceive
		   * trailing blank lines in a file as a format error when
		   * reading in multi-record mode.
		   */
		  do {
		    if ((s = afp.MSAFileGetLine()) == null) {
		      return null;
		    }
		  } while (!s.trim().equals(""));

		  if (s.startsWith("# STOCKHOLM 1."))
		    throw new InvalidSequenceFormat();

		  /* Read the alignment file one line at a time.
		   */
		  while ((s = afp.MSAFileGetLine()) != null) 
		    {
		      s = s.trim();

		      if (s.startsWith("#")) {
		    	  if      (s.startsWith("#=GF"))   status = parse_gf(msa, s);
		    	  else if (s.startsWith("#=GS"))   ; // status = parse_gs(msa, s);
		    	  else if (s.startsWith("#=GC"))   ; // status = parse_gc(msa, s);
		    	  else if (s.startsWith("#=GR"))   ; // status = parse_gr(msa, s);
		    	  else							   ; // status = parse_comment(msa, s);
		      }
		      else if (s.startsWith("//"))   break;
		      else if (s.equals(""))                   continue;
		      else          status = parse_sequence(msa, s);

		      if (status == 0)  
		    	  throw new InvalidSequenceFormat();
		    }

		  if (s == null && msa.nseq != 0)
			  throw new InvalidSequenceFormat();

		  if (s == null && msa.nseq == 0) {
		    				/* probably just some junk at end of file */
		      return null; 
		    }
		  
		  // MSAVerifyParse(msa);
		  return msa;
	}
}
