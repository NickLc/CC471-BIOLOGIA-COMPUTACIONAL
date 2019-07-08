package edu.ucf.cs.hmm.squid;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
/**
 * Adapted from Hammer squid library
 *
 * 
 */
public class SQFILE {
	public enum SQFILE_FORMAT {
		SQFILE_UNKNOWN, SQFILE_IG, SQFILE_GENBANK, SQFILE_EMBL, SQFILE_GCG, SQFILE_STRIDER, SQFILE_FASTA, SQFILE_ZUKER, SQFILE_IDRAW, SQFILE_PIR, SQFILE_RAW, SQFILE_SQUID, SQFILE_GCGDATA, MSAFILE_UNKNOWN, MSAFILE_STOCKHOLM, MSAFILE_SELEX, MSAFILE_MSF, MSAFILE_CLUSTAL, MSAFILE_A2M, MSAFILE_PHYLIP, MSAFILE_EPS
	}

	static final int kStartLength = 500;

	BufferedReader reader; /* open file pointer                  */

	String fname; /* name of file; used for diagnostics */

	int linenumber; /* what line are we on in the file    */

	String buf; /* dynamically allocated buffer */

	int buflen; /* allocation length for buf                */

	int ssimode; /* SSI_OFFSET_I32 or SSI_OFFSET_I64        */

	long ssioffset; /* disk offset to last line read into buf  */

	long r_off; /* offset to start of record               */

	long d_off; /* offset to start of sequence data        */

	int rpl; /* residues per data line for this file; -1 if unset, 0 if invalid */

	int lastrpl; /* rpl on last line seen */

	int maxrpl; /* max rpl on any line of the file */

	int bpl; /* bytes per data line; -1 if unset, 0 if invalid */

	int lastbpl; /* bpl on last line seen */

	int maxbpl; /* max bpl on any line of the file */

	String seq; /* growing sequence during parse */

	SQINFO sqinfo; /* name, id, etc, gathered during parse */

	String sp;

	int seqlen; /* current sequence length */

	SQFILE_FORMAT format; /* format of seqfile we're reading. */

	MSA msa;

	MSAFILE afp;

	/** Read a line from a sequence file into V->buf
	 *           If the fgets() is NULL, sets V->buf[0] to '\0'.
	 *
	 */
	void SeqfileGetLine() throws IOException {
		buf = reader.readLine();
		linenumber++;
	}

	/**  Add a line of sequence to the growing string in V.
	 *
	 *   @param s sequence to add
	 */
	void addseq(String s) {
		int lineRpl = 0; /* valid residues per line */
		int lineBpl = 0; /* characters per line     */

		if (ssimode == -1) { /* Normal mode: keeping the seq */
			for (int i = 0; i < s.length(); i++) {
				if (!Util.isdigit(s.charAt(i)) && !Util.isspace(s.charAt(i))) {
					seq += s.charAt(i);
					seqlen += seq.length();
				}
			}
		} else /* else: indexing mode, discard the seq */
		{
			lineRpl = 0;
			for (int i = 0; i < s.length(); i++) {
				if (!Util.isdigit(s.charAt(i)) && !Util.isspace(s.charAt(i))) {
					lineRpl++;
				}
			}
			seqlen += lineRpl;
			bpl = s.length();

			if (rpl != 0) { /* 0 means we already know rpl is invalid       */
				if (lastrpl > 0) { /* we're on something that's not the first line */
					if (rpl > 0 && lastrpl != rpl)
						rpl = 0;
					else if (rpl == -1)
						rpl = lastrpl;
				}
				lastrpl = lineRpl;
				if (lineRpl > maxrpl)
					maxrpl = lineRpl; /* make sure we check max length of final lines */
			}
			if (bpl != 0) { /* 0 means we already know bpl is invalid       */
				if (lastbpl > 0) { /* we're on something that's not the first line */
					if (bpl > 0 && lastbpl != bpl)
						bpl = 0;
					else if (bpl == -1)
						bpl = lastbpl;
				}
				lastbpl = lineBpl;
				if (lineBpl > maxbpl)
					maxbpl = lineBpl; /* make sure we check max length of final lines */
			}
		} /* end of indexing mode of addseq(). */

	}

	/** Determine format of an open file.
	 *           Returns format code.
	 *           Rewinds the file.
	 *           
	 *           Autodetects the following unaligned formats:
	 *              SQFILE_FASTA
	 *              SQFILE_GENBANK
	 *              SQFILE_EMBL
	 *              SQFILE_GCG
	 *              SQFILE_GCGDATA
	 *              SQFILE_PIR
	 *           Also autodetects the following alignment formats:
	 *              MSAFILE_STOCKHOLM
	 *              MSAFILE_MSF
	 *              MSAFILE_CLUSTAL
	 *              MSAFILE_SELEX
	 *              MSAFILE_PHYLIP
	 *
	 *           Can't autodetect MSAFILE_A2M, calls it SQFILE_FASTA.
	 *           MSAFileFormat() does the opposite.
	 *
	 * @param     reader SQFILE
	 *           
	 * Return:   format code, or SQFILE_UNKNOWN if unrecognized
	 */
	SQFILE_FORMAT SeqfileFormat(BufferedReader reader) throws IOException {
		String buf;
		SQFILE_FORMAT fmt = SQFILE_FORMAT.SQFILE_UNKNOWN;
		int ndataline;
		boolean has_junk;

		buf = null;
		ndataline = 0;
		has_junk = false;
		while ((buf = reader.readLine()) != null) {
			if (buf.trim().equals(""))
				continue;

			/* Well-behaved formats identify themselves in first nonblank line.
			 */
			if (ndataline == 0) {
				if (buf.startsWith(">>>>") && buf.contains("Len: ")) {
					fmt = SQFILE_FORMAT.SQFILE_GCGDATA;
					return fmt;
				}

				if (buf.charAt(0) == '>') {
					fmt = SQFILE_FORMAT.SQFILE_FASTA;
					return fmt;
				}

				if (buf.startsWith("!!AA_SEQUENCE")
						|| buf.startsWith("!!NA_SEQUENCE")) {
					fmt = SQFILE_FORMAT.SQFILE_GCG;
					return fmt;
				}

				if (buf.startsWith("# STOCKHOLM 1.")) {
					fmt = SQFILE_FORMAT.MSAFILE_STOCKHOLM;
					return fmt;
				}

				if (buf.startsWith("CLUSTAL")
						&& buf.contains("multiple sequence alignment")) {
					fmt = SQFILE_FORMAT.MSAFILE_CLUSTAL;
					return fmt;
				}

				if (buf.startsWith("!!AA_MULTIPLE_ALIGNMENT")
						|| buf.startsWith("!!NA_MULTIPLE_ALIGNMENT")) {
					fmt = SQFILE_FORMAT.MSAFILE_MSF;
					return fmt;
				}

				/* PHYLIP id: also just a good bet */
				String[] tokens = buf.split("[ |\t|\n]+");
				if (tokens.length >= 2) {
					try {
						Integer.parseInt(tokens[0]);
						Integer.parseInt(tokens[1]);
						fmt = SQFILE_FORMAT.MSAFILE_PHYLIP;
						return fmt;
					} catch (NumberFormatException e) {
					}
				}
			}
			/* We trust that other formats identify themselves soon.
			 */
			/* dead giveaways for extended SELEX */
			if (buf.startsWith("#=AU") || buf.startsWith("#=ID")
					|| buf.startsWith("#=AC") || buf.startsWith("#=DE")
					|| buf.startsWith("#=GA") || buf.startsWith("#=TC")
					|| buf.startsWith("#=NC") || buf.startsWith("#=SQ")
					|| buf.startsWith("#=SS") || buf.startsWith("#=CS")
					|| buf.startsWith("#=RF")) {
				fmt = SQFILE_FORMAT.MSAFILE_SELEX;
				return fmt;
			}

			if (buf.startsWith("///") || buf.startsWith("ENTRY ")) {
				fmt = SQFILE_FORMAT.SQFILE_PIR;
				return fmt;
			}

			/* a ha, diagnostic of an (old) MSF file */
			if (buf.contains("..") && buf.contains("MSF:")
					&& buf.contains("Check:")) {
				fmt = SQFILE_FORMAT.MSAFILE_MSF;
				return fmt;
			}

			/* unaligned GCG (must follow MSF test!) */
			if (buf.contains(" Check: ") && buf.contains("..")) {
				fmt = SQFILE_FORMAT.SQFILE_GCG;
				return fmt;
			}

			if (buf.startsWith("LOCUS ") || buf.startsWith("ORIGIN ")) {
				fmt = SQFILE_FORMAT.SQFILE_GENBANK;
				return fmt;
			}

			if (buf.startsWith("ID   ") || buf.startsWith("SQ   ")) {
				fmt = SQFILE_FORMAT.SQFILE_EMBL;
				return fmt;
			}

			/* But past here, we're being desperate. A simple SELEX file is
			 * very difficult to detect; we can only try to disprove it.
			 */
			if (buf.contains("#") || buf.contains("%"))
				continue; /* skip comment lines */

			/* Disproof 1. Noncomment, nonblank lines in a SELEX file
			 * must have at least two space-delimited fields (name/seq)
			 */
			String tokens[] = buf.split("[ |\t|\n]+");
			if (tokens.length < 2)
				has_junk = true;

			/* Disproof 2. 
			 * The sequence field should look like a sequence.
			 */
			if (tokens.length > 1 && Seqtype(tokens[1]) == SQINFO.SEQTYPE.kOtherSeq)
				has_junk = true;

			ndataline++;
			if (ndataline == 300)
				break; /* only look at first 300 lines */
		}

		if (ndataline == 0)
			throw new IOException();

		/* If we've made it this far, we've run out of data, but there
		 * was at least one line of it; check if we've
		 * disproven SELEX. If not, cross our fingers, pray, and guess SELEX. 
		 */
		if (has_junk == true)
			fmt = SQFILE_FORMAT.SQFILE_UNKNOWN;
		else
			fmt = SQFILE_FORMAT.MSAFILE_SELEX;

		return fmt;
	}

	
	/** Open a sequence database file and prepare for reading
	 *           sequentially.
	 *           
	 * @param filename - name of file to open
	 * @param format   - format of file
	 */

	void SeqfileOpen(String filename, SQFILE_FORMAT _format) throws IOException
	{
		  rpl      = -1;		/* flag meaning "unset" */
		  lastrpl  = 0;
		  maxrpl   = 0;
		  bpl      = -1;		/* flag meaning "unset" */
		  lastbpl  = 0;
		  maxbpl   = 0;
		  ssimode = -1;

		  format = _format;
		  reader = new BufferedReader(new FileReader(filename));

	      fname    = filename;
		  

		  /* Invoke autodetection if we haven't already been told what
		   * to expect.
		   */
		  if (format == SQFILE_FORMAT.SQFILE_UNKNOWN)
		  {
		      format = SeqfileFormat(reader);
		      reader.close();
		      reader = new BufferedReader(new FileReader(filename));
		      if (format == SQFILE_FORMAT.SQFILE_UNKNOWN)
				throw new UnknownFormatException();
		  }

		  /* The hack for sequential access of an interleaved alignment file:
		   * read the alignment in, we'll copy sequences out one at a time.
		   */
		  msa        = null;
		  afp        = null;
		  linenumber = 0;
		  buf        = null;
		  if (IsAlignmentFormat(format))	
		  {
		      /* We'll be reading from the MSA interface. Copy our data
		       * to the MSA afp's structure.
		       */
		      afp           = new MSAFILE();
		      afp.reader        = reader;            /* just a ptr, don't close */
		      afp.fname    = fname;        /* just a ptr, don't free */
		      afp.format   = format;       /* e.g. format */
		      afp.linenumber = linenumber;	/* e.g. 0 */
		      afp.buf      = null;

		      if ((msa = afp.Read()) == null)
		    	  throw new IOException();
						/* hack: overload/reuse msa->lastidx; indicates
						   next seq to return upon a ReadSeq() call */
		      msa.lastidx = 0;

		      return;
		    }

		  /* Load the first line.
		   */
		  SeqfileGetLine(); 
		  return;
	}
	/** Returns a (very good) guess about type of sequence:
	 *           kDNA, kRNA, kAmino, or kOtherSeq.
	 *
	 */
	SQINFO.SEQTYPE Seqtype(String seq) {
		int saw; /* how many non-gap characters I saw */
		char c;
		String cs;
		int po = 0; /* count of protein-only */
		int nt = 0; /* count of t's */
		int nu = 0; /* count of u's */
		int na = 0; /* count of nucleotides */
		int aa = 0; /* count of amino acids */
		int no = 0; /* count of others */

		seq = seq.toUpperCase();

		/* Look at the first 300 non-gap characters
		 */
		for (saw = 0; saw < seq.length() && saw < 300; saw++) {
			c = seq.charAt(saw);
			cs = new String(new char[] { c });
			if (!Util.isgap(c)) {
				if ("EFIPQZ".contains(cs))
					po++;
				else if ("ACGTUN".contains(cs)) {
					na++;
					if (c == 'T')
						nt++;
					else if (c == 'U')
						nu++;
				} else if ("ABCDEFGHIKLMNPQRSTVWXYZ*".contains(cs))
					aa++;
				else if (Util.isalpha(c))
					no++;
				saw++;
			}
		}

		if (no > 0)
			return SQINFO.SEQTYPE.kOtherSeq;
		else if (po > 0)
			return SQINFO.SEQTYPE.kAmino;
		else if (na > aa) {
			if (nu > nt)
				return SQINFO.SEQTYPE.kRNA;
			else
				return SQINFO.SEQTYPE.kDNA;
		} else
			return SQINFO.SEQTYPE.kAmino; /* ooooh. risky. */
	}
	boolean IsAlignmentFormat(SQFILE_FORMAT format)
	{
		if (format == SQFILE_FORMAT.MSAFILE_STOCKHOLM ||
				format == SQFILE_FORMAT.MSAFILE_SELEX || format == SQFILE_FORMAT.MSAFILE_MSF ||
				format == SQFILE_FORMAT.MSAFILE_CLUSTAL || format == SQFILE_FORMAT.MSAFILE_A2M  ||
				format == SQFILE_FORMAT.MSAFILE_PHYLIP || format == SQFILE_FORMAT.MSAFILE_EPS)
			return true;
		return false;
	}
	
	/**  Given an aligned string of some type (either sequence or
	 *           secondary structure, for instance), dealign it relative
	 *           to a given aseq. Return the new string.
	 *
	 * @param aseq template alignment
	 * @param ss   string to make dealigned copy of; same length as aseq
	 * @return     dealigned copy of ss
	 */
	String MakeDealignedString(String aseq, String ss) throws SQParameterException
	{
	  String newString;
	  int   apos, rpos;

	  newString = new String();
	  for (apos = rpos = 0; apos < aseq.length(); apos++)
		  if (!Util.isgap(aseq.charAt(apos)))
	      {
	        newString += ss.charAt(apos);
	        rpos++;
	      }
	  if (aseq.length() != ss.length())
	  {
		  throw new SQParameterException();
	  }
	  return newString;
	}
	
	/** Read next sequence from an open database file.
	 *           Return the sequence
	 *           
	 * @return sequence
	 *                     
	 */           
	Sequence ReadSeq() throws InvalidSequenceFormat, SQParameterException, IOException, InvalidFlagException
	{
		Sequence result = new Sequence();
		boolean gotuw;

		/* Here's the hack for sequential access of sequences from
		 * the multiple sequence alignment formats
		 */
		if (IsAlignmentFormat(format))
		{
		    if (msa.lastidx >= msa.nseq) 
		    { /* out of data. try to read another alignment */				
		    	if ((msa = afp.Read()) == null)
		    		return null;
		    	msa.lastidx = 0;
		    }
					/* copy and dealign the appropriate aligned seq */
		    seq = MakeDealignedString(msa.aseq.get(msa.lastidx).toString(), msa.aseq.get(msa.lastidx).toString());
		    seqlen = seq.length();

		      /* Extract sqinfo stuff for this sequence from the msa.
		       * Tedious; code that should be cleaned.
		       */
		    result.sqInfo.flags = 0;
		    if (msa.sqname.get(msa.lastidx) != null) 
		    	result.sqInfo.SetSeqinfoString(msa.sqname.get(msa.lastidx), SQINFO.SQINFO_NAME);
		    if (msa.sqacc != null && msa.sqacc.get(msa.lastidx) != null) 
		    	result.sqInfo.SetSeqinfoString(msa.sqacc.get(msa.lastidx), SQINFO.SQINFO_ACC);
		    if (msa.sqdesc != null && msa.sqdesc.get(msa.lastidx) != null) 
		    	result.sqInfo.SetSeqinfoString(msa.sqdesc.get(msa.lastidx), SQINFO.SQINFO_DESC);
		    if (msa.ss != null && msa.ss.size()>msa.lastidx && msa.ss.get(msa.lastidx) != null && 
		    		!msa.ss.get(msa.lastidx).equals("")) {
		    	result.sqInfo.ss = MakeDealignedString(msa.aseq.get(msa.lastidx).toString(), 
		    			msa.ss.get(msa.lastidx).toString());
		    	result.sqInfo.flags |= SQINFO.SQINFO_SS;
		    }
		    if (msa.sa != null && msa.sa.size()>msa.lastidx && msa.sa.get(msa.lastidx) != null &&
		    		!msa.sa.get(msa.lastidx).equals("")) {
		    	result.sqInfo.sa = MakeDealignedString(msa.aseq.get(msa.lastidx).toString(), 
		    			msa.sa.get(msa.lastidx).toString());
		    	result.sqInfo.flags |= SQINFO.SQINFO_SA;
		    }
		    msa.lastidx++;
	    }
	  else {
	    if (buf==null) return null;

	    if (ssimode == -1) {	/* normal mode */
	      seq           = new String();
	    } else {			/* index mode: discarding seq */
	      seq           = null;
	    }
	    seqlen        = 0;
	    sqinfo        = result.sqInfo;
	    sqinfo.flags = 0;

	    SequenceReader r;
	    if (format==SQFILE_FORMAT.SQFILE_IG)
	    {
	    	r = new IGReader();
	    	r.read(this);
	    }
	    else if (format==SQFILE_FORMAT.SQFILE_STRIDER)
	    {
	    	//readStrider(V);
	    }
	    else if (format==SQFILE_FORMAT.SQFILE_GENBANK)
	    {
	    	//readGenBank(V);
	    }
	    else if (format==SQFILE_FORMAT.SQFILE_FASTA)
	    {
	    	r = new PearsonReader();
	    	r.read(this);
	    }
	    else if (format==SQFILE_FORMAT.SQFILE_EMBL)
	    {
	    	//readEMBL(V);
	    }
	    else if (format==SQFILE_FORMAT.SQFILE_ZUKER)
	    {
	    	//readZuker(V);
	    }
	    else if (format==SQFILE_FORMAT.SQFILE_PIR)
	    {
	    	//readPIR(V);
	    }
	    else if (format==SQFILE_FORMAT.SQFILE_GCGDATA)
	    {
	    	//readGCGdata(V);
	    } 
	    else if (format==SQFILE_FORMAT.SQFILE_GCG)
	    {
//	    	do {			/* skip leading comments on GCG file */
//	    		gotuw = (V.buf.contains("..")));
//	    		if (gotuw) readUWGCG(V);
//	    		V.SeqfileGetLine();
//	    	} while (V.buf != null);
	    }
	    else if (format==SQFILE_FORMAT.SQFILE_IDRAW)   /* SRE: no attempt to read idraw postscript */
	    {
	    }
	    else
	    {
	      throw new InvalidSequenceFormat();
	    }
	  }

	  /* Cleanup
	   */
	  result.sqInfo.len    = seqlen; 
	  result.sqInfo.flags |= SQINFO.SQINFO_LEN;
	  result.seqs = seq;
	  return result;
	}
	
	/**  Open a data file and
	 *           parse it into an array of rseqs (raw, unaligned
	 *           sequences).
	 *
	 *           Weights are currently only supported for MSF format.
	 *           Sequences read from all other formats will be assigned
	 *           weights of 1.0. If the caller isn't interested in
	 *           weights, it passes NULL as ret_weights.
	 *
	 */
	public ArrayList<Sequence> ReadMultipleRseqs(String seqfile, SQFILE_FORMAT fformat) throws InvalidSequenceFormat, SQParameterException, IOException, InvalidFlagException
	{
		ArrayList<Sequence> result = new ArrayList<Sequence>();
		SeqfileOpen(seqfile, fformat);

		Sequence s;
		while ((s=ReadSeq())!=null)
	    {
			result.add(s);
	    }
	  return result;
	}
	/** Try to guess whether an alignment is protein 
	 *           or nucleic acid; return a code for the
	 *           type (kRNA, kDNA, or kAmino).
	 *
	 * @param aseq array of aligned sequences. (Could also
	 *                   be an rseq unaligned sequence array)
	 * @param nseq  - number of aseqs
	 *
	 * Returns:  kRNA, kDNA, kAmino;
	 *           kOtherSeq if inconsistency is detected.
	 */
	static SQINFO.SEQTYPE
	GuessAlignmentSeqtype(ArrayList<StringBuffer> aseq, int nseq)
	{
	  int idx;
	  int nrna   = 0;
	  int ndna   = 0;
	  int namino = 0;
	  int nother = 0;

	  for (idx = 0; idx < nseq; idx++)
	  {
		  SQINFO.SEQTYPE atype = Seqtype(aseq.get(idx));
	    if (atype==SQINFO.SEQTYPE.kRNA)
	    	nrna++;
	    else if (atype==SQINFO.SEQTYPE.kDNA)
	    	ndna++;
	    else if (atype==SQINFO.SEQTYPE.kAmino)
	    	namino++;
	    else
	    	nother++;
	  }
	  /* Unambiguous decisions:
	   */
	  if (nother!=0)         return SQINFO.SEQTYPE.kOtherSeq;
	  if (namino == nseq) return SQINFO.SEQTYPE.kAmino;
	  if (ndna   == nseq) return SQINFO.SEQTYPE.kDNA;
	  if (nrna   == nseq) return SQINFO.SEQTYPE.kRNA;

	  /* Ambiguous decisions:
	   */
	  if (namino == 0)    return SQINFO.SEQTYPE.kRNA; /* it's nucleic acid, but seems mixed RNA/DNA */
	  return SQINFO.SEQTYPE.kAmino;		   /* some amino acid seen; others probably short seqs, some 
					      of which may be entirely ACGT (ala,cys,gly,thr). We 
					      could be a little more sophisticated: U would be a giveaway
					      that we're not in protein seqs */
	}
	/** Returns a (very good) guess about type of sequence:
	 *           kDNA, kRNA, kAmino, or kOtherSeq.
	 *           
	 */
	static SQINFO.SEQTYPE Seqtype(StringBuffer seq)
	{
	  int  saw=0;			/* how many non-gap characters I saw */
	  char c;
	  String cs;
	  int  po = 0;			/* count of protein-only */
	  int  nt = 0;			/* count of t's */
	  int  nu = 0;			/* count of u's */
	  int  na = 0;			/* count of nucleotides */
	  int  aa = 0;			/* count of amino acids */
	  int  no = 0;			/* count of others */
	  
	  /* Look at the first 300 non-gap characters
	   */
	  for (int i=0; i<seq.length() && saw < 300; i++)
	    {
	      c = seq.substring(i, i+1).toUpperCase().charAt(0);
	      if (! Util.isgap(c)) 
		{
	    	  cs = new String(new char[] { c });
		  if ("EFIPQZ".contains(cs)) po++;
		  else if ("ACGTUN".contains(cs)) {
		    na++;
		    if (c == 'T') nt++;
		    else if (c == 'U') nu++;
		  }
		  else if ("ABCDEFGHIKLMNPQRSTVWXYZ*".contains(cs)) aa++;
		  else if (Util.isalpha(c)) no++;
		  saw++;
		}
	    }

	  if (no > 0) return SQINFO.SEQTYPE.kOtherSeq;
	  else if (po > 0) return SQINFO.SEQTYPE.kAmino;
	  else if (na > aa) {
	    if (nu > nt) return SQINFO.SEQTYPE.kRNA;
	    else return SQINFO.SEQTYPE.kDNA;
	    }
	  else return SQINFO.SEQTYPE.kAmino;		/* ooooh. risky. */
	}
	
	/** GCG checksum for a multiple alignment: sum of
	 *           individual sequence checksums (including their
	 *           gap characters) modulo 10000.
	 *
	 *           Implemented using spec provided by Steve Smith of
	 *           Genetics Computer Group.
	 *           
	 * @param seqs sequences to be checksummed; aligned or not
	 * @param nseq number of sequences
	 *           
	 * @return     the checksum, a number between 0 and 9999
	 */                      
	static int GCGMultchecksum(StringBuffer[] seqs, int nseq)
	{
	  int chk = 0;
	  int idx;

	  for (idx = 0; idx < nseq; idx++)
	    chk = (chk + GCGchecksum(seqs[idx], seqs[idx].length())) % 10000;
	  return chk;
	}
	
	/** Calculate a GCG checksum for a sequence.
	 *           Code provided by Steve Smith of Genetics
	 *           Computer Group.
	 *
	 * @param seq sequence to calculate checksum for.
	 *                  may contain gap symbols.
	 * @param len length of sequence (usually known,
	 *                  so save a strlen() call)       
	 *
	 * @return  GCG checksum.
	 */
	static int GCGchecksum(StringBuffer seq, int len)
	{
	  int i;			/* position in sequence */
	  int chk = 0;			/* calculated checksum  */

	  for (i = 0; i < len; i++) 
	    chk = (chk + (i % 57 + 1) * ((int)seq.charAt(i))) % 10000;
	  return chk;
	}
}
