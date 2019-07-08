package edu.ucf.cs.hmm.squid;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Adapted from Hammer squid library
 *
 * 
 */

public class MSA {
	  
	  static final int MSA_SET_WGT=(1 << 0);  /* track whether wgts were set, or left at default 1.0 */
	  /* Mandatory information associated with the alignment.
	   */
	  public ArrayList<StringBuffer> aseq;                  /* the alignment itself, [0..nseq-1][0..alen-1] */
	  public ArrayList<String> sqname;                /* names of sequences, [0..nseq-1][0..alen-1]   */
	  public ArrayList<Double> wgt;                   /* sequence weights [0..nseq-1]                 */
	  public int    alen;                  /* length of alignment (columns)                */
	  public int    nseq;                  /* number of seqs in alignment                  */

	  /* Optional information that we understand, and might have.
	   */
	  int    flags;                 /* flags for what optional info is valid    */
	  SQINFO.SEQTYPE    type;                  /* kOtherSeq, kRNA/hmmNUCLEIC, or kAmino/hmmAMINO */
	  String name;                  /* name of alignment, or NULL */
	  String desc;                  /* description of alignment, or NULL */
	  String acc;                   /* accession of alignment, or NULL */
	  public String au;                    /* "author" information, or NULL */
	  String ss_cons;               /* consensus secondary structure string, or NULL */
	  String sa_cons;               /* consensus surface accessibility string, or NULL */
	  public StringBuffer rf;                    /* reference coordinate system, or NULL */
	  ArrayList<String> sqacc = new ArrayList<String>();                 /* accession numbers for individual sequences */
	  ArrayList<String> sqdesc = new ArrayList<String>();                /* description lines for individual sequences */
	  public ArrayList<String> ss = new ArrayList<String>();                    /* per-seq secondary structure annotation, or NULL */
	  public ArrayList<String> sa = new ArrayList<String>();                    /* per-seq surface accessibility annotation, or NULL */
	  double  cutoff[];       /* NC, TC, GA cutoffs propagated to Pfam/Rfam */
	  boolean cutoff_is_set[];/* TRUE if a cutoff is set; else FALSE */

	  /* Optional information that we don't understand.
	   * That is, we know what type of information it is, but it's
	   * either (interpreted as) free-text comment, or it's Stockholm
	   * markup with unfamiliar tags.
	   */
	  ArrayList<String> comment = new ArrayList<String>();              /* free text comments, or NULL      */
	  int     ncomment;             /* number of comment lines          */
	  int     alloc_ncomment;       /* number of comment lines alloc'ed */

	  ArrayList<String> gf_tag = new ArrayList<String>();               /* markup tags for unparsed #=GF lines  */
	  ArrayList<String> gf = new ArrayList<String>();                   /* annotations for unparsed #=GF lines  */
	  int     ngf;                  /* number of unparsed #=GF lines        */
	  int     alloc_ngf;            /* number of gf lines alloc'ed          */

	  ArrayList<String> gs_tag = new ArrayList<String>();               /* markup tags for unparsed #=GS lines     */
	  ArrayList<ArrayList<String>> gs = new ArrayList<ArrayList<String>>();                   /* [0..ngs-1][0..nseq-1][free text] markup */
	  Map<String, Integer> gs_idx = new HashMap<String, Integer>();               /* hash of #=GS tag types                  */
	  int     ngs;                  /* number of #=GS tag types                */

	  ArrayList<String> gc_tag = new ArrayList<String>();               /* markup tags for unparsed #=GC lines  */
	  ArrayList<ArrayList<String>>  gc = new ArrayList<ArrayList<String>>();                   /* [0..ngc-1][0..alen-1] markup         */
	  Map<String, Integer> gc_idx = new HashMap<String, Integer>();               /* hash of #=GC tag types               */
	  int     ngc;                  /* number of #=GC tag types             */

	  ArrayList<String>  gr_tag = new ArrayList<String>();               /* markup tags for unparsed #=GR lines   */
	  ArrayList<ArrayList<String>>  gr = new ArrayList<ArrayList<String>>();                   /* [0..ngr][0..nseq-1][0..alen-1] markup */
	  Map<String, Integer> gr_idx = new HashMap<String, Integer>();               /* hash of #=GR tag types                */
	  int     ngr;                  /* number of #=GR tag types              */
	  public double scores[];
	  
	  /* Stuff we need for our own maintenance of the data structure
	   */
	  Map<String, Integer> index = new HashMap<String, Integer>();                 /* name ->seqidx hash table */
	  int   nseqalloc;             /* number of seqs currently allocated for   */
	  int   nseqlump;              /* lump size for dynamic expansions of nseq */
	  int   sqlen;                 /* individual sequence lengths during parsing */
	  int   sslen;                 /* individual ss lengths during parsing       */
	  int   salen;                 /* individual sa lengths during parsing       */
	  int   lastidx;               /* last index we saw; use for guessing next   */
	  
	  static final int MSA_CUTOFF_TC1 = 0;
	  static final int MSA_CUTOFF_TC2 = 1;
	  static final int MSA_CUTOFF_GA1 = 2;
	  static final int MSA_CUTOFF_GA2 = 3;
	  static final int MSA_CUTOFF_NC1 = 4;
	  static final int MSA_CUTOFF_NC2 = 5;
	  
	  /**  Add an unparsed #=GF markup line to the MSA
	   *           structure, allocating as necessary.
	   *
	   * @param    tag   markup tag (e.g. "AU")
	   * @param    value free text markup (e.g. "Alex Bateman")
	   * 
	   */
	  void MSAAddGF(String tag, String value)
	  {
	    /* If this is our first recorded unparsed #=GF line, we need to malloc();
	     * if we've filled availabl space If we already have a hash index, and the GF
	     * Note the arbitrary lumpsize of 10 lines per allocation...
	     */
	    if (gf_tag == null) {
	    	gf_tag = new ArrayList<String>();
	    	gf = new ArrayList<String>();
	    }
	    gf_tag.add(tag);
	    gf.add(value);

	    return;
	  }
	  /**  From a sequence name, return seqidx appropriate
	   *           for an MSA structure.
	   *
	   *           1) try to guess the index. (pass -1 if you can't guess)
	   *           2) Look up name in msa's hashtable.
	   *           3) If it's a new name, store in msa's hashtable;
	   *                                  expand allocs as needed;
	   *                                  save sqname.
	   *
	   * @param	name   a sequence name
	   * @param guess  a guess at the right index, or -1 if no guess.
	   *
	   * @return  seqidx
	   */
	  int MSAGetSeqidx(String name, int guess)
	  {
	                                  /* can we guess? */
	    if (guess >= 0 && guess < nseq && name.equals(sqname.get(guess)))
	      return guess;
	                                  /* else, a lookup in the index */
	    if (index.containsKey(name))
	      return index.get(name);
	                                  /* else, it's a new name */
	    index.put(name, index.size());
	    MSAExpand(name);
	    nseq++;
	    return index.size()-1;
	  }
	  
	  void MSAExpand(String name)
	  {
		  aseq.add(new StringBuffer());
		  sqname.add(name);
		  wgt.add(-1.0);
		  ss.add("");
		  sa.add("");
		  sqacc.add("");
		  sqdesc.add("");
	  }
	  /** Add an (unparsed) comment line to the MSA structure,
	   *           allocating as necessary.
	   *
	   * @param s  comment line to add
	   *
	   */
	  void
	  MSAAddComment(String s)
	  {
	    /* If this is our first recorded comment, we need to malloc();
	     * and if we've filled available space, we need to realloc().
	     * Note the arbitrary lumpsize of 10 lines per allocation...
	     */
	    comment.add(s);
	    ncomment++;
	    return;
	  }
	  
	  /** Last function called after a multiple alignment is
	   *           parsed. Checks that parse was successful; makes sure
	   *           required information is present; makes sure required
	   *           information is consistent. Some fields that are
	   *           only use during parsing may be freed (sqlen, for
	   *           example).
	   *           
	   *           Some fields in msa may be modified (msa->alen is set,
	   *           for example).
	   *
	   */
	  void MSAVerifyParse() throws InvalidSequenceFormat
	  {
	    int idx;

	    if (nseq == 0) 
	    	throw new InvalidSequenceFormat("Parse error: no sequences were found for alignment %s",
	  			  name != null ? name : "");

	    /* We can rely on msa->sqname[] being valid for any index,
	     * because of the way the line parsers always store any name
	     * they add to the index.
	     */
	    for (idx = 0; idx < nseq; idx++)
	      {
	  				/* aseq is required. */
	        if (aseq.get(idx) == null) 
	        	throw new InvalidSequenceFormat("Parse error: No sequence for %s in alignment %s", sqname.get(idx),
	        			name != null ? name : "");
	  				/* either all weights must be set, or none of them */
//	        if ((flags & MSA_SET_WGT)!=0 && wgt.get(idx) == -1.0)
//	        	throw new InvalidSequenceFormat("Parse error: some weights are set, but %s doesn't have one in alignment %s", 
//	        			sqname.get(idx), name != null ? name : "");
	      }
	  	    
	  				/* all aseq must be same length. */
//	        if (msa.sqlen[idx] != alen)
//	  	Die("Parse error: sequence %s: length %d, expected %d in alignment %s",
//	  	    msa->sqname[idx], msa->sqlen[idx], msa->alen,
//	  	    msa->name != NULL ? msa->name : "");
//	  				/* if SS is present, must have length right */
//	        if (msa->ss != NULL && msa->ss[idx] != NULL && msa->sslen[idx] != msa->alen) 
//	  	Die("Parse error: #=GR SS annotation for %s: length %d, expected %d in alignment %s",
//	  	    msa->sqname[idx], msa->sslen[idx], msa->alen,
//	  	    msa->name != NULL ? msa->name : "");
//	  				/* if SA is present, must have length right */
//	        if (msa->sa != NULL && msa->sa[idx] != NULL && msa->salen[idx] != msa->alen) 
//	  	Die("Parse error: #=GR SA annotation for %s: length %d, expected %d in alignment %s",
//	  	    msa->sqname[idx], msa->salen[idx], msa->alen,
//	  	    msa->name != NULL ? msa->name : "");
//	      }

	  			/* if cons SS is present, must have length right */
//	    if (msa->ss_cons != NULL && strlen(msa->ss_cons) != msa->alen) 
//	      Die("Parse error: #=GC SS_cons annotation: length %d, expected %d in alignment %s",
//	  	strlen(msa->ss_cons), msa->alen,
//	  	msa->name != NULL ? msa->name : "");
//
//	  			/* if cons SA is present, must have length right */
//	    if (msa->sa_cons != NULL && strlen(msa->sa_cons) != msa->alen) 
//	      Die("Parse error: #=GC SA_cons annotation: length %d, expected %d in alignment %s",
//	  	strlen(msa->sa_cons), msa->alen,
//	  	msa->name != NULL ? msa->name : "");
//
//	  				/* if RF is present, must have length right */
//	    if (msa->rf != NULL && strlen(msa->rf) != msa->alen) 
//	      Die("Parse error: #=GC RF annotation: length %d, expected %d in alignment %s",
//	  	strlen(msa->rf), msa->alen,
//	  	msa->name != NULL ? msa->name : "");

	  				/* Check that all or no weights are set */
	    if ((flags & MSA_SET_WGT)==0)
	    {
	    	for (int i=0;i<nseq;i++)
	    	{
	    		wgt.set(i, 1.0);
	    	}
	    }
	  				/* Clean up a little from the parser */
//	    if (msa->sqlen != NULL) { free(msa->sqlen); msa->sqlen = NULL; }
//	    if (msa->sslen != NULL) { free(msa->sslen); msa->sslen = NULL; }
//	    if (msa->salen != NULL) { free(msa->salen); msa->salen = NULL; }

	    return;
	  }
	  public MSA()
	  {
		  aseq   = new ArrayList<StringBuffer>();
		  sqname = new ArrayList<String>();
		  wgt = new ArrayList<Double>();
	  }
	  public MSA(int _nseq, int _alen)
	  {
		  alen = _alen;
		  nseq = _nseq;
		  aseq   = new ArrayList<StringBuffer>();
		  sqname = new ArrayList<String>();
		  
		  aseq   = new ArrayList<StringBuffer>();
		  sqname = new ArrayList<String>();
		  wgt = new ArrayList<Double>();

		  for (int i = 0; i < nseq; i++)
		    {
		      sqname.add("");
		      wgt.add(-1.0);
		      aseq.add(new StringBuffer());
		    }
		  scores = new double[nseq];

		  nseq      = 0;
		  nseqalloc = nseq;
		  nseqlump  = nseq;

		  flags   = 0;
		  type    = SQINFO.SEQTYPE.kOtherSeq;

	  }
}
