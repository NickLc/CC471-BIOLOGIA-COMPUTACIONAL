package edu.ucf.cs.hmm.squid;

/**
 * Adapted from Hammer squid library
 *
 * 
 */

public class SQINFO {
	static public final int SQINFO_NAME	=	(1 << 0);
	static public final int SQINFO_ID		=	(1 << 1);
	static public final int SQINFO_ACC   	=	(1 << 2);
	static public final int SQINFO_DESC  	=	(1 << 3);
	static public final int SQINFO_START 	=	(1 << 4);
	static public final int SQINFO_STOP  	=	(1 << 5);
	static public final int SQINFO_LEN   	=	(1 << 6);
	static public final int SQINFO_TYPE  	=	(1 << 7);
	static public final int SQINFO_OLEN  	=	(1 << 8);
	static public final int SQINFO_SS    	=	(1 << 9);
	static public final int SQINFO_SA    	=	(1 << 10);
	
	static public final int SQINFO_NAMELEN = 64;
	static public final int SQINFO_DESCLEN = 128;
	
	enum SEQTYPE {
		kOtherSeq, kDNA, kRNA, kAmino
	}
	
	public int flags; /* what extra data are available         */

	public String name; /* up to 63 characters of name           */

	String id; /* up to 63 char of database identifier  */

	public String acc; /* up to 63 char of database accession # */

	public String desc;/* up to 127 char of description         */

	int len; /* length of this seq                    */

	int start; /* (1..len) start position on source seq */

	int stop; /* (1..len) end position on source seq   */

	int olen; /* original length of source seq         */

	public SEQTYPE type; /* kRNA, kDNA, kAmino, or kOther         */

	public String ss; /* 0..len-1 secondary structure string   */

	public String sa; /* 0..len-1 % side chain surface access. */
	
	int SetSeqinfoString(String sptr, int flag) throws InvalidFlagException
	{
					/* silently ignore NULL. */
	  if (sptr == null) return 1;

	  sptr.trim();	/* ignore whitespace */

	  switch (flag) {
	  case SQINFO.SQINFO_NAME:
	    if (sptr.length()>0 && sptr.charAt(0) != '-')
	    { 
	    	name = sptr;
	    	flags   |= SQINFO.SQINFO_NAME;
	     }
	    break;

	  case SQINFO.SQINFO_ID:
	    if (sptr.length()>0 && sptr.charAt(0) != '-')
	    { 
	    	id = sptr;
	    	flags |= SQINFO.SQINFO_ID;
	    }
	    break;

	  case SQINFO.SQINFO_ACC:
	    if (sptr.length()>0 && sptr.charAt(0) != '-')
	    {
	    	acc = sptr;
	    	flags   |= SQINFO.SQINFO_ACC;
	    }
	    break;

	  case SQINFO.SQINFO_DESC:
	    if (sptr.length()>0 && sptr.charAt(0) != '-')
	    { 
		  desc += sptr;
		  flags   |= SQINFO.SQINFO_DESC;
	    }
	    break;

	  case SQINFO.SQINFO_START:
		  try
		  {
			  start = Integer.parseInt(sptr);
		  }
		  catch (NumberFormatException e)
		  {
			  return 0;
		  }
		  if (start != 0) flags |= SQINFO.SQINFO_START;
	    break;

	  case SQINFO.SQINFO_STOP:
	    	try
	    	{
	    		stop = Integer.parseInt(sptr);
	    	}
	    	catch (NumberFormatException e)
	    	{
	    		return 0;
	    	}
	    	if (stop != 0) flags |= SQINFO.SQINFO_STOP;
	    break;

	  case SQINFO.SQINFO_OLEN:
	    try
	    {
	    	olen = Integer.parseInt(sptr);
	    }
	    catch (NumberFormatException e)
	    {
	    	return 0;
	    }
	    if (olen != 0) flags |= SQINFO.SQINFO_OLEN;
	    break;

	  default:
	    throw new InvalidFlagException();
	  }
	  return 1;
	}
}
