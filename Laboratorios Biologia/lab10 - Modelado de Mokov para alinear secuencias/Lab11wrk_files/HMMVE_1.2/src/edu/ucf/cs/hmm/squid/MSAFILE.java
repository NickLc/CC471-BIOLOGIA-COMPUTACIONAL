package edu.ucf.cs.hmm.squid;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;

public class MSAFILE {
	  BufferedReader reader;                      /* open file pointer                         */
	  String fname;                  /* name of file. used for diagnostic output  */
	  int   linenumber;             /* what line are we on in the file           */

	  String buf;                    /* buffer for line input w/ sre_fgets() */
	  int   buflen;                 /* current allocated length for buf     */

	  SSIFILE ssi;                 /* open SSI index file; or NULL, if none. */

	  int   do_gzip;                /* TRUE if f is a pipe from gzip -dc (need pclose(f))  */
	  int   do_stdin;               /* TRUE if f is stdin (don't close f, not our problem) */
	  SQFILE.SQFILE_FORMAT format;                 /* format of alignment file we're reading */
		/* Function: MSAFileRead()
		 * Date:     SRE, Fri May 28 16:01:43 1999 [St. Louis]
		 *
		 * Purpose:  Read the next msa from an open alignment file.
		 *           This is a wrapper around format-specific calls.
		 *
		 * Args:     afp     - open alignment file
		 *
		 * Returns:  next alignment, or NULL if out of alignments
		 */
		MSA Read() throws IOException
		{
		  MSA msa = null;
		  MSAReader r;
		  
		  if (format == SQFILE.SQFILE_FORMAT.MSAFILE_STOCKHOLM)
		  {
			  r = new StockholmReader();
		  } else if (format == SQFILE.SQFILE_FORMAT.MSAFILE_MSF)
		  {
			  r = new MSFReader();
		  } else
			  throw new NotImplementedException();
//		  switch (format) {
//		  case MSAFILE_STOCKHOLM: msa = ReadStockholm(afp); break;
//		  case MSAFILE_MSF:       msa = ReadMSF(afp);       break;
//		  case MSAFILE_A2M:       msa = ReadA2M(afp);       break;
//		  case MSAFILE_CLUSTAL:   msa = ReadClustal(afp);   break;
//		  case MSAFILE_SELEX:     msa = ReadSELEX(afp);     break;
//		  case MSAFILE_PHYLIP:    msa = ReadPhylip(afp);    break;
//		  default:
//		    throw new IOException();
		  msa = r.Read(this);
		  return msa;
		}
		String MSAFileGetLine() throws IOException
		{
			buf = reader.readLine();
			linenumber++;
			return buf;
		}
		public void write(Writer w, MSA msa, SQFILE.SQFILE_FORMAT outfmt, boolean do_oneline)
		{
//			  if (outfmt==SQFILE.SQFILE_FORMAT.MSAFILE_A2M)
//			  {
//				  WriteA2M(w, msa);
//			  }
//			  else if (outfmt==SQFILE.SQFILE_FORMAT.MSAFILE_CLUSTAL)
//			  {
//				  WriteClustal(w, msa);
//			  }
//			  else if (outfmt==SQFILE.SQFILE_FORMAT.MSAFILE_MSF)
			  if (outfmt==SQFILE.SQFILE_FORMAT.MSAFILE_MSF)
			  {
				  MSFReader r = new MSFReader();
				  r.write(w, msa);
			  }
//			  {      WriteMSF(w, msa);
//			  }
//			  else if (outfmt==SQFILE.SQFILE_FORMAT.MSAFILE_PHYLIP)
//			  {    WritePhylip(w, msa);
//			  }
//			  else if (outfmt==SQFILE.SQFILE_FORMAT.MSAFILE_SELEX)
//			  {
//				  WriteSELEX(w, msa);
//			  }
//			  else if (outfmt==SQFILE.SQFILE_FORMAT.MSAFILE_STOCKHOLM)
//			  {
//			    if (do_oneline) WriteStockholmOneBlock(w, msa);
//			    else            WriteStockholm(w, msa);
//			  }
//			  else
//			    throw new RuntimeException("can't write. no such alignment format");
//			  }
//			}
		}
}
