package edu.ucf.cs.hmm.squid;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.Calendar;
import java.util.Date;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MSFReader extends MSAReader {
	/* Function: ReadMSF()
	 * Date:     SRE, Tue Jun  1 08:07:22 1999 [St. Louis]
	 *
	 * Purpose:  Parse an alignment read from an open MSF format
	 *           alignment file. (MSF is a single-alignment format.)
	 *           Return the alignment, or NULL if we've already
	 *           read the alignment.
	 *           
	 * Args:     afp  - open alignment file
	 *
	 * Returns:  MSA * - an alignment object
	 *                   caller responsible for an MSAFree()
	 *           NULL if no more alignments
	 *
	 * Diagnostics: 
	 *           Will Die() here with a (potentially) useful message
	 *           if a parsing error occurs.
	 */
	@Override
	public MSA Read(MSAFILE afp) throws IOException {
		MSA msa = new MSA();
		String s;
		int alleged_alen;
		SQINFO.SEQTYPE alleged_type;
		int alleged_checksum;
		int slen;
		int sqidx;

		if ((s = afp.MSAFileGetLine()) == null)
			return null;

		/* The first line is the header.
		 * This is a new-ish GCG feature. Don't count on it, so
		 * we can be a bit more tolerant towards non-GCG software
		 * generating "MSF" files.
		 */
		if (s.startsWith("!!AA_MULTIPLE_ALIGNMENT")) {
			msa.type = SQINFO.SEQTYPE.kAmino;
			if ((s = afp.MSAFileGetLine()) == null)
				return null;
		} else if (s.startsWith("!!NA_MULTIPLE_ALIGNMENT")) {
			msa.type = SQINFO.SEQTYPE.kRNA;
			if ((s = afp.MSAFileGetLine()) == null)
				return null;
		}

		/* Now we're in the free text comment section of the MSF file.
		 * It ends when we see the "MSF: Type: Check: .." line.
		 * This line must be present. 
		 */
		do {
			if (s.contains("..") && s.contains("MSF:")) {
				Pattern p = Pattern
						.compile("^.+MSF: +([0-9]+) +Type: +([PNX]).+Check: +([0-9]+) +\\.\\.");
				Matcher m = p.matcher(s);
				if (m.find()) {
					alleged_alen = Integer.parseInt(m.group(1));
					String t = m.group(2);
					if (t.equals("N")) {
						alleged_type = SQINFO.SEQTYPE.kRNA;
					} else if (t.equals("P")) {
						alleged_type = SQINFO.SEQTYPE.kAmino;
					} else if (t.equals("X")) {
						alleged_type = SQINFO.SEQTYPE.kOtherSeq;
					} else
						alleged_type = SQINFO.SEQTYPE.kOtherSeq;
					alleged_checksum = Integer.parseInt(m.group(3));
					if (msa.type == SQINFO.SEQTYPE.kOtherSeq)
						msa.type = alleged_type;
					break;
				}
			}
			if (s != null && !s.equals(""))
				msa.MSAAddComment(s);
		} while ((s = afp.MSAFileGetLine()) != null);

		/* Now we're in the name section.
		 * GCG has a relatively poorly documented feature: only sequences that
		 * appear in this list will be read from the alignment section. Commenting
		 * out sequences in the name list (by preceding them with "!") is
		 * allowed as a means of manually defining subsets of sequences in
		 * the alignment section. We can support this feature reasonably
		 * easily because of the hash table for names in the MSA: we
		 * only add names to the hash table when we see 'em in the name section.
		 */
		while ((s = afp.MSAFileGetLine()) != null) {
			s=s.trim();

			if (s.equals(""))
				continue; /* skip blank lines */
			else if (s.startsWith("!"))
				msa.MSAAddComment(s);
			else if (s.contains("Name:")) {
				/* We take the name and the weigh, and that's it */
				Pattern p = Pattern
						.compile("^Name:\\s+?(\\S+).+?Len:\\s+?(\\S+)\\s+Check:\\s+(\\S+)\\s+?Weight:\\s+(\\S+)");
				Matcher m = p.matcher(s);
				if (m.find()) {
					String name;
					double weight;
					name = m.group(1).trim();
					weight = Double.parseDouble(m.group(4));
					msa.index.put(name, msa.index.size());
					sqidx = msa.index.size();
					msa.MSAExpand(name);
					msa.nseq++;
					msa.wgt.add(weight);
					msa.flags |= MSA.MSA_SET_WGT;
				}
			} else if (s.startsWith("//"))
				break;
			else {
				throw new InvalidSequenceFormat();
			}

		}

		/* And now we're in the sequence section. 
		 * As discussed above, if we haven't seen a sequence name, then we
		 * don't include the sequence in the alignment.
		 * Also, watch out for coordinate-only lines.
		 */
		while ((s = afp.MSAFileGetLine()) != null) {
			s = s.trim();
			if (s.equals(""))
				continue;
			if (s.indexOf(" ")==-1)
				continue;
			String name = s.substring(0, s.indexOf(" "));
			String seq = s.substring(s.indexOf(" ")).trim();

			/* The test for a coord line: digits starting both fields
			 */
			try {
				Integer.parseInt(name);
				Integer.parseInt(seq);
				continue;
			} catch (NumberFormatException e) {
			}

			/* It's not blank, and it's not a coord line: must be sequence
			 */
			sqidx = msa.index.get(name);
			if (sqidx < 0)
				continue; /* not a sequence we recognize */

			msa.aseq.set(sqidx, new StringBuffer(msa.aseq.get(sqidx)+seq));
		}

		/* We've left blanks in the aseqs; take them back out.
		 */
		for (sqidx = 0; sqidx < msa.nseq; sqidx++) {
			String sequence = msa.aseq.get(sqidx).toString();
			if (sequence == null || sequence.equals(""))
				throw new InvalidSequenceFormat();

			String newSeq = "";
			for (int i = 0; i < sequence.length(); i++) {
				if (sequence.charAt(i) == ' ' || sequence.charAt(i) == '\t')
					continue;
				newSeq += sequence.charAt(i);
			}
			msa.aseq.set(sqidx, new StringBuffer(newSeq));
		}

		msa.MSAVerifyParse(); /* verifies, and also sets alen and wgt. */
		return msa;
	}
		/* Function: WriteMSF()
		 * Date:     SRE, Mon May 31 11:25:18 1999 [St. Louis]
		 *
		 * Purpose:  Write an alignment in MSF format to an open file.
		 *
		 * Args:     fp    - file that's open for writing.
		 *           msa   - alignment to write. 
		 *
		 *                   Note that msa->type, usually optional, must be
		 *                   set for WriteMSF to work. If it isn't, a fatal
		 *                   error is generated.
		 *
		 * Returns:  (void)
		 */
		void write(Writer w, MSA msa)
		{
		  Date now;			/* current time as Date */
		  StringBuffer gcg_aseq[];              /* aligned sequences with gaps converted to GCG format */
		  StringBuffer gcg_sqname[];		/* sequence names with GCG-valid character sets */
		  int    idx;			/* counter for sequences         */
		  int    len;			/* tmp variable for name lengths */
		  int    namelen;		/* maximum name length used      */
		  int    pos;			/* position counter              */
		  String buffer;		/* buffer for writing seq        */
		  int    i;			/* another position counter */

		  /*****************************************************************
		   * Make copies of sequence names and sequences.
		   *   GCG recommends that name characters should only contain
		   *   alphanumeric characters, -, or _
		   *   Some GCG and GCG-compatible software is sensitive to this.
		   *   We silently convert all other characters to '_'.
		   *   
		   *   For sequences, GCG allows only ~ and . for gaps.
		   *   Otherwise, everthing is interpreted as a residue;
		   *   so squid's IUPAC-restricted chars are fine. ~ means
		   *   an external gap. . means an internal gap.
		   *****************************************************************/ 
		   
						/* make copies that we can edit */
		  PrintWriter writer = new PrintWriter(w);
		   gcg_aseq   = new StringBuffer[msa.nseq];
		   gcg_sqname = new StringBuffer[msa.nseq];
		   for (idx = 0; idx < msa.nseq; idx++)
		     {
		       gcg_aseq[idx]   = new StringBuffer(msa.aseq.get(idx));
		       gcg_sqname[idx] = new StringBuffer(msa.sqname.get(idx));
		     }
						/* alter names as needed  */
		   for (idx = 0; idx < msa.nseq; idx++)
		     for (i=0; i<gcg_sqname[idx].length(); i++)
		       if (! Util.isalnum(gcg_sqname[idx].charAt(i)) && gcg_sqname[idx].charAt(i) != '-' && gcg_sqname[idx].charAt(i) != '_')
		    	   gcg_sqname[idx].setCharAt(i, '_');
						/* alter gap chars in seq  */
		   for (idx = 0; idx < msa.nseq; idx++)
		     {
		       for (i=0;i<gcg_aseq[idx].length()&&Util.isgap(gcg_aseq[idx].charAt(i));i++)
		    	   gcg_aseq[idx].setCharAt(i, '~');
		       for (; i<gcg_aseq[idx].length(); i++)
		    	   if (Util.isgap(gcg_aseq[idx].charAt(i))) gcg_aseq[idx].setCharAt(i, '.');
		       for (pos = msa.alen-1; pos > 0 && Util.isgap(gcg_aseq[idx].charAt(pos)); pos--)
		    	   gcg_aseq[idx].setCharAt(pos, '~');
		     }
						/* calculate max namelen used */
		  namelen = 0;
		  for (idx = 0; idx < msa.nseq; idx++)
		    if ((len = msa.sqname.get(idx).length()) > namelen) 
		      namelen = len;

		  /*****************************************************
		   * Write the MSF header
		   *****************************************************/
						/* required file type line */
		  if (msa.type == SQINFO.SEQTYPE.kOtherSeq)
		    msa.type = SQFILE.GuessAlignmentSeqtype(msa.aseq, msa.nseq);

		  if      (msa.type == SQINFO.SEQTYPE.kRNA)   writer.println("!!NA_MULTIPLE_ALIGNMENT 1.0");
		  else if (msa.type == SQINFO.SEQTYPE.kDNA)   writer.println("!!NA_MULTIPLE_ALIGNMENT 1.0");
		  else if (msa.type == SQINFO.SEQTYPE.kAmino) writer.println("!!AA_MULTIPLE_ALIGNMENT 1.0");
		  else if (msa.type == SQINFO.SEQTYPE.kOtherSeq) 
		    throw new RuntimeException("WriteMSF(): couldn't guess whether that alignment is RNA or protein.%n"); 
		  else    
		    throw new RuntimeException("Invalid sequence type in WriteMSF()"); 

						/* free text comments */
		  if (msa.ncomment > 0)
		    {
		      for (idx = 0; idx < msa.ncomment; idx++)
		    	  writer.println(msa.comment.get(idx));
		      writer.println();
		    }
						/* required checksum line */
		  now = new Date();
		  writer.format(" %s  MSF: %d  Type: %c  %s  Check: %d  ..%n", 
			  msa.name != null ? msa.name : "squid.msf",
			  msa.alen,
			  msa.type == SQINFO.SEQTYPE.kRNA ? 'N' : 'P',
			  String.format("%tB %td, %tY %tH:%tM", now,now,now,now,now),
			  SQFILE.GCGMultchecksum(gcg_aseq, msa.nseq));
		  writer.println();

		  /*****************************************************
		   * Names/weights section
		   *****************************************************/

		  for (idx = 0; idx < msa.nseq; idx++)
		    {
		      writer.format(" Name: %-"+namelen+"s  Len:  %5d  Check: %4d  Weight: %.2f%n",
			      gcg_sqname[idx],
			      msa.alen,
			      SQFILE.GCGchecksum(gcg_aseq[idx], msa.alen),
			      msa.wgt.get(idx));
		    }
		  writer.println();
		  writer.println("//");

		  /*****************************************************
		   * Write the sequences
		   *****************************************************/

		  for (pos = 0; pos < msa.alen; pos += 50)
		    {
		      writer.println();	/* Blank line between sequence blocks */

						/* Coordinate line */
		      len = (pos + 50) > msa.alen ? msa.alen - pos : 50;
		      if (len > 10)
			writer.format("%"+namelen+"s  %-6d%"+(len + ((len-1)/10) - 12)+"s%6d%n", "", 
				pos+1,
				"", pos + len);
		      else
			writer.format("%"+namelen+"s  %-6d%n", "", pos+1);

		      for (idx = 0; idx < msa.nseq; idx++)
			{
			  writer.format("%-"+namelen+"s ", gcg_sqname[idx]);
						/* get next line's worth of 50 from seq */
			  buffer = gcg_aseq[idx].substring(pos, pos+50<gcg_aseq[idx].length()? pos+50:gcg_aseq[idx].length());
			  		/* draw the sequence line */
			  for (i = 0; i < len; i++)
			    {
			      if ((i % 10)==0) writer.print(' ');
			      writer.print(buffer.charAt(i));
			    }
			  writer.println();
			}
		    }
		  writer.close();
		  return;
	}
}
