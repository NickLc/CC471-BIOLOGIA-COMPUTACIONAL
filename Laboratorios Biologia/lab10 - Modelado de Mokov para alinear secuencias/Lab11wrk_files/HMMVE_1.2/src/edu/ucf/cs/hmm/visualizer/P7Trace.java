package edu.ucf.cs.hmm.visualizer;

import java.util.ArrayList;

import edu.ucf.cs.hmm.squid.MSA;
import edu.ucf.cs.hmm.squid.Sequence;
import edu.ucf.cs.hmm.squid.Util;

/**
 * Adapted from Hmmer
 */
public class P7Trace {
	/*
	 * Plan 7 model state types used in traceback structure
	 */
	public static final int STBOGUS = 0;

	public static final int STM = 1;

	public static final int STD = 2;

	public static final int STI = 3;

	public static final int STS = 4;

	public static final int STN = 5;

	public static final int STB = 6;

	public static final int STE = 7;

	public static final int STC = 8;

	public static final int STT = 9;

	public static final int STJ = 10;

	int tlen; /* length of traceback */

	public int[] statetype; /* state type used for alignment */

	public int[] nodeidx; /* index of aligned node, 1..M (if M,D,I), or 0 */

	public int[] pos; /* position in dsq, 1..L, or 0 if none */

	public double score;

	/**
	 * Allocation and freeing of traceback structures
	 */
	static P7Trace P7AllocTrace(int tlen) {
		P7Trace tr = new P7Trace();
		tr.statetype = new int[tlen];
		tr.nodeidx = new int[tlen];
		tr.pos = new int[tlen];
		return tr;
	}

	/**
	 * Resize traceback structure
	 */
	void P7ReallocTrace(int tlen) {
		int[] newStatetype = new int[tlen];
		int[] newNodeidx = new int[tlen];
		int[] newPos = new int[tlen];
		for (int i = 0; i < statetype.length; i++) {
			newStatetype[i] = statetype[i];
			newNodeidx[i] = nodeidx[i];
			newPos[i] = pos[i];
		}
		statetype = newStatetype;
		nodeidx = newNodeidx;
		pos = newPos;
	}

	/**
	 * Reverse the arrays in a traceback structure. Tracebacks from
	 * Forward() and Viterbi() are collected backwards, and call this function
	 * when they're done.
	 * 
	 * It's possible to reverse the arrays in place more efficiently; but the
	 * realloc/copy strategy has the advantage of reallocating the trace into
	 * the right size of memory. (Tracebacks overallocate.)
	 * 
	 */
	void P7ReverseTrace() {
		int[] statetype_n;
		int[] nodeidx_n;
		int[] pos_n;
		int opos, npos;

		/*
		 * Allocate
		 */
		statetype_n = new int[tlen];
		nodeidx_n = new int[tlen];
		pos_n = new int[tlen];

		/*
		 * Reverse the trace.
		 */
		for (opos = tlen - 1, npos = 0; npos < tlen; npos++, opos--) {
			statetype_n[npos] = statetype[opos];
			nodeidx_n[npos] = nodeidx[opos];
			pos_n[npos] = pos[opos];
		}

		/*
		 * Swap old, new arrays.
		 */
		statetype = statetype_n;
		nodeidx = nodeidx_n;
		pos = pos_n;
	}

	/**
	 * Given a gap-containing string of length n, pull all the non-gap
	 * characters as far as possible to the right, leaving gaps on the left
	 * side. Used to rearrange the positions of insertions in HMMER alignments.
	 */
	static void rightjustify(StringBuffer s, int start, int n) {
		int npos;
		int opos;

		npos = start + n - 1;
		opos = start + n - 1;
		while (opos >= start) {
			if (Util.isgap(s.charAt(opos)))
				opos--;
			else
				s.setCharAt(npos--, s.charAt(opos--));
		}
		while (npos >= start)
			s.setCharAt(npos--, '.');
	}

	/**
	 * Convert an array of traceback structures for a set of sequences
	 * into a new multiple alignment.
	 * 
	 * Insertions are put into lower case and are not aligned; instead, Nterm is
	 * right-justified, Cterm is left-justified, and internal insertions are
	 * split in half and the halves are justified in each direction (the
	 * objective being to increase the chances of getting insertions aligned
	 * well enough for them to become a match). SAM gap char conventions are
	 * used: - in match columns, . in insert columns
	 * 
	 * NOTE: Does not recognize J state.
	 * 
	 * Can handle traces with D->I and I->D transitions; though these are
	 * disallowed by Plan7, they might be generated by aligning an alignment to
	 * a model, in the ImposeMasterTrace() step. Thus, --withali might generate
	 * alignments that are inconsistent with Plan7, that would have to be
	 * trace_doctor()'ed. xref STL6/p.117.
	 * 
	 * @param dsq digitized unaligned sequences 
	 * @param sequences array of sequences
	 * @param mlen length of model (number of match states)
	 * @param tr array of tracebacks
	 * @param alphabet alphabet of the hmm
	 * @param matchonly TRUE if we don't print insert-generated symbols at all
	 * @param MSA structure; null on failure
	 */
	static MSA P7Traces2Alignment(int[][] dsq, ArrayList<Sequence> sequences,
			int mlen, P7Trace[] tr, String alphabet, boolean matchonly) {
		int nseq = sequences.size();
		// SQINFO[] sqinfo
		int idx; /* counter for sequences */
		int alen; /* width of alignment */
		int[] inserts; /* array of max gaps between aligned columns */
		int[] matmap; /* matmap[k] = apos of match k [1..M] */
		int nins; /* counter for inserts */
		int apos; /* position in aligned sequence (0..alen-1) */
		int rpos; /* position in raw digital sequence (1..L) */
		int tpos; /* position counter in traceback */
		int statetype; /* type of current state, e.g. STM */
		int k; /* counter over states in model */

		/*
		 * Here's the problem. We want to align the match states in columns, but
		 * some sequences have inserted symbols in them; we need some sort of
		 * overall knowledge of where the inserts are and how long they are in
		 * order to create the alignment.
		 * 
		 * Here's our trick. inserts[] is a 0..hmm->M array; inserts[i] stores
		 * the maximum number of times insert substate i was used. This is the
		 * maximum number of gaps to insert between canonical column i and i+1.
		 * inserts[0] is the N-term tail; inserts[M] is the C-term tail.
		 * 
		 * Remember that N and C emit on transition, hence the check for an N->N
		 * or C->C transition before bumping nins.
		 */
		inserts = new int[mlen + 1];
		for (k = 0; k <= mlen; k++)
			inserts[k] = 0;
		for (idx = 0; idx < nseq; idx++) {
			nins = 0;
			for (tpos = 0; tpos < tr[idx].tlen; tpos++) {
				switch (tr[idx].statetype[tpos]) {
				case STI:
					nins++;
					break;
				case STN:
					if (tr[idx].statetype[tpos - 1] == STN)
						nins++;
					break;
				case STC:
					if (tr[idx].statetype[tpos - 1] == STC)
						nins++;
					break;
				case STM:
				case STD: /* M,D: record max. reset ctr. */
					if (nins > inserts[tr[idx].nodeidx[tpos] - 1])
						inserts[tr[idx].nodeidx[tpos] - 1] = nins;
					nins = 0;
					break;
				case STB: /* B; record N-tail max, reset ctr */
					if (nins > inserts[0])
						inserts[0] = nins;
					nins = 0;
					break;
				case STT: /* T: record C-tail max */
					if (nins > inserts[mlen])
						inserts[mlen] = nins;
					break;
				case STS:
				case STE:
					break; /* ignore other states */
				case STJ:
					throw new RuntimeException(
							"yo! you don't support J in Traces2Alignment(), remember?");
				default:
					throw new RuntimeException(
							"Traces2Alignment reports unrecognized statetype");
				}
			}
		}

		/* Insert compression option. */
		if (matchonly)
			for (k = 0; k <= mlen; k++)
				if (inserts[k] > 1)
					inserts[k] = 1;

		/***********************************************************************
		 * Construct the alignment
		 **********************************************************************/
		/* calculate alignment length and matmap */
		matmap = new int[mlen + 1];
		matmap[0] = -1;
		alen = inserts[0];
		for (k = 1; k <= mlen; k++) {
			matmap[k] = alen;
			alen += inserts[k] + 1;
		}

		MSA msa = new MSA(nseq, alen); /* RETURN: new alignment */

		for (idx = 0; idx < nseq; idx++) {
			msa.scores[idx] = tr[idx].score;
			/* blank an aseq */
			for (int i = 0; i < alen; i++)
				msa.aseq.get(idx).append(" ");
			for (apos = 0; apos < alen; apos++)
				msa.aseq.get(idx).setCharAt(apos, '.');
			for (k = 1; k <= mlen; k++)
				msa.aseq.get(idx).setCharAt(matmap[k], '-');

			/* align the sequence */
			apos = 0;
			for (tpos = 0; tpos < tr[idx].tlen; tpos++) {
				statetype = tr[idx].statetype[tpos]; /* just for clarity */
				rpos = tr[idx].pos[tpos];
				k = tr[idx].nodeidx[tpos];

				if (statetype == STM) {
					apos = matmap[k];
					msa.aseq.get(idx).setCharAt(apos,
							alphabet.charAt(dsq[idx][rpos]));
					apos++;
				} else if (statetype == STD) {
					apos = matmap[k] + 1; /*
											 * need for handling D->I; xref
											 * STL6/p.117
											 */
				} else if (statetype == STI) {
					if (matchonly)
						msa.aseq.get(idx).setCharAt(apos, '*'); /*
																 * insert
																 * compression
																 * option
																 */
					else {
						msa.aseq.get(idx).setCharAt(apos,
								alphabet.toLowerCase().charAt(dsq[idx][rpos]));
						apos++;
					}
				} else if ((statetype == STN || statetype == STC) && rpos > 0) {
					if (matchonly)
						msa.aseq.get(idx).setCharAt(apos, '*'); /*
																 * insert
																 * compression
																 * option
																 */
					else {
						msa.aseq.get(idx).setCharAt(apos,
								alphabet.toLowerCase().charAt(dsq[idx][rpos]));
						apos++;
					}
				} else if (statetype == STE)
					apos = matmap[mlen] + 1; /* set position for C-term tail */
			}

			/*
			 * N-terminal extension is right-justified. Internal inserts are
			 * split in half, and C-term is right-justified. C-terminal
			 * extension remains left-justified.
			 */
			if (!matchonly) {
				rightjustify(msa.aseq.get(idx), 0, inserts[0]);

				for (k = 1; k < mlen; k++)
					if (inserts[k] > 1) {
						for (nins = 0, apos = matmap[k] + 1; Util
								.islower(msa.aseq.get(idx).charAt(apos)); apos++)
							nins++;
						nins /= 2; /* split the insertion in half */
						rightjustify(msa.aseq.get(idx), matmap[k] + 1 + nins,
								inserts[k] - nins);
					}
			}
		}

		/***********************************************************************
		 * Build the rest of the MSA annotation.
		 **********************************************************************/

		msa.nseq = nseq;
		msa.alen = alen;
		msa.au = "";
		msa.au = String.format("HMMER %s", "2.3.2");
		/* copy sqinfo array and weights */
		for (idx = 0; idx < nseq; idx++) {
			msa.sqname.set(idx, sequences.get(idx).sqInfo.name);
			// if ((sqinfo[idx].flags & SQINFO.SQINFO_ACC)!=0)
			// MSASetSeqAccession(msa, idx, sqinfo[idx].acc);
			// if ((sqinfo[idx].flags & SQINFO.SQINFO_DESC)!=0)
			// MSASetSeqDescription(msa, idx, sqinfo[idx].desc);
			//	
			// if ((sqinfo[idx].flags & SQINFO.SQINFO_SS)!=0) {
			// MakeAlignedString(msa.aseq.get(idx), alen,
			// sqinfo[idx].ss, msa.ss.get(idx)));
			// }
			// if ((sqinfo[idx].flags & SQINFO.SQINFO_SA)!=0) {
			// MakeAlignedString(msa.aseq.get(idx), alen,
			// sqinfo[idx].sa, msa.sa.get(idx)));
			// }
			msa.wgt.set(idx, sequences.get(idx).weight);
		}

		/*
		 * #=RF annotation: x for match column, . for insert column
		 */
		msa.rf = new StringBuffer();
		for (apos = 0; apos < alen; apos++)
			msa.rf.append('.');
		for (k = 1; k <= mlen; k++)
			msa.rf.setCharAt(matmap[k], 'x');

		/*
		 * Currently, we produce no consensus structure. #=CS, generated from
		 * HMM structural annotation, would go here.
		 */
		return msa;
	}
}