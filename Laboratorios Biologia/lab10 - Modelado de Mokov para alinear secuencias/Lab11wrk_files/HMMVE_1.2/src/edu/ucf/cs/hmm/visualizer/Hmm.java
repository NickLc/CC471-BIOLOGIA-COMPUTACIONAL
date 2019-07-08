package edu.ucf.cs.hmm.visualizer;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;

import edu.ucf.cs.hmm.squid.MSA;
import edu.ucf.cs.hmm.squid.Sequence;
import java.util.Arrays;

/**
 * 
 * Part of the code is adapted from Hmmer
 *
 */
public class Hmm {
	/*
	 * Flags for plan7->flags. Note: Some models have scores but no
	 * probabilities (for instance, after reading from an HMM save file). Other
	 * models have probabilities but no scores (for instance, during training or
	 * building). Since it costs time to convert either way, I use PLAN7_HASBITS
	 * and PLAN7_HASPROB flags to defer conversion until absolutely necessary.
	 * This means I have to be careful about keeping these flags set properly
	 * when I fiddle a model.
	 */
	public static final int PLAN7_HASBITS = (1 << 0); /*
														 * raised if model has
														 * log-odds scores
														 */

	public static final int PLAN7_DESC = (1 << 1); /*
													 * raised if description
													 * exists
													 */

	public static final int PLAN7_RF = (1 << 2); /*
													 * raised if #RF annotation
													 * available
													 */

	public static final int PLAN7_CS = (1 << 3); /*
													 * raised if #CS annotation
													 * available
													 */

	public static final int PLAN7_XRAY = (1 << 4); /*
													 * raised if structural data
													 * available
													 */

	public static final int PLAN7_HASPROB = (1 << 5); /*
														 * raised if model has
														 * probabilities
														 */

	public static final int PLAN7_HASDNA = (1 << 6); /*
														 * raised if protein
														 * HMM->DNA seq params
														 * set
														 */

	//renamed from PLAN7_STATS to PLAN7_EVD
	// since hmmer 3 has a field named 'STATS'
	public static final int PLAN7_EVD = (1 << 7); /*
													 * raised if EVD parameters
													 * are available
													 */

	public static final int PLAN7_MAP = (1 << 8); /*
													 * raised if alignment map
													 * is available
													 */

	public static final int PLAN7_ACC = (1 << 9); /*
													 * raised if accession
													 * number is available
													 */

	public static final int PLAN7_GA = (1 << 10); /*
													 * raised if gathering
													 * thresholds available
													 */

	public static final int PLAN7_TC = (1 << 11); /*
													 * raised if trusted cutoffs
													 * available
													 */

	public static final int PLAN7_NC = (1 << 12); /*
													 * raised if noise cutoffs
													 * available
													 */

	public static final int PLAN7_CA = (1 << 13); /*
													 * raised if surface
													 * accessibility avail.
													 */

	public static final String PACKAGE_VERSION = "2.3.2";

	public static final int INTSCALE = 1000;

	public static final int MAXABET = 20;

	public static final int MAXCODE = 24;

	public static final int INFTY = 987654321;

	/*
	 * Indices for Plan7 main model state transitions. Used for indexing
	 * hmm->t[k][] mnemonic: Transition from Match to Match = TMM
	 */
	public static final int TMM = 0;

	public static final int TMI = 1;

	public static final int TMD = 2;

	public static final int TIM = 3;

	public static final int TII = 4;

	public static final int TDM = 5;

	public static final int TDD = 6;

	/*
	 * Indices for special state types, I: used for dynamic programming xmx[][]
	 * mnemonic: eXtra Matrix for B state = XMB
	 */
	public static final int XMB = 0;

	public static final int XME = 1;

	public static final int XMC = 2;

	public static final int XMJ = 3;

	public static final int XMN = 4;

	/*
	 * Indices for special state types, II: used for hmm->xt[] indexing
	 * mnemonic: eXtra Transition for N state = XTN
	 */
	public static final int XTN = 0;

	public static final int XTE = 1;

	public static final int XTC = 2;

	public static final int XTJ = 3;

	/*
	 * Indices for extra state transitions Used for indexing hmm->xt[][].
	 */
	public static final int MOVE = 0; /* trNB, trEC, trCT, trJB */

	public static final int LOOP = 1; /* trNN, trEJ, trCC, trJJ */

	private enum Atype {
		hmmAMINO, hmmNUCLEIC
	}

	Atype alphabetType;

	int alphabetSize = 0;

	int alphabetIupac = 0;

	char Degenerate[][] = new char[MAXCODE][MAXABET];

	int DegenCount[] = new int[MAXCODE];

	/*
	 * Annotation on the model. A name is mandatory. Other fields are optional;
	 * whether they are present is flagged in the stateflags bit array.
	 */
	public String name; // name of the model

	public String acc; /* accession number of model (Pfam) + */

	String desc; /* brief description of model + */

	String[] rf; /* reference line from alignment 0..M + */

	String[] cs; /* consensus structure line 0..M + */

	String[] ca; /* consensus accessibility line 0..M */

	ArrayList<String> comlog = new ArrayList<String>(); /*
														 * command line(s) that
														 * built model +
														 */

	public int nseq; /* number of training sequences + */

	String ctime; /* creation date + */

	int map[]; // map of alignment cols onto model 1..M

	long checksum; /* checksum of training sequences + */

	double ga1, ga2; /* per-seq/per-domain gathering thresholds (bits) + */

	double tc1, tc2; /* per-seq/per-domain trusted cutoff (bits) + */

	double nc1, nc2; /* per-seq/per-domain noise cutoff (bits) + */

	String alphabets[];

	String alphabet;

	/*
	 * The main model in probability form: data-dependent probabilities. This is
	 * the core Krogh/Haussler model. Transition probabilities are usually
	 * accessed as a two-D array: hmm->t[k][TMM], for instance. They are
	 * allocated such that they can also be stepped through in 1D by pointer
	 * manipulations, for efficiency in DP algorithms.
	 */
	public int M; // length of the model (# nodes)

	public double t[][]; // transition prob's. t[1..M-1][0..6]

	public double mat[][]; // match emissions. mat[1..M][0..alphabet-1]

	double ins[][]; // insert emissions. ins[1..M-1][0..alphabet-1]

	double tbd1; // B->D1 prob (data dependent)

	/*
	 * The unique states of Plan 7 in probability form. These are the
	 * algorithm-dependent, data-independent probabilities. Some parts of the
	 * code may briefly use a trick of copying tbd1 into begin[0]; this makes it
	 * easy to call FChoose() or FNorm() on the resulting vector. However, in
	 * general begin[0] is not a valid number.
	 */
	double xt[][] = new double[4][2]; // N,E,C,J extra states: 2 transitions

	public double begin[]; // 1..M B->M state transitions

	double end[]; // 1..M M->E state transitions (!= a dist!)

	/*
	 * The null model probabilities.
	 */
	public double nullModel[] = new double[MAXABET]; // random sequence" emission
												// prob's

	double p1; // null model loop probability

	/*
	 * The model in log-odds score form. These are created from the
	 * probabilities by LogoddsifyHMM(). By definition, null[] emission scores
	 * are all zero. Note that emission distributions are over 26 upper-case
	 * letters, not just the unambiguous protein or DNA alphabet: we
	 * precalculate the scores for all IUPAC degenerate symbols we may see.
	 * Non-IUPAC symbols simply have a -INFTY score.
	 * 
	 * Note the reversed indexing on msc, isc, tsc -- for efficiency reasons.
	 * They're not probability vectors any more so we can reorder them without
	 * wildly complicating our life.
	 * 
	 * The _mem ptrs are where the real memory is alloc'ed and free'd, as
	 * opposed to where it is accessed. This came in with Erik Lindahl's altivec
	 * port; it allows alignment on 16-byte boundaries. In the non-altivec code,
	 * this is just a little redundancy; tsc and tsc_mem point to the same
	 * thing, for example.
	 * 
	 * Only valid if PLAN7_HASBITS is set.
	 */
	int tsc[][]; /* transition scores [0.6][1.M-1] - */

	int msc[][]; /* match emission scores [0.MAXCODE-1][1.M] - */

	int isc[][]; /* ins emission scores [0.MAXCODE-1][1.M-1] - */

	int xsc[][] = new int[4][2]; /* N,E,C,J transitions - */

	int bsc[]; /* begin transitions [1.M] - */

	int esc[]; /* end transitions [1.M] - */

	// int *tsc_mem, *msc_mem, *isc_mem, *bsc_mem, *esc_mem;
	/*
	 * DNA translation scoring parameters For aligning protein Plan7 models to
	 * DNA sequence. Lookup value for a codon is calculated by pos1 * 16 + pos2 *
	 * 4 + pos3, where 'pos1' is the digitized value of the first nucleotide
	 * position; if any of the positions are ambiguous codes, lookup value 64 is
	 * used (which will generally have a score of zero)
	 * 
	 * Only valid if PLAN7_HASDNA is set.
	 */
	int dnam[][]; /* triplet match scores [0.64][1.M] - */

	int dnai[][]; /* triplet insert scores [0.64][1.M] - */

	int dna2; /* -1 frameshift, doublet emission, M or I - */

	int dna4; /* +1 frameshift, doublet emission, M or I - */

	/*
	 * P-value and E-value statistical parameters Only valid if PLAN7_EVD is
	 * set.
	 */
	double mu; /* EVD mu + */

	double lambda; /* EVD lambda + */

	long flags; /* bit flags indicating state of HMM, valid data + */

	public String consensus;

	double consensusP[];

	public boolean modified = false;



	// HMMER3 specific variables
	static final int PLAN7_EFFN = (1 << 14); /*  */
	static final int PLAN7_STATS = (1 << 15); /*  */
	static final int PLAN7_COMPO = (1 << 16); /*  */
	int hmmerVersion;

	double effn; /* effective number of training sequences + */

	double[][] stats=new double[3][2];  /* STATS LOCAL values. only valid if PLAN7_STATS is set */

	double[] compo;
	/////////// end HMMER3 variables ////

	/**
	 *  convert ascii string to probability. Assumes string
	 * is a negagive natural log. Used for HMMER3
	 * @param s
	 * @return
	 */
	static double asciiNegLn2Prob(String s){
		if(s.equals("*")){
			return 0;
		}
		return negLn2Prob(Double.parseDouble(s));
	}
	/**
	 * convert from a negative natrual log to a probability
	 * @param d
	 * @return
	 */
	static double negLn2Prob(double d)
	{
		return  Math.exp(-d);
	}
	/**
	 *  format a probabily for output in a HMMER3 file.
	 * 0 probabilyt is represented as '*'. 5 decimals are printed.
	 * @param d
	 * @return
	 */
	static String prob2NegLnAscii(double d)
	{
		if(d==0)
			return "*";
		double logd=prob2NegLn(d);
		if(logd == 0 )
			logd= Math.abs(logd); // avoid negative zero, messes up file format
		return String.format("%.5f",logd);
	}
	/**
	 * convert from a proabilty to a negative nattural log.
	 * @param d
	 * @return.
	 */
	static double prob2NegLn(double d)
	{
		return -1 * Math.log(d);
	}

/**
	 * Conversion function from ascii to probability
	 *
	 * @param s ascii string
	 * @param n odds in null model
	 */
	static double ascii2prob(String s, double n) {
		if (s.equals("*")) {
			return 0;
		}
		return score2Prob(Integer.parseInt(s), n);
	}

	/**
	 * Format a probability for output to an ASCII save file. 
	 * 
	 */
	static String prob2ascii(double p, double nullModel) {
		String buffer = "";

		if (p == 0.0)
			return "*";
		buffer = String.format("%6d", Prob2Score(p, nullModel));
		return buffer;
	}

	/**
	 * Convert a probability to a scaled integer log_2 odds score.
	 * Round to nearest integer (i.e. note use of +0.5 and floor()) Return the
	 * score.
	 */
	static int Prob2Score(double p, double nullModel) {
		if (p == 0.0)
			return -INFTY;
		else
			return (int) Math.floor(0.5 + INTSCALE * sreLOG2(p / nullModel));
	}

	/**
	 * Math function to convert score to probability
	 * @param sc score
	 * @param n odds in null model
	 * @return probability
	 */
	public static double score2Prob(int sc, double n) {
		return (n * sreEXP2((double) sc / INTSCALE));
	}

	/**
	 * Math function to calculate exp2
	 */
	static double sreEXP2(double sc) {
		return Math.exp((sc) * 0.69314718);
	}

	/**
	 * Convert a probability to a scaled integer log_2 odds score.
	 * Round to nearest integer (i.e. note use of +0.5 and floor()) Return the
	 * score.
	 */
	static int prob2Score(double p, double n) {
		if (p == 0.0)
			return -INFTY;
		else
			return (int) Math.floor(0.5 + INTSCALE * sreLOG2(p / n));
	}

	/**
	 * Math function to calculate log2
	 */
	static double sreLOG2(double sc) {
		return ((sc) > 0 ? Math.log(sc) * 1.44269504 : -9999.);
	}

	/**
	 * Returns the log of the sum of two log probabilities.
	 * log(exp(p1)+exp(p2)) = p1 + log(1 + exp(p2-p1)) for p1 > p2 Note that
	 * this is in natural log space, not log_2.
	 */
	double LogSum(double p1, double p2) {
		if (p1 > p2)
			return (p1 - p2 > 50.) ? p1 : p1 + Math.log(1. + Math.exp(p2 - p1));
		else
			return (p2 - p1 > 50.) ? p2 : p2 + Math.log(1. + Math.exp(p1 - p2));
	}

	/**
	 * Convenience function for setting up Degenerate[][] global for
	 * the alphabet.
	 */
	void set_degenerate(char iupac, String syms) {
		DegenCount[alphabet.indexOf(iupac)] = syms.length();
		for (int i = 0; i < syms.length(); i++) {
			Degenerate[alphabet.indexOf(iupac)][alphabet
					.indexOf(syms.charAt(i))] = 1;
		}
	}

	/**
	 * Set the alphabet globals, given an alphabet type of either
	 * hmmAMINO or hmmNUCLEIC.
	 */
	void setAlphabet(Atype type) {
		int x;

		if (type == Atype.hmmAMINO) {
			alphabetType = type;
			alphabet = "ACDEFGHIKLMNPQRSTVWYUBZX";
			alphabetSize = 20;
			alphabetIupac = 24;
			for (x = 0; x < alphabetIupac; x++) {
				for (int i = 0; i < alphabetSize; i++) {
					Degenerate[x][i] = 0;
				}
			}
			for (x = 0; x < alphabetSize; x++) {
				Degenerate[x][x] = 1;
				DegenCount[x] = 1;
			}
			set_degenerate('U', "S"); /* selenocysteine is treated as serine */
			set_degenerate('B', "ND");
			set_degenerate('Z', "QE");
			set_degenerate('X', "ACDEFGHIKLMNPQRSTVWY");
		} else if (type == Atype.hmmNUCLEIC) {
			alphabetType = type;
			alphabet = "ACGTUNRYMKSWHBVDX";
			alphabetSize = 4;
			alphabetIupac = 17;
			for (x = 0; x < alphabetIupac; x++) {
				for (int i = 0; i < alphabetSize; i++)
					Degenerate[x][i] = 0;
			}
			for (x = 0; x < alphabetSize; x++) {
				Degenerate[x][x] = 1;
				DegenCount[x] = 1;
			}
			set_degenerate('U', "T");
			set_degenerate('N', "ACGT");
			set_degenerate('X', "ACGT");
			set_degenerate('R', "AG");
			set_degenerate('Y', "CT");
			set_degenerate('M', "AC");
			set_degenerate('K', "GT");
			set_degenerate('S', "CG");
			set_degenerate('W', "AT");
			set_degenerate('H', "ACT");
			set_degenerate('B', "CGT");
			set_degenerate('V', "ACG");
			set_degenerate('D', "AGT");
		} else
			throw new RuntimeException(
					"No support for non-nucleic or protein alphabets");
	}

	void initializeArrays()
	{
		// initial arrays
		t = new double[hmmerVersion==3?M+1:M][7];
		mat = new double[M + 1][alphabetSize];
		ins = new double[hmmerVersion==3?M+1:M][alphabetSize];
		begin = new double[M + 1];
		end = new double[M + 1];
		map = new int[M + 1];

		tsc = new int[7][M];
		msc = new int[MAXCODE][M + 1];
		isc = new int[MAXCODE][M];
		bsc = new int[M + 1];
		esc = new int[M + 1];

		dnam = new int[65][M + 1];
		dnai = new int[65][M + 1];

		rf = new String[M+1];
		cs = new String[M+1];
	}

	/**
	 * Read from hmm file
	 * @param hmmFilename
	 * @throws IOException
	 * @throws InvalidFormatException
	 */
	public void read(String hmmFilename) throws IOException,
			InvalidFormatException {
		BufferedReader reader = new BufferedReader(new FileReader(hmmFilename));
		read(reader);
		reader.close();
	}

	/**
	 * Read from hmm file
	 * @param reader hmm reader
	 * @throws IOException
	 * @throws InvalidFormatException
	 */
	public void read(BufferedReader reader) throws IOException, InvalidFormatException 
	{
		String buffer;
		buffer = reader.readLine();
		if(buffer.startsWith("HMMER2.0"))
			hmmerVersion=2;
		else if(buffer.startsWith("HMMER3/b"))
			hmmerVersion=3;
		else
			throw new InvalidFormatException();


		while ((buffer = reader.readLine()) != null) {
			if( ! readCommonHeader(buffer))
				break;

			if(hmmerVersion == 2 &&  ! readHmmer2Header(buffer))
				break;
			else if(hmmerVersion == 3 && ! readHmmer3Header(buffer))
				break;
		}

		// partial check for mandatory fields
		if (M < 1) {
			throw new InvalidFormatException();
		}
		if (name == null || name.trim().equals("")) {
			throw new InvalidFormatException();
		}

		initializeArrays();

		// skip an annotation line
		buffer = reader.readLine();

		if(hmmerVersion==2)
			readHmmer2Model(reader);
		else if(hmmerVersion==3)
			readHmmer3Model(reader);

		while ((buffer = reader.readLine()) != null) {
			if (buffer.startsWith("//"))
				break;
		}

		getRepresentativeString();

	}
	/**
	 *  Reads the fields common to all HMMER formats. It will return true until
	 * it reads the HMM field and then it will return false.
	 * @param buffer  the current line to parse
	 * @return false if the HMM field is read, true otherwise
	 * @throws IOException
	 * @throws InvalidFormatException
	 */
	public boolean readCommonHeader(String buffer) throws IOException, InvalidFormatException
	{
		String[] tokens;

		if (buffer.startsWith("NAME ")) {
			name = buffer.substring(6);
		} else if (buffer.startsWith("ACC  ")) {
			flags |= PLAN7_ACC;
			acc = buffer.substring(6);
		} else if (buffer.startsWith("DESC ")) {
			flags |= PLAN7_DESC;
			desc = buffer.substring(6);
		} else if (buffer.startsWith("NSEQ ")) {
			try {
				nseq = Integer.parseInt(buffer.substring(6));
			} catch (NumberFormatException e) {
				e.printStackTrace();
				throw new InvalidFormatException();
			}
		} else if (buffer.startsWith("LENG ")) {
			try {
				M = Integer.parseInt(buffer.substring(6));
			} catch (NumberFormatException e) {
				e.printStackTrace();
				throw new InvalidFormatException();
			}
		} else if (buffer.startsWith("ALPH ")) {
			Atype atype;
			if (buffer.substring(6).toUpperCase().equals("AMINO")) {
				atype = Atype.hmmAMINO;
				alphabetSize = 20;
				alphabetIupac = 24;
			} else if (buffer.substring(6).toUpperCase().equals("NUCLEIC")) {
				atype = Atype.hmmNUCLEIC;
				alphabetSize = 4;
				alphabetIupac = 17;
			} else {
				throw new InvalidFormatException();
			}
			setAlphabet(atype);
			alphabets = new String[alphabetSize];
		} else if (buffer.startsWith("RF   ")) { /*
														* Reference annotation
														* present?
														*/
			if (buffer.substring(6, 7).toUpperCase().equals("Y"))
				flags |= PLAN7_RF;
		} else if (buffer.startsWith("CS   ")) { /*
														* Consensus annotation
														* present?
														*/
			if (buffer.substring(6, 7).toUpperCase().equals("Y"))
				flags |= PLAN7_CS;
		} else if (buffer.startsWith("MAP  ")) { /*
														* Map annotation
														* present?
														*/
			if (buffer.substring(6, 7).toUpperCase().equals("Y"))
				flags |= PLAN7_MAP;
		} else if (buffer.startsWith("COM  ")) { /* Command line log */
			buffer = buffer.substring(6);
			comlog.add(buffer);
		} else if (buffer.startsWith("DATE ")) { /* Date file created */
			ctime = buffer.substring(6);
		} else if (buffer.startsWith("GA   ")) {
			buffer = buffer.substring(6);
			tokens = buffer.split("[\t|\n| ]+");
			if (tokens.length < 2)
				throw new InvalidFormatException();
			ga1 = Double.parseDouble(tokens[0]);
			ga2 = Double.parseDouble(tokens[1]);
			flags |= PLAN7_GA;
		} else if (buffer.startsWith("TC   ")) {
			buffer = buffer.substring(6);
			tokens = buffer.split("[\t|\n| ]+");
			if (tokens.length < 2)
				throw new InvalidFormatException();
			tc1 = Double.parseDouble(tokens[0]);
			tc2 = Double.parseDouble(tokens[1]);
			flags |= PLAN7_TC;
		} else if (buffer.startsWith("NC   ")) {
			buffer = buffer.substring(6);
			tokens = buffer.split("[\t|\n| ]+");
			if (tokens.length < 2)
				throw new InvalidFormatException();
			nc1 = Double.parseDouble(tokens[0]);
			nc2 = Double.parseDouble(tokens[1]);
			flags |= PLAN7_NC;
		}else if (buffer.startsWith("CKSUM")) {
			checksum = Long.parseLong(buffer.substring(6));
		} else if(buffer.startsWith("HMM ")) {
			tokens = buffer.split("\\s+");
			if(tokens.length < alphabetSize + 1)
				throw new InvalidFormatException();
			for(int x = 0; x < alphabetSize; x++)
				alphabets[x] = tokens[x + 1];
			return false;
		}
		return true;
	}
	/**
	 *  read header fields specific to HMMER2. Will return false if it
	 * finds the end of the header (usually the HMM field, though this 
	 * currently only apears in the common section)
	 * @param buffer
	 * @return false if end of header is found, true otherwise
	 * @throws IOException
	 * @throws InvalidFormatException
	 */
	public boolean readHmmer2Header(String buffer) throws IOException, InvalidFormatException
	{
		String[] tokens;
		if(buffer.startsWith("XT"))
		{
			tokens = buffer.split("\\s+");
			if(tokens.length < 9)
				throw new InvalidFormatException();
			int i = 1;
			for(int k = 0; k < 4; k++)
				for(int x = 0; x < 2; x++)
				{
					xt[k][x] = ascii2prob(tokens[i], 1.0);
					i++;
				}
		} else if(buffer.startsWith("NULT"))
		{
			tokens = buffer.split("\\s+");
			if(tokens.length < 3)
				throw new InvalidFormatException();
			p1 = ascii2prob(tokens[1], 1);
			p1 = p1 / (p1 + ascii2prob(tokens[2], 1.0));
		} else if(buffer.startsWith("NULE"))
		{
			tokens = buffer.split("\\s+");
			for(int x = 0; x < alphabetSize; x++)
				nullModel[x] = ascii2prob(tokens[x + 1],
						1. / (float) alphabetSize);
		} else if(buffer.startsWith("EVD  "))
		{ /* EVD parameters */
			flags |= PLAN7_EVD;
			tokens = buffer.split("[\t|\n| ]+");
			if(tokens.length < 3)
				throw new InvalidFormatException();
			mu = Double.parseDouble(tokens[1]);
			lambda = Double.parseDouble(tokens[2]);
		}
		return true;
	}
	/**
	 *  read header fields specific to HMMER2. Will return false if it
	 * finds the end of the header (usually the HMM field, though this
	 * currently only apears in the common section)
	 * @param reader false if end of header is found, true otherwise
	 * @throws IOException
	 * @throws InvalidFormatException
	 */
	public boolean readHmmer3Header(String buffer) throws IOException, InvalidFormatException
	{
		String[] tokens;

		if(buffer.startsWith("EFFN")){
			try {
				effn = Double.parseDouble(buffer.substring(6));
				flags |= PLAN7_EFFN;
			} catch (NumberFormatException e) {
				e.printStackTrace();
				throw new InvalidFormatException();
			}
		} else if(buffer.startsWith("STATS")){
			tokens = buffer.split("\\s+");
			if(tokens.length < 5 || ! tokens[1].trim().equals("LOCAL"))
				throw new InvalidFormatException();

			//TODO: check somewhere that all three of these are present, or none are
			int section;
			if(tokens[2].trim().equals("MSV"))
				section = 0;
				else if(tokens[2].trim().equals("VITERBI"))
				section = 1;
				else if(tokens[2].trim().equals("FORWARD"))
				section = 2;
				else
					throw new InvalidFormatException();

			try{
				stats[section][0] = asciiNegLn2Prob(tokens[3]);
				stats[section][1] = asciiNegLn2Prob(tokens[4]);
			}catch(NumberFormatException e){
				e.printStackTrace();
				throw new InvalidFormatException();
			}
			flags |= PLAN7_STATS;
		}
		return true;
	}
	/**
	 * parse the model section of HMMER3
	 * @param reader
	 * @throws IOException
	 * @throws InvalidFormatException
	 */
	void readHmmer3Model(BufferedReader reader) throws IOException, InvalidFormatException
	{
		String buffer;
		String[] tokens;

		buffer=reader.readLine().trim();
		if(buffer.startsWith("COMPO"))  //optional line
		{
			flags |= PLAN7_COMPO;
			compo = new double[alphabetSize];
			tokens = buffer.split("\\s+");
			if(tokens.length < alphabetSize + 1)
				throw new InvalidFormatException();
			for (int x = 0; x < alphabetSize; x++) {
				compo[x] = asciiNegLn2Prob(tokens[x+1]);
			}

			buffer=reader.readLine().trim();
		}
		// parse i_0 emission line
		tokens = buffer.split("\\s+");
		if(tokens.length < alphabetSize)
			throw new InvalidFormatException("tokens: "+tokens.length+" buffer: "+buffer);
		for (int x = 0; x < alphabetSize; x++) {
			ins[0][x] = asciiNegLn2Prob(tokens[x]);
		}

		// parse 0 state transitions
		buffer=reader.readLine().trim();
		tokens = buffer.split("\\s+");
		if(tokens.length < 7)
			throw new InvalidFormatException();
		for (int x = 0; x < 7; x++) 
			t[0][x] = asciiNegLn2Prob(tokens[x]);

		for (int k = 1; k <= M; k++) {

			// Line 1:  match emissions, map, rf, and cs
			buffer=reader.readLine().trim();
			tokens = buffer.split("\\s+");
			if(tokens.length < alphabetSize+3)
				throw new InvalidFormatException();
			for (int x = 0; x < alphabetSize; x++)
				mat[k][x] = asciiNegLn2Prob(tokens[x+1]);
			if( ( flags & PLAN7_MAP) != 0)
				map[k] = Integer.parseInt(tokens[alphabetSize + 1]);
			if( ( flags & PLAN7_RF) != 0)
				rf[k] = tokens[alphabetSize + 2];
			if( ( flags & PLAN7_CS) != 0)
				cs[k] = tokens[alphabetSize + 3];

			// Line 2:  insert emissions
			buffer=reader.readLine().trim();
			tokens = buffer.split("\\s+");
			if(tokens.length < alphabetSize)
				throw new InvalidFormatException();
			// hmmer3 does not use given emissons for insert states
			//loadDefaultBackground(ins[k]);
			for (int x = 0; x < alphabetSize; x++)
				ins[k][x] = asciiNegLn2Prob(tokens[x]);

			// Line 3: 7 transitions
			buffer=reader.readLine().trim();
			tokens = buffer.split("\\s+");
			if(tokens.length < 7)
				throw new InvalidFormatException();
			for (int x = 0; x < 7; x++)
				t[k][x] = asciiNegLn2Prob(tokens[x]);
		}
		// fill in other stuff
		xt[XTE][MOVE]=1.0;
		xt[XTE][LOOP]=0.0;
		loadDefaultBackground(nullModel);

		for(int k=1; k <= M; k++)
			begin[k]=end[k]=0.0;

		// hmmer3 includes a i_0 state which we remove and renomalize
	    double norm = t[0][TMM] + t[0][TMD];
		begin[1] = t[0][TMM] / norm; //M_0 -> M_1
		tbd1 = t[0][TMD] / norm;  // M_0 ->  D_1

		//hmmer3 includes a i_M state which we remove and renomalize
		//t[M][TMM]=t[M][TMI]=0.0;
		end[M]=1.0;

		//null transition and N,C,J transitions are set per sequence

	}
	/**
	 * parse the model section of HMMER2
	 * @param reader
	 * @throws IOException
	 * @throws InvalidFormatException
	 */
	void readHmmer2Model(BufferedReader reader) throws IOException, InvalidFormatException
	{
		String buffer;
		String[] tokens;

		/*
		 * Main model section. Read as integer log odds, convert to
		 * probabilities
		 */

		// parse tbd1 line
		buffer = reader.readLine();
		buffer = buffer.trim();
		tokens = buffer.split("\\s+");
		if (tokens.length < 3) {
			throw new InvalidFormatException();
		}
		double p = ascii2prob(tokens[0], 1.0);
		tbd1 = ascii2prob(tokens[2], 1.0);
		tbd1 = tbd1 / (p + tbd1);

		// main model
		for (int k = 1; k <= M; k++) {
			// Line 1: k, match emissions, map
			buffer = reader.readLine();
			buffer = buffer.trim();
			tokens = buffer.split("\\s+");
			if (tokens.length < alphabetSize + 2) {
				throw new InvalidFormatException();
			}
			for (int x = 0; x < alphabetSize; x++) {
				mat[k][x] = ascii2prob(tokens[x + 1], nullModel[x]);
			}
			map[k] = Integer.parseInt(tokens[alphabetSize + 1]);

			// Line 2: RF and insert emissions
			buffer = reader.readLine();
			buffer = buffer.trim();
			tokens = buffer.split("\\s+");
			if (tokens.length < alphabetSize + 1) {
				throw new InvalidFormatException();
			}
			if ((flags & PLAN7_RF) != 0)
				rf[k] = tokens[0];
			if (k < M) {
				for (int x = 0; x < alphabetSize; x++) {
					ins[k][x] = ascii2prob(tokens[x + 1], nullModel[x]);
				}
			}
			// Line 3: CS and transitions
			buffer = reader.readLine();
			buffer = buffer.trim();
			tokens = buffer.split("\\s+");
			if (tokens.length < 10) {
				throw new InvalidFormatException();
			}
			if ((flags & PLAN7_CS) != 0)
				cs[k] = tokens[0];
			for (int x = 0; x < 7; x++) {
				if (k < M)
					t[k][x] = ascii2prob(tokens[x + 1], 1.0);
			}
			begin[k] = ascii2prob(tokens[8], 1.0);
			end[k] = ascii2prob(tokens[9], 1.0);
		}

	}


	/**
	 * write the model to the given file
	 * @param hmmFilename
	 * @throws IOException
	 */
	public void write(String hmmFilename) throws IOException {
		FileWriter fileWriter = new FileWriter(hmmFilename);
		write(fileWriter);
		fileWriter.close();
	}
	/**
	 * write the model to the given stream. The caller is responsble 
	 * for closing the stream.
	 * @param w
	 * @throws IOException
	 */
	public void write(Writer w) throws IOException 
	{
		PrintWriter writer = new PrintWriter(w);

		if(hmmerVersion==2)
			writer.println("HMMER2.0  [" + PACKAGE_VERSION + "]"); /* magic header */
		else if(hmmerVersion==3)
			writer.println("HMMER3/b "); /* magic header */

		writeCommonHeader(writer);

		if(hmmerVersion==2)
			writeHmmer2Header(writer);
		else if(hmmerVersion==3)
			writeHmmer3Header(writer);


		if(hmmerVersion==2)
			writeHmmer2Model(writer);
		else if(hmmerVersion==3)
			writeHmmer3Model(writer);

		writer.println("//");
	}
	/**
	 *  write header feilds common to all formats
	 * @param writer
	 * @throws IOException
	 */
	void writeCommonHeader(PrintWriter writer) throws IOException
	{
		writer.println("NAME  " + name);
		if (acc != null)
			writer.println("ACC   " + acc);
		if (desc != null)
			writer.println("DESC  " + desc);
		writer.println("LENG  " + M);
		if (alphabetType == Atype.hmmAMINO) {
			writer.println("ALPH  " + "Amino");
		} else {
			writer.println("ALPH  " + "Nucleic");
		}
		if ((flags & PLAN7_RF) != 0)
			writer.println("RF    " + "yes");
		else
			writer.println("RF    " + "no");
		if ((flags & PLAN7_CS) != 0)
			writer.println("CS    " + "yes");
		else
			writer.println("CS    " + "no");
		if ((flags & PLAN7_MAP) != 0)
			writer.println("MAP   " + "yes");
		else
			writer.println("MAP   " + "no");
		for (int x = 0; x < comlog.size(); x++) {
			writer.println("COM   " + comlog.get(x));
		}
		writer.println("NSEQ  " + nseq);
		writer.println("DATE  " + ctime);
		writer.println("CKSUM " + checksum);
		if ((flags & PLAN7_GA) != 0)
			writer.format("GA    %.2f %.2f%n", ga1, ga2);
		if ((flags & PLAN7_TC) != 0)
			writer.format("TC    %.2f %.2f%n", tc1, tc2);
		if ((flags & PLAN7_NC) != 0)
			writer.format("NC    %.2f %.2f%n", nc1, nc2);


	}
	/**
	 *  write hmmer2 specific headers
	 * @param writer
	 * @throws IOException
	 */
	void writeHmmer2Header(PrintWriter writer) throws IOException
	{
		int k; /* counter for nodes */
		int x; /* counter for symbols */
		/*
		 * Specials
		 */
		writer.print("XT     ");
		for (k = 0; k < 4; k++)
			for (x = 0; x < 2; x++)
				writer.format("%6s ", prob2ascii(xt[k][x], 1.0));
		writer.println();

		/*
		 * Save the null model first, so HMM readers can decode log odds scores
		 * on the fly. Save as log odds probabilities relative to
		 * 1/Alphabet_size (flat distribution)
		 */
		writer.print("NULT  ");
		writer.format("%6s ", prob2ascii(p1, 1.0)); /* p1 */
		writer.format("%6s%n", prob2ascii(1.0 - p1, 1.0)); /* p2 */
		writer.print("NULE  ");
		for ( x = 0; x < alphabetSize; x++)
			writer.format("%6s ", prob2ascii(nullModel[x],
					1 / (float) (alphabetSize)));
		writer.println();

		/*
		 * EVD statistics
		 */
		if ((flags & PLAN7_EVD) != 0 && !modified) // invalidate EVD when
														// hmm is modified
			writer.format("EVD   %10f %10f%n", mu, lambda);


	}
	/**
	 *  write hmmer3 specific headers
	 * @param writer
	 * @throws IOException
	 */
	void writeHmmer3Header(PrintWriter writer) throws IOException
	{
		if( (flags & PLAN7_EFFN) != 0)
			writer.format("EFFN  %5f%n",effn);
		if( (flags & PLAN7_STATS) != 0)
		{
			writer.format("%5s %5s %-7s %9s %8s%n","STATS","LOCAL","MSV", prob2NegLnAscii(stats[0][0]), prob2NegLnAscii(stats[0][1]));
			writer.format("%5s %5s %-7s %9s %8s%n",	"STATS","LOCAL","VITERBI",   prob2NegLnAscii(stats[1][0]), prob2NegLnAscii(stats[1][1]));
			writer.format("%5s %5s %-7s %9s %8s%n","STATS","LOCAL","FORWARD",prob2NegLnAscii(stats[2][0]), prob2NegLnAscii(stats[2][1]));
		}
	}
	/**
	 *  write model for hmmer2 format
	 * @param writer
	 * @throws IOException
	 */
	void writeHmmer2Model(PrintWriter writer) throws IOException
	{
		int k; /* counter for nodes */
		int x; /* counter for symbols */
		int ts; /* counter for state transitions */

		writer.print("HMM  ");
		for (x = 0; x < alphabetSize; x++)
			writer.format("%7s", alphabets[x]);
		writer.println();
		writer.format("       %6s %6s %6s %6s %6s %6s %6s %6s %6s%n", "m->m",
				"m->i", "m->d", "i->m", "i->i", "d->m", "d->d", "b->m", "m->e");

		/*
		 * Print HMM parameters (main section of the save file)
		 */
		writer.format("       %6s %6s ", prob2ascii(1 - tbd1, 1.0), "*");
		writer.format("%6s%n", prob2ascii(tbd1, 1.0));
		for (k = 1; k <= M; k++) {
			/* Line 1: k, match emissions, map */
			writer.format(" %5d ", k);
			for (x = 0; x < alphabetSize; x++)
				writer.format("%6s ", prob2ascii(mat[k][x], nullModel[x]));
			if ((flags & PLAN7_MAP) != 0)
				writer.format("%5d", map[k]);
			writer.println();
			/* Line 2: RF and insert emissions */
			if ((flags & PLAN7_RF) != 0)
				writer.format(" %5s ", rf[k]);
			else
				writer.format(" %5s ", "-");
			for (x = 0; x < alphabetSize; x++) {
				if (k < M)
					writer.format("%6s ", prob2ascii(ins[k][x], nullModel[x]));
				else
					writer.format("%6s ", "*");
			}

			writer.println();
			/* Line 3: CS and transition probs */
			if ((flags & PLAN7_CS) != 0)
				writer.format(" %5s ", cs[k]);
			else
				writer.format(" %5s ", "-");
			for (ts = 0; ts < 7; ts++) {
				if (k < M)
					writer.format("%6s ", prob2ascii(t[k][ts], 1.0));
				else
					writer.format("%6s ", "*");
			}

			writer.format("%6s ", prob2ascii(begin[k], 1.0));
			writer.format("%6s ", prob2ascii(end[k], 1.0));

			writer.println();
		}

	}
	/**
	 * write model for hmmer3 format
	 * @param writer
	 * @throws IOException
	 */
	void writeHmmer3Model(PrintWriter writer) throws IOException
	{
		int k; /* counter for nodes */
		int x; /* counter for symbols */
		int ts; /* counter for state transitions */

		writer.print("HMM  ");
		for (x = 0; x < alphabetSize; x++)
			writer.format("%9s", alphabets[x]);
		writer.println();
		writer.format("       %9s%9s%9s%9s%9s%9s%9s%n",
				"m->m", "m->i", "m->d", "i->m", "i->i", "d->m", "d->d");

		if( (flags & PLAN7_COMPO) != 0){
			writer.format("%7s ","COMPO");
			for(x=0; x < alphabetSize; x++)
				writer.format("%9s",prob2NegLnAscii(compo[x]));
			writer.println();
		}

		// write i_0 emissions
		writer.print("        ");
		for(x=0; x < alphabetSize; x++)
			writer.format("%9s",prob2NegLnAscii(ins[0][x]));
		writer.println();

		// write state 0 transitions
		writer.print("        ");
		for(x=0; x < 7; x++)
			writer.format("%9s",prob2NegLnAscii(t[0][x]));
		writer.println();

		for(k=1; k <= M; k++) {
			// Line 1 match emissions
			writer.format("%7d ",k);
			for(x=0; x < alphabetSize; x++)
				writer.format("%9s",prob2NegLnAscii(mat[k][x]));
			writer.format("%7d %s %s%n",
					(flags & PLAN7_MAP ) !=0 ? map[k]:"-",
					(flags & PLAN7_RF) !=0 ? rf[k]:"-",
					(flags & PLAN7_CS) !=0 ? cs[k]:"-");

			//Line 2 insert emissions
			writer.format("%7s ","");
			for(x=0; x < alphabetSize; x++)
				writer.format("%9s",prob2NegLnAscii(ins[k][x]));
			writer.println();

			// Line 3 transitions
			writer.format("%7s ","");
			for(ts=0; ts < 7; ts++)
				writer.format("%9s",prob2NegLnAscii(t[k][ts]));
			writer.println();
		}
	}

	static double FSum(double vec[], int n) {
		return FSum(vec, 0, n);
	}

	static void FNorm(double vec[], int n) {
		FNorm(vec, 0, n);
	}

	static double FSum(double vec[], int s, int n) {
		double sum = 0;
		for (int x = s; x < s + n; x++)
			sum += vec[x];
		return sum;
	}

	static void FNorm(double vec[], int s, int n) {
		int x;
		double sum;

		sum = FSum(vec, n);
		if (sum != 0.0)
			for (x = s; x < s + n; x++)
				vec[x] /= sum;
		else
			for (x = s; x < s + n; x++)
				vec[x] = 1. / (double) n;
	}

	static void FScale(double vec[], int n, double scale) {
		FScale(vec, 0, n, scale);
	}

	static void FScale(double vec[], int s, int n, double scale) {
		int x;
		for (x = s; x < s + n; x++)
			vec[x] *= scale;
	}

	/**
	 * Take an HMM in counts form, and renormalize all of its
	 * probability vectors. Also enforces Plan7 restrictions on nonexistent
	 * transitions.
	 * 
	 */
	void Plan7Renormalize() {
		int k; /* counter for model position */
		int st; /* counter for special states */
		double d; /* denominator */

		/* match emissions */
		for (k = 1; k <= M; k++)
			FNorm(mat[k], alphabetSize);
		/* insert emissions */
		for (k = 1; k < M; k++)
			FNorm(ins[k], alphabetSize);
		/* begin transitions */
		d = FSum(begin, 1, M) + tbd1;
		FScale(begin, 1, M, 1. / d);
		tbd1 /= d;
		/* main model transitions */
		for (k = 1; k < M; k++) {
			d = FSum(t[k], 3) + end[k];
			FScale(t[k], 3, 1. / d);
			end[k] /= d;

			FNorm(t[k], 3, 2); /* insert */
			FNorm(t[k], 5, 2); /* delete */
		}
		/* null model emissions */
		FNorm(nullModel, alphabetSize);
		/* special transitions */
		for (st = 0; st < 4; st++)
			FNorm(xt[st], 2);
		/* enforce nonexistent transitions */
		/* (is this necessary?) */
		t[0][TDM] = t[0][TDD] = 0.0;

		flags &= ~PLAN7_HASBITS; /* clear the log-odds ready flag */
		flags |= PLAN7_HASPROB; /* set the probabilities OK flag */
	}

	/**
	 * Generate consensus string of the hmm
	 *
	 */
	void getRepresentativeString() {
		consensusP = new double[M + 1];
		String result = "";
		for (int i = 1; i <= M; i++) {
			double p = 0;
			int max = -1;
			for (int j = 0; j < alphabetSize; j++) {
				// use >= here so if mat[i] is all zeros we still get a valid max index
				if (mat[i][j] >= p) {
					p = mat[i][j];
					max = j;
				}
			}
			if (alphabetSize == 20) {
				if (p > 0.5) {
					result += alphabets[max].toUpperCase();
				} else {
					result += alphabets[max].toLowerCase();
				}
				consensusP[i] = p;
			}
		}
		consensus = result;
	}

	/**
	 * Given a sequence character x and an hmm emission probability
	 * vector, calculate the log-odds (base 2) score of the symbol.
	 * 
	 * Easy if x is in the emission alphabet, but not so easy is x is a
	 * degenerate symbol. The "correct" Bayesian philosophy is to calculate
	 * score(X) by summing over p(x) for all x in the degenerate symbol X to get
	 * P(X), doing the same sum over the prior to get F(X), and doing log_2
	 * (P(X)/F(X)). This gives an X a zero score, for instance.
	 * 
	 * Though this is correct in a formal Bayesian sense -- we have no
	 * information on the sequence, so we can't say if it's random or model, so
	 * it scores zero -- it sucks, big time, for scoring biological sequences.
	 * Sequences with lots of X's score near zero, while real sequences have
	 * average scores that are negative -- so the X-laden sequences appear to be
	 * lifted out of the noise of a full histogram of a database search. Correct
	 * or not, this is highly undesirable.
	 * 
	 * So therefore we calculated the expected score of the degenerate symbol by
	 * summing over all x in X: e_x log_2 (p(x)/f(x)) where the expectation of
	 * x, e_x, is calculated from the random model.
	 * 
	 * Empirically, this works; it also has a wooly hand-waving probabilistic
	 * justification that I'm happy enough about.
	 * 
	 * Args: p - probabilities of normal symbols null - null emission model
	 * ambig - index of the degenerate character in Alphabet[]
	 * 
	 * Return: the integer log odds score of x given the emission vector and the
	 * null model, scaled up by INTSCALE.
	 */
	int DegenerateSymbolScore(double[] p, double nullModel[], int ambig) {
		int x;
		double numer = 0.0;
		double denom = 0.0;

		for (x = 0; x < alphabetSize; x++) {
			if (Degenerate[ambig][x] != 0) {
				numer += nullModel[x] * sreLOG2(p[x] / nullModel[x]);
				denom += nullModel[x];
			}
		}
		return (int) (INTSCALE * numer / denom);
	}

	/**
	 * Take an HMM with valid probabilities, and fill in the integer
	 * log-odds score section of the model.
	 * 
	 * Notes on log-odds scores: type of parameter probability score
	 * ----------------- ----------- ------ any emission p_x log_2 p_x/null_x
	 * N,J,C /assume/ p_x = null_x so /always/ score zero. transition to
	 * emitters t_x log_2 t_x/p1 (M,I; N,C; J) NN and CC loops are often equal
	 * to p1, so usu. score zero. C->T transition t_x log_2 t_x/p2 often zero,
	 * usu. C->T = p2. all other transitions t_x log_2 t_x (no null model
	 * counterpart, so null prob is 1)
	 * 
	 * Notes on entry/exit scores, B->M and M->E: The probability form model
	 * includes delete states 1 and M. these states are removed from a search
	 * form model to prevent B->D...D->E->J->B mute cycles, which would
	 * complicate dynamic programming algorithms. The data-independent S/W B->M
	 * and M->E transitions are folded together with data-dependent B->D...D->M
	 * and M->D...D->E paths.
	 * 
	 * This process is referred to in the code as "wing folding" or "wing
	 * retraction"... the analogy is to a swept-wing fighter in landing vs. high
	 * speed flight configuration.
	 * 
	 * Note on Viterbi vs. forward flag: Wing retraction must take forward vs.
	 * Viterbi into account. If forward, sum two paths; if Viterbi, take max. I
	 * tried to slide this by as a sum, without the flag, but Alex detected it
	 * as a bug, because you can then find cases where the Viterbi score doesn't
	 * match the P7TraceScore().
	 * 
	 */
	void P7Logoddsify(int viterbi_mode) {
		int k; /* counter for model position */
		int x; /* counter for symbols */
		double accum;
		double tbm, tme;

		if ((flags & PLAN7_HASBITS) != 0)
			return;

		/*
		 * Symbol emission scores
		 */
		for (k = 1; k <= M; k++) {
			/* match/insert emissions in main model */
			for (x = 0; x < alphabetSize; x++) {
				msc[x][k] = prob2Score(mat[k][x], nullModel[x]);
				if (k < M)
					isc[x][k] = prob2Score(ins[k][x], nullModel[x]);
			}
			/* degenerate match/insert emissions */
			for (x = alphabetSize; x < alphabetIupac; x++) {
				msc[x][k] = DegenerateSymbolScore(mat[k], nullModel, x);
				if (k < M)
					isc[x][k] = DegenerateSymbolScore(ins[k], nullModel, x);
			}
		}

		/*
		 * State transitions.
		 * 
		 * A note on "folding" of D_1 and D_M. These two delete states are
		 * folded out of search form models in order to prevent null cycles in
		 * the dynamic programming algorithms (see code below). However, we use
		 * their log transitions when we save the model! So the following log
		 * transition probs are used *only* in save files, *never* in search
		 * algorithms: log (tbd1), D1 -> M2, D1 -> D2 Mm-1 -> Dm, Dm-1 -> Dm
		 * 
		 * In a search algorithm, these have to be interpreted as -INFTY because
		 * their contributions are folded into bsc[] and esc[] entry/exit
		 * scores. They can't be set to -INFTY here because we need them in save
		 * files.
		 */
		for (k = 1; k < M; k++) {
			tsc[TMM][k] = prob2Score(t[k][TMM], p1);
			tsc[TMI][k] = prob2Score(t[k][TMI], p1);
			tsc[TMD][k] = prob2Score(t[k][TMD], 1.0);
			tsc[TIM][k] = prob2Score(t[k][TIM], p1);
			tsc[TII][k] = prob2Score(t[k][TII], p1);
			tsc[TDM][k] = prob2Score(t[k][TDM], p1);
			tsc[TDD][k] = prob2Score(t[k][TDD], 1.0);
		}

		/*
		 * B->M entry transitions. Note how D_1 is folded out. M1 is just B->M1
		 * M2 is sum (or max) of B->M2 and B->D1->M2 M_k is sum (or max) of
		 * B->M_k and B->D1...D_k-1->M_k These have to be done in log space,
		 * else you'll get underflow errors; and we also have to watch for
		 * log(0). A little sloppier than it probably has to be; historically,
		 * doing in this in log space was in response to a bug report.
		 */
		accum = tbd1 > 0.0 ? Math.log(tbd1) : -9999.;
		for (k = 1; k <= M; k++) {
			tbm = begin[k] > 0. ? Math.log(begin[k]) : -9999.; /* B->M_k part */

			/* B->D1...D_k-1->M_k part we get from accum */
			if (k > 1 && accum > -9999.) {
				if (t[k - 1][TDM] > 0.0) {
					if (viterbi_mode != 0)
						tbm = Math.max(tbm, accum + Math.log(t[k - 1][TDM]));
					else
						tbm = LogSum(tbm, accum + Math.log(t[k - 1][TDM]));
				}

				accum = (t[k - 1][TDD] > 0.0) ? accum + Math.log(t[k - 1][TDD])
						: -9999.;
			}
			/* Convert from log_e to scaled integer log_2 odds. */
			if (tbm > -9999.)
				bsc[k] = (int) Math.floor(0.5 + INTSCALE * 1.44269504
						* (tbm - Math.log(p1)));
			else
				bsc[k] = -INFTY;
		}

		/*
		 * M->E exit transitions. Note how D_M is folded out. M_M is 1 by
		 * definition M_M-1 is sum of M_M-1->E and M_M-1->D_M->E, where D_M->E
		 * is 1 by definition M_k is sum of M_k->E and M_k->D_k+1...D_M->E Must
		 * be done in log space to avoid underflow errors. A little sloppier
		 * than it probably has to be; historically, doing in this in log space
		 * was in response to a bug report.
		 */
		esc[M] = 0;
		accum = 0.;
		for (k = M - 1; k >= 1; k--) {
			tme = end[k] > 0. ? Math.log(end[k]) : -9999.;
			if (accum > -9999.) {
				if (t[k][TMD] > 0.0) {
					if (viterbi_mode != 0)
						tme = Math.max(tme, accum + Math.log(t[k][TMD]));
					else
						tme = LogSum(tme, accum + Math.log(t[k][TMD]));
				}
				accum = (t[k][TDD] > 0.0) ? accum + Math.log(t[k][TDD])
						: -9999.;
			}
			/* convert from log_e to scaled integer log odds. */
			esc[k] = (tme > -9999.) ? (int) Math.floor(0.5 + INTSCALE
					* 1.44269504 * tme) : -INFTY;
		}

		/* special transitions */
		xsc[XTN][LOOP] = prob2Score(xt[XTN][LOOP], p1);
		xsc[XTN][MOVE] = prob2Score(xt[XTN][MOVE], 1.0);
		xsc[XTE][LOOP] = prob2Score(xt[XTE][LOOP], 1.0);
		xsc[XTE][MOVE] = prob2Score(xt[XTE][MOVE], 1.0);
		xsc[XTC][LOOP] = prob2Score(xt[XTC][LOOP], p1);
		xsc[XTC][MOVE] = prob2Score(xt[XTC][MOVE], 1. - p1);
		xsc[XTJ][LOOP] = prob2Score(xt[XTJ][LOOP], p1);
		xsc[XTJ][MOVE] = prob2Score(xt[XTJ][MOVE], 1.0);

		flags |= PLAN7_HASBITS; /* raise the log-odds ready flag */
	}

	/**
	 * Internal representation of a sequence in HMMER is as a char
	 * array. 1..L are the indices of seq symbols in Alphabet[]. 0,L+1 are
	 * sentinel bytes, set to be Alphabet_iupac -- i.e. one more than the
	 * maximum allowed index.
	 * 
	 * Assumes that 'X', the fully degenerate character, is the last character
	 * in the allowed alphabet.
	 * 
	 * Args: seq - sequence to be digitized (0..L-1) L - length of sequence
	 * 
	 * Return: digitized sequence, dsq. dsq is allocated here and must be free'd
	 * by caller.
	 */
	public int[] DigitizeSequence(String seq) {
		seq = seq.toUpperCase();
		int[] dsq = new int[seq.length() + 2];

		dsq[0] = dsq[seq.length() + 1] = alphabetIupac;
		for (int i = 1; i <= seq.length(); i++)
			dsq[i] = alphabet.indexOf(seq.charAt(i - 1));
		return dsq;
	}

	/**
	 * Convert a scaled integer log-odds score to a floating point
	 * score for output. (could be a macro but who cares.)
	 */
	double Scorify(int sc) {
		return ((double) sc / INTSCALE);
	}

	/**
	 * The Viterbi dynamic programming algorithm. Identical to
	 * Forward() except that max's replace sum's.
	 * 
	 * This is the slower, more understandable version of P7Viterbi(). The
	 * default version in fast_algorithms.c is portably optimized and more
	 * difficult to understand; the ALTIVEC version in fast_algorithms.c is
	 * vectorized with Altivec-specific code, and is pretty opaque.
	 * 
	 * This function is not enabled by default; it is only activated by -DSLOW
	 * at compile time.
	 * 
	 * @param dsq sequence in digitized form 
	 * @param L length of dsq hmm - the model
	 * @param mx reused DP matrix ret_tr - RETURN: traceback; pass NULL if it's not
	 * wanted
	 * 
	 * Return: trace of the sequence
	 */
	P7Trace P7Viterbi(int[] dsq, int L, DPMatrix mx) {
		int i, k;
		int sc;

		if(hmmerVersion==3) //hmmer3
			throw new IllegalStateException("The viterbi algorthm does not work with hmmer3 models");

		/*
		 * Allocate a DP matrix with 0..L rows, 0..M-1 columns.
		 */
		mx.ResizePlan7Matrix(L, M);
		int xmx[][] = mx.xmx;
		int mmx[][] = mx.mmx;
		int imx[][] = mx.imx;
		int dmx[][] = mx.dmx;

		/*
		 * Initialization of the zero row.
		 */
		xmx[0][XMN] = 0; /* S->N, p=1 */
		xmx[0][XMB] = xsc[XTN][MOVE]; /* S->N->B, no N-tail */
		xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY; /*
															 * need seq to get
															 * here
															 */
		for (k = 0; k <= M; k++)
			mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY; /* need seq to get here */

		/*
		 * Recursion. Done as a pull. Note some slightly wasteful boundary
		 * conditions: tsc[0] = -INFTY for all eight transitions (no node 0) D_M
		 * and I_M are wastefully calculated (they don't exist)
		 */
		for (i = 1; i <= L; i++) {
			mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;

			for (k = 1; k <= M; k++) {
				/* match state */
				mmx[i][k] = -INFTY;
				if ((sc = mmx[i - 1][k - 1] + tsc[TMM][k - 1]) > mmx[i][k])
					mmx[i][k] = sc;
				if ((sc = imx[i - 1][k - 1] + tsc[TIM][k - 1]) > mmx[i][k])
					mmx[i][k] = sc;
				if ((sc = xmx[i - 1][XMB] + bsc[k]) > mmx[i][k])
					mmx[i][k] = sc;
				if ((sc = dmx[i - 1][k - 1] + tsc[TDM][k - 1]) > mmx[i][k])
					mmx[i][k] = sc;
				if (msc[dsq[i]][k] != -INFTY)
					mmx[i][k] += msc[dsq[i]][k];
				else
					mmx[i][k] = -INFTY;

				/* delete state */
				dmx[i][k] = -INFTY;
				if ((sc = mmx[i][k - 1] + tsc[TMD][k - 1]) > dmx[i][k])
					dmx[i][k] = sc;
				if ((sc = dmx[i][k - 1] + tsc[TDD][k - 1]) > dmx[i][k])
					dmx[i][k] = sc;

				/* insert state */
				if (k < M) {
					imx[i][k] = -INFTY;
					if ((sc = mmx[i - 1][k] + tsc[TMI][k]) > imx[i][k])
						imx[i][k] = sc;
					if ((sc = imx[i - 1][k] + tsc[TII][k]) > imx[i][k])
						imx[i][k] = sc;
					if (isc[dsq[i]][k] != -INFTY)
						imx[i][k] += isc[dsq[i]][k];
					else
						imx[i][k] = -INFTY;
				}
			}

			/*
			 * Now the special states. Order is important here. remember, C and
			 * J emissions are zero score by definition,
			 */
			/* N state */
			xmx[i][XMN] = -INFTY;
			if ((sc = xmx[i - 1][XMN] + xsc[XTN][LOOP]) > -INFTY)
				xmx[i][XMN] = sc;

			/* E state */
			xmx[i][XME] = -INFTY;
			for (k = 1; k <= M; k++)
				if ((sc = mmx[i][k] + esc[k]) > xmx[i][XME])
					xmx[i][XME] = sc;
			/* J state */
			xmx[i][XMJ] = -INFTY;
			if ((sc = xmx[i - 1][XMJ] + xsc[XTJ][LOOP]) > -INFTY)
				xmx[i][XMJ] = sc;
			if ((sc = xmx[i][XME] + xsc[XTE][LOOP]) > xmx[i][XMJ])
				xmx[i][XMJ] = sc;

			/* B state */
			xmx[i][XMB] = -INFTY;
			if ((sc = xmx[i][XMN] + xsc[XTN][MOVE]) > -INFTY)
				xmx[i][XMB] = sc;
			if ((sc = xmx[i][XMJ] + xsc[XTJ][MOVE]) > xmx[i][XMB])
				xmx[i][XMB] = sc;

			/* C state */
			xmx[i][XMC] = -INFTY;
			if ((sc = xmx[i - 1][XMC] + xsc[XTC][LOOP]) > -INFTY)
				xmx[i][XMC] = sc;
			if ((sc = xmx[i][XME] + xsc[XTE][MOVE]) > xmx[i][XMC])
				xmx[i][XMC] = sc;
		}
		/* T state (not stored) */
		sc = xmx[L][XMC] + xsc[XTC][MOVE];

		P7Trace tr;
		tr = P7ViterbiTrace(dsq, L, mx);

		tr.score = Scorify(sc); /* the total Viterbi score. */
		return tr;
	}

	/**
	 * Traceback of a Viterbi matrix: i.e. retrieval of optimum
	 * alignment.
	 * 
	 * @param dsq sequence aligned to (digital form) 1..N 
	 * @param N length of seq
	 * @param mx - the matrix to trace back in
	 * @return traceback.
	 * 
	 */
	P7Trace P7ViterbiTrace(int[] dsq, int N, DPMatrix mx) {
		P7Trace tr;
		int curralloc; /* current allocated length of trace */
		int tpos; /* position in trace */
		int i; /* position in seq (1..N) */
		int k = 0; /* position in model (1..M) */
		int xmx[][], mmx[][], imx[][], dmx[][];
		int sc=-1; /* temp var for pre-emission score */

		/*
		 * Overallocate for the trace. S-N-B- ... - E-C-T : 6 states + N is
		 * minimum trace; add N more as buffer.
		 */
		curralloc = N * 2 + 6;
		tr = P7Trace.P7AllocTrace(curralloc);

		xmx = mx.xmx;
		mmx = mx.mmx;
		imx = mx.imx;
		dmx = mx.dmx;

		/*
		 * Initialization of trace We do it back to front; ReverseTrace() is
		 * called later.
		 */
		tr.statetype[0] = P7Trace.STT;
		tr.nodeidx[0] = 0;
		tr.pos[0] = 0;
		tr.statetype[1] = P7Trace.STC;
		tr.nodeidx[1] = 0;
		tr.pos[1] = 0;
		tpos = 2;
		i = N; /* current i (seq pos) we're trying to assign */

		/*
		 * Traceback
		 */
		while (tr.statetype[tpos - 1] != P7Trace.STS) {
			switch (tr.statetype[tpos - 1]) {
			case P7Trace.STM: /* M connects from i-1,k-1, or B */
				sc = mmx[i + 1][k + 1] - msc[dsq[i + 1]][k + 1];
				if (sc <= -INFTY) {
					return null;
				} else if (sc == xmx[i][XMB] + bsc[k + 1]) {
					/* Check for wing unfolding */
					if (prob2Score(begin[k + 1], p1) + 1 * INTSCALE <= bsc[k + 1])
						while (k > 0) {
							tr.statetype[tpos] = P7Trace.STD;
							tr.nodeidx[tpos] = k--;
							tr.pos[tpos] = 0;
							tpos++;
							if (tpos == curralloc) { /*
														 * grow trace if
														 * necessary
														 */
								curralloc += N;
								tr.P7ReallocTrace(curralloc);
							}
						}

					tr.statetype[tpos] = P7Trace.STB;
					tr.nodeidx[tpos] = 0;
					tr.pos[tpos] = 0;
				} else if (sc == mmx[i][k] + tsc[TMM][k]) {
					tr.statetype[tpos] = P7Trace.STM;
					tr.nodeidx[tpos] = k--;
					tr.pos[tpos] = i--;
				} else if (sc == imx[i][k] + tsc[TIM][k]) {
					tr.statetype[tpos] = P7Trace.STI;
					tr.nodeidx[tpos] = k;
					tr.pos[tpos] = i--;
				} else if (sc == dmx[i][k] + tsc[TDM][k]) {
					tr.statetype[tpos] = P7Trace.STD;
					tr.nodeidx[tpos] = k--;
					tr.pos[tpos] = 0;
				} else
					throw new RuntimeException("traceback failed");
				break;

			case P7Trace.STD: /* D connects from M,D */
				if (dmx[i][k + 1] <= -INFTY) {
					return null;
				} else if (dmx[i][k + 1] == mmx[i][k] + tsc[TMD][k]) {
					tr.statetype[tpos] = P7Trace.STM;
					tr.nodeidx[tpos] = k--;
					tr.pos[tpos] = i--;
				} else if (dmx[i][k + 1] == dmx[i][k] + tsc[TDD][k]) {
					tr.statetype[tpos] = P7Trace.STD;
					tr.nodeidx[tpos] = k--;
					tr.pos[tpos] = 0;
				} else
					throw new RuntimeException("traceback failed");
				break;

			case P7Trace.STI: /* I connects from M,I */
				sc = imx[i + 1][k] - isc[dsq[i + 1]][k];
				if (sc <= -INFTY) {
					return null;
				} else if (sc == mmx[i][k] + tsc[TMI][k]) {
					tr.statetype[tpos] = P7Trace.STM;
					tr.nodeidx[tpos] = k--;
					tr.pos[tpos] = i--;
				} else if (sc == imx[i][k] + tsc[TII][k]) {
					tr.statetype[tpos] = P7Trace.STI;
					tr.nodeidx[tpos] = k;
					tr.pos[tpos] = i--;
				} else
					throw new RuntimeException("traceback failed");
				break;

			case P7Trace.STN: /* N connects from S, N */
				if (i == 0 && xmx[i][XMN] == 0) {
					tr.statetype[tpos] = P7Trace.STS;
					tr.nodeidx[tpos] = 0;
					tr.pos[tpos] = 0;
				} else if (i > 0
						&& xmx[i + 1][XMN] == xmx[i][XMN] + xsc[XTN][LOOP]) {
					tr.statetype[tpos] = P7Trace.STN;
					tr.nodeidx[tpos] = 0;
					tr.pos[tpos] = 0; /* note convention adherence: */
					tr.pos[tpos - 1] = i--; /* first N doesn't emit */
				} else
					throw new RuntimeException("traceback failed");
				break;

			case P7Trace.STB: /* B connects from N, J */
				if (xmx[i][XMB] <= -INFTY) {
					return null;
				} else if (xmx[i][XMB] == xmx[i][XMN] + xsc[XTN][MOVE]) {
					tr.statetype[tpos] = P7Trace.STN;
					tr.nodeidx[tpos] = 0;
					tr.pos[tpos] = 0;
				} else if (xmx[i][XMB] == xmx[i][XMJ] + xsc[XTJ][MOVE]) {
					tr.statetype[tpos] = P7Trace.STJ;
					tr.nodeidx[tpos] = 0;
					tr.pos[tpos] = 0;
				}

				else
					throw new RuntimeException("traceback failed");
				break;

			case P7Trace.STE: /* E connects from any M state. k set here */
				if (xmx[i][XME] <= -INFTY) {
					return null;
				}
				for (k = M; k >= 1; k--)
					if (xmx[i][XME] == mmx[i][k] + esc[k]) {
						/* check for wing unfolding */
						if (prob2Score(end[k], 1.) + 1 * INTSCALE <= esc[k]) {
							int dk; /*
									 * need a tmp k while moving thru delete
									 * wing
									 */
							for (dk = M; dk > k; dk--) {
								tr.statetype[tpos] = P7Trace.STD;
								tr.nodeidx[tpos] = dk;
								tr.pos[tpos] = 0;
								tpos++;
								if (tpos == curralloc) { /*
															 * grow trace if
															 * necessary
															 */
									curralloc += N;
									tr.P7ReallocTrace(curralloc);
								}
							}
						}

						tr.statetype[tpos] = P7Trace.STM;
						tr.nodeidx[tpos] = k--;
						tr.pos[tpos] = i--;
						break;
					}
				if (k < 0)
					throw new RuntimeException("traceback failed");
				break;

			case P7Trace.STC: /* C comes from C, E */
				if (xmx[i][XMC] <= -INFTY) {
					return null;
				} else if (xmx[i][XMC] == xmx[i - 1][XMC] + xsc[XTC][LOOP]) {
					tr.statetype[tpos] = P7Trace.STC;
					tr.nodeidx[tpos] = 0;
					tr.pos[tpos] = 0; /* note convention adherence: */
					tr.pos[tpos - 1] = i--; /* first C doesn't emit */
				} else if (xmx[i][XMC] == xmx[i][XME] + xsc[XTE][MOVE]) {
					tr.statetype[tpos] = P7Trace.STE;
					tr.nodeidx[tpos] = 0;
					tr.pos[tpos] = 0; /* E is a nonemitter */
				}

				else
					throw new RuntimeException("Traceback failed.");
				break;

			case P7Trace.STJ: /* J connects from E, J */
				if (xmx[i][XMJ] <= -INFTY) {
					return null;
				} else if (xmx[i][XMJ] == xmx[i - 1][XMJ] + xsc[XTJ][LOOP]) {
					tr.statetype[tpos] = P7Trace.STJ;
					tr.nodeidx[tpos] = 0;
					tr.pos[tpos] = 0; /* note convention adherence: */
					tr.pos[tpos - 1] = i--; /* first J doesn't emit */
				} else if (xmx[i][XMJ] == xmx[i][XME] + xsc[XTE][LOOP]) {
					tr.statetype[tpos] = P7Trace.STE;
					tr.nodeidx[tpos] = 0;
					tr.pos[tpos] = 0; /* E is a nonemitter */
				}

				else
					throw new RuntimeException("Traceback failed.");
				break;

			default:
				throw new RuntimeException("traceback failed");

			} /* end switch over statetype[tpos-1] */

			tpos++;
            
			if (tpos == curralloc) { /* grow trace if necessary */
				curralloc += N;
				tr.P7ReallocTrace(curralloc);
			}
		} /* end traceback, at S state; tpos == tlen now */
		tr.tlen = tpos;
		tr.P7ReverseTrace();
		// dsq is always 2 longer than actual sequence, since is stores other stuff in first and last cells.
		sc = xmx[dsq.length-2][XMC] + xsc[XTC][MOVE];
		tr.score = Scorify(sc);
		return tr;
	}

	/**
	 * Delete a state from hmm
	 * @param state to delete
	 * @param command the delelted node will be saved in command if not null
	 */
	public void deleteState(int d, HmmCommand command) {
		// modify necessary flags
		modified = true;
		M--;

		// generate new data structure
		// We use t.length and ins.length instead of M
		// since there length differs with hmmer format.
		double[][] newt = new double[t.length-1][7];
		double[][] newmat = new double[M + 1][alphabetSize];
		double[][] newins = new double[ins.length-1][alphabetSize];
		double[] newbegin = new double[M + 1];
		double[] newend = new double[M + 1];
		int[] newmap = new int[M + 1];
		int[][] newtsc = new int[7][M];
		int[][] newmsc = new int[MAXCODE][M + 1];
		int[][] newisc = new int[MAXCODE][M];
		int[] newbsc = new int[M + 1];
		int[] newesc = new int[M + 1];

		int[][] newdnam = new int[65][M + 1];
		int[][] newdnai = new int[65][M + 1];

		String[] newrf = new String[M+1]; // should be M+1 since we have one for every state
		String[] newcs = new String[M+1];

		// copy node upto d
		// start at 0 to include state 0 in hmmer3
		for (int x = 0; x < d; x++) {
			newmat[x] = mat[x];
			newbegin[x] = begin[x];
			newend[x] = end[x];

			newmap[x] = map[x];
			newrf[x] = rf[x];
			newcs[x] = cs[x];

			if (x < newt.length) {
				newt[x] = t[x];
				newins[x] = ins[x];
			}
		}

		// save node d to command history
		if (command != null) {
			command.oldState.mat = mat[d];
			command.oldState.begin = begin[d];
			command.oldState.end = end[d];

			if (d != M + 1) {
				command.oldState.t = t[d];
				command.oldState.ins = ins[d];

				command.oldState.rf = rf[d];
				command.oldState.cs = cs[d];
			} else {
				command.oldState.t = t[d - 1];
				command.oldState.ins = ins[d - 1];
				command.oldState.rf = rf[d - 1];
				command.oldState.cs = cs[d - 1];
			}
		}

		// skip node d, copy the rest
		for (int x = d; x <= M; x++) {
			newmat[x] = mat[x + 1];
			newbegin[x] = begin[x + 1];
			newend[x] = end[x + 1];

			newmap[x] = map[x+1];
			newrf[x] = rf[x + 1];
			newcs[x] = cs[x + 1];

			if (x < newt.length) {
				newt[x] = t[x + 1];
				newins[x] = ins[x + 1];
			}
		}

		if( d==1 && hmmerVersion==2) //removing first state in a hmmer2 model
		{
			newbegin[1]=begin[1];
		}
		else if( d==M+1 && hmmerVersion==2) //removing last state in a hmmer2 model
		{
			newend[M]=1;
		}
		else if( d==M+1 && hmmerVersion==3) //removing last state in a hmmer3 model
		{
			newt[M][TMD]=0;
			double norm= newt[M][TMM]+newt[M][TMI];
			newt[M][TMM] = newt[M][TMM] / norm;
			newt[M][TMI] = newt[M][TMI] / norm;

			newt[M][TDM]=1;
			newt[M][TDD]=0;
		}

		t = newt;
		mat = newmat;
		ins = newins;
		begin = newbegin;
		end = newend;
		map = newmap;
		tsc = newtsc;
		msc = newmsc;
		isc = newisc;
		bsc = newbsc;
		esc = newesc;

		dnam = newdnam;
		dnai = newdnai;

		rf = newrf;
		cs = newcs;
		getRepresentativeString();
	}

	/**
	 * Insert a state to hmm
	 * @param d the state to insert into 
	 * @param before if command is null, new state data will be copied from the previous state or from next state
	 * @param command if not null, the new node take the data inside command
	 */
	public void insertState(int d, boolean before, HmmCommand command) {
		// modify necessary flags
		modified = true;
		M++;

		// expand data structure
		double[][] newt = new double[t.length+1][7];
		double[][] newmat = new double[M + 1][alphabetSize];
		double[][] newins = new double[ins.length+1][alphabetSize];
		double[] newbegin = new double[M + 1];
		double[] newend = new double[M + 1];
		int[] newmap = new int[M + 1];
		int[][] newtsc = new int[7][M];
		int[][] newmsc = new int[MAXCODE][M + 1];
		int[][] newisc = new int[MAXCODE][M];
		int[] newbsc = new int[M + 1];
		int[] newesc = new int[M + 1];

		int[][] newdnam = new int[65][M + 1];
		int[][] newdnai = new int[65][M + 1];

		String[] newrf = new String[M+1];
		String[] newcs = new String[M+1];

		// copy from old data structure up to d, unchanged nodes
		for (int x = 0; x < d; x++) {
			newmat[x] = mat[x];
			newbegin[x] = begin[x];
			newend[x] = end[x];
			newmap[x] = map[x];
			for (int i = 0; i < MAXCODE; i++)
				newmsc[i][x] = msc[i][x];
			newbsc[x] = bsc[x];
			newesc[x] = esc[x];

			for (int i = 0; i < 65; i++) {
				newdnam[i][x] = dnam[i][x];
				newdnai[i][x] = dnai[i][x];
			}

			newrf[x] = rf[x];
			newcs[x] = cs[x];

			if (x  < t.length) {
				newt[x] = t[x];
				newins[x] = ins[x];
			}
			if (x  < tsc.length) {
				for (int i = 0; i < 7; i++)
					newtsc[i][x] = tsc[i][x];
				for (int i = 0; i < MAXCODE; i++)
					newisc[i][x] = isc[i][x];
			}
		}

		// new node, if command is not null (pass as a reference), use the value
		// in command, other wise,
		// copy from the node next to it, if "before" is true, copy from
		// previous node, otherwise, copy from next node
		if (command != null) {
			newmat[d] = command.oldState.mat;
			newbegin[d] = command.oldState.begin;
			newend[d] = command.oldState.end;

			if (d != M) {
				newt[d] = command.oldState.t;
				newins[d] = command.oldState.ins;

				newrf[d] = command.oldState.rf;
				newcs[d] = command.oldState.cs;
			} else { // inserting new last state
				newt[d - 1] = command.oldState.t;
				newins[d - 1] = command.oldState.ins;
				newrf[d - 1] = command.oldState.rf;
				newcs[d - 1] = command.oldState.cs;
			}
		} else {
			int copyFrom;
			if (before)
				copyFrom = d - 1;
			else
				copyFrom = d;

			newmat[d] = mat[copyFrom].clone();
			newrf[d] = rf[copyFrom];
			newcs[d] = cs[copyFrom];

			if (copyFrom == M - 1) // copying from last state
				copyFrom--;

			if (d != M) {
				newt[d] = t[copyFrom].clone();
				newins[d] = ins[copyFrom].clone();

				// these are different in the last state, so don't copy from there
				newbegin[d] = begin[copyFrom];
				newend[d] = end[copyFrom];

			}

			if (d == M && before) {
				// m->e on last state is 1, but should not be 1 on the second to last state.
				// also copy the b->m for consistancy
				newbegin[M-1] = begin[M-2];
				newend[M-1] = end[M-2];

				newt[M - 1] = t[M - 2].clone();
				newins[M - 1] = ins[M - 2].clone();

				if(M-1 < t.length)
				{
					newt[M] = t[M-1].clone();
					newins[M] = ins[M-1].clone();
				}

			}
		}

		// copy the rest nodes
		for (int x = d + 1; x <= M; x++) {
			newmat[x] = mat[x - 1];
			newbegin[x] = begin[x - 1];
			newend[x] = end[x - 1];
			newmap[x] = map[x-1];
			newrf[x] = rf[x - 1];
			newcs[x] = cs[x - 1];

			if (x-1  < t.length) {
				newt[x] = t[x - 1];
				newins[x] = ins[x - 1];
			}
		}

		// set the member variable to new data
		t = newt;
		mat = newmat;
		ins = newins;
		begin = newbegin;
		end = newend;
		map = newmap;
		tsc = newtsc;
		msc = newmsc;
		isc = newisc;
		bsc = newbsc;
		esc = newesc;

		dnam = newdnam;
		dnai = newdnai;

		rf = newrf;
		cs = newcs;

		// we add a node, so consensus string changes, easy to regenerate it
		getRepresentativeString();
	}


	/**
	 *  update sequence specific parameters of model. Only needed for HMMER3.
	 * @param L
	 */
	private void modifyTransitions(int L)
	{
		p1 = (double)L/L+1;
		xt[XTN][MOVE]=xt[XTC][MOVE]=xt[XTJ][MOVE]=2.0/(L+2);
		xt[XTN][LOOP]=xt[XTC][LOOP]=xt[XTJ][LOOP]=(double)L/(L+2); // 1-(2/(L+2))
	}
	  
	/**
	 * Align digitized sequence dsq to this model useing mx for storage.
	 * This is a more general version of doAlign to allow the user to precompute
	 * dsq and reuse the memory in mx.
	 * @param dsq a digitized sequece, as gotten from DigitizeSequence
	 * @param mx storage area, as gotten by  getNewPlan7Matrix
	 * @return trace of dsq
	 */
	public P7Trace doAlign(int[] dsq, DPMatrix mx) {  
		double saved1 = xt[XTE][MOVE];
		double saved2 = xt[XTE][LOOP];
		double tmmm=1,tmmi=0; //init values should not be used, just to satisfy javac

		if(hmmerVersion==3)
		{
			tmmm=t[M][TMM]; //save original values 
			tmmi=t[M][TMI];
			// In last state, M->D is always 0. We further
			// set M->M and M->I to zero to be consistant
			// with end[M] being 1.
			t[M][TMM]=t[M][TMI]=0.0; 
			modifyTransitions(dsq.length-2);
		}

		xt[XTE][MOVE] = 1.; /* only 1 domain/sequence ("global" alignment) */
		xt[XTE][LOOP] = 0.;

		P7Logoddsify(1);
		P7Trace tr = P7Viterbi(dsq, dsq.length-2, mx);

		xt[XTE][MOVE] = saved1;
		xt[XTE][LOOP] = saved2;

		if(hmmerVersion==3)
		{ //restore values so if we write to a file we don't lose the real values.
			t[M][TMM]=tmmm;
			t[M][TMI]=tmmi;
		}
		return tr;
	}
	/**
	 * Digitalize sequence, then create trace of that sequence
	 * @param seqs sequence to trace
	 * @return trace of the sequence
	 */
	public P7Trace doAlign(String seqs){
        
		int[] dsq = DigitizeSequence(seqs);
        return doAlign(dsq);
    }
	/**
	 * Align digitized sequence dsq to this model.
	 * This allows one to reuse a previously digitized sequence
	 * @param dsq a digitized sequece, as gotten from DigitizeSequence
	 * @return trace of dsq
	 */
	public P7Trace doAlign(int[] dsq) {
        //return doAlign(dsq,DPMatrix.CreatePlan7Matrix(1, M, 25, 0));
        return doAlign(dsq,DPMatrix.CreatePlan7Matrix(dsq.length, M, 25, 0));
    }
	/**
	 * Allocate storage used for running the alignment. This
	 * object can be reused for any sequence, but only with this model.
	 * This can save a lot of time by not re-allocateing memory.
	 * @return
	 */
    public  DPMatrix getNewPlan7Matrix()
    {
        return DPMatrix.CreatePlan7Matrix(1, M, 26, 1);
    }

	/**
	 * Align multiple sequence against the hmm
	 * @param sequences sequences to align
	 * @return MSA file structure for the align result
	 */
	public MSA alignAll(ArrayList<Sequence> sequences) {
		double saved1 = xt[XTE][MOVE];
		double saved2 = xt[XTE][LOOP];
		double tmmm=1,tmmi=0; //init values should not be used, just to satisfy javac

		if(hmmerVersion==3)
		{
			tmmm=t[M][TMM];
			tmmi=t[M][TMI];
			t[M][TMM]=t[M][TMI]=0.0;
		}

		xt[XTE][MOVE] = 1.; /* only 1 domain/sequence ("global" alignment) */
		xt[XTE][LOOP] = 0.;

		if(hmmerVersion==2)
			P7Logoddsify(1);
		P7Trace tr[] = new P7Trace[sequences.size()];
		int[][] dsq = new int[sequences.size()][];
		DPMatrix mx = DPMatrix.CreatePlan7Matrix(1, M, 25, 0);
		for (int i = 0; i < sequences.size(); i++) {
			dsq[i] = DigitizeSequence(sequences.get(i).seqs);

			if(hmmerVersion==3)
			{
				modifyTransitions(dsq[i].length-2);

				// this depends on p1 and xt, which in turn now
				//depend on the sequence length. So we must
				// regenerate it for each sequence
				P7Logoddsify(1);
			}

			tr[i] = P7Viterbi(dsq[i], sequences.get(i).seqs.length(), mx);
		}

		xt[XTE][MOVE] = saved1;
		xt[XTE][LOOP] = saved2;

		if(hmmerVersion==3)
		{
			t[M][TMM]=tmmm;
			t[M][TMI]=tmmi;
		}

		MSA msa = P7Trace.P7Traces2Alignment(dsq, sequences, M, tr, alphabet,
				false);
		return msa;
	}

	/**
	 * Creates a uniform distro for Nucleic sequences, and
	 * a pre-computed distro for Amino Acids. Hmmer3 uses
	 * these distros since it does not store the background distro
	 * in the file.
	 * @param d
	 */
	private void loadDefaultBackground(double[] d)
	{
		if(alphabetType != Atype.hmmNUCLEIC)
			for(int x=0; x < alphabetSize; x++)
				d[x]= 1.0/alphabetSize;
		else
		{
			d[0] = 0.0787945;     /* A */
			d[1] = 0.0151600;     /* C */
			d[2] = 0.0535222;     /* D */
			d[3] = 0.0668298;     /* E */
			d[4] = 0.0397062;     /* F */
			d[5] = 0.0695071;     /* G */
			d[6] = 0.0229198;     /* H */
			d[7] = 0.0590092;     /* I */
			d[8] = 0.0594422;     /* K */
			d[9] = 0.0963728;     /* L */
			d[10]= 0.0237718;     /* M */
			d[11]= 0.0414386;     /* N */
			d[12]= 0.0482904;     /* P */
			d[13]= 0.0395639;     /* Q */
			d[14]= 0.0540978;     /* R */
			d[15]= 0.0683364;     /* S */
			d[16]= 0.0540687;     /* T */
			d[17]= 0.0673417;     /* V */
			d[18]= 0.0114135;     /* W */
			d[19]= 0.0304133;     /* Y */
		}
	}
}
