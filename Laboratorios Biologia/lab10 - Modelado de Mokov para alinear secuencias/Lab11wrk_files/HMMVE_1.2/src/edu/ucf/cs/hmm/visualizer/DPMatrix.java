package edu.ucf.cs.hmm.visualizer;

/**
 * 
 * Adapted from Hmmer
 * Declaration of Plan7 dynamic programming matrix structure.
 *
 */
public class DPMatrix {
	int xmx[][]; /* special scores [0.1..N][BECJN] */

	int mmx[][]; /* match scores [0.1..N][0.1..M] */

	int imx[][]; /* insert scores [0.1..N][0.1..M-1.M] */

	int dmx[][]; /* delete scores [0.1..N][0.1..M-1.M] */

	/*
	 * Hidden ptrs where the real memory is kept; this trick was introduced by
	 * Erik Lindahl with the Altivec port; it's used to align xmx, etc. on
	 * 16-byte boundaries for cache optimization.
	 */
	// int xmx_mem[][];
	// int mmx_mem[][];
	// int imx_mem[][];
	// int dmx_mem[][];
	/*
	 * The other trick brought in w/ the Lindahl Altivec port; dp matrix is
	 * retained and grown, rather than reallocated for every HMM or sequence.
	 * Keep track of current allocated-for size in rows (sequence length N) and
	 * columns (HMM length M). Also keep track of pad sizes: how much we should
	 * overallocate rows or columns when we reallocate. If pad = 0, then we're
	 * not growable in this dimension.
	 */
	int maxN; /* alloc'ed for seq of length N; N+1 rows */

	int maxM; /* alloc'ed for HMM of length M; M+1 cols */

	int padN; /* extra pad in sequence length/rows */

	int padM; /* extra pad in HMM length/columns */

	/*
	 * Function: CreatePlan7Matrix()
	 * 
	 * Purpose: Create a dynamic programming matrix for standard Forward,
	 * Backward, or Viterbi, with scores kept as scaled log-odds integers. Keeps
	 * 2D arrays compact in RAM in an attempt to maximize cache hits.
	 * 
	 * The mx structure can be dynamically grown, if a new HMM or seq exceeds
	 * the currently allocated size. Dynamic growing is more efficient than an
	 * alloc/free of a whole matrix for every new target. The
	 * ResizePlan7Matrix() call does this reallocation, if needed. Here, in the
	 * creation step, we set up some pads - to inform the resizing call how much
	 * to overallocate when it realloc's.
	 * 
	 * Args: N - N+1 rows are allocated, for sequence. M - size of model in
	 * nodes padN - over-realloc in seq/row dimension, or 0 padM - over-realloc
	 * in HMM/column dimension, or 0
	 * 
	 * Return: mx mx is allocated here. Caller frees with FreePlan7Matrix(mx).
	 */
	static public DPMatrix CreatePlan7Matrix(int N, int M, int padN, int padM) {
		DPMatrix mx = new DPMatrix();

		mx.xmx = new int[N + 3][5];
		mx.mmx = new int[N + 1][M + 2];
		mx.imx = new int[N + 1][M + 2];
		mx.dmx = new int[N + 1][M + 2];
		// mx.xmx_mem = new int [N+1][5];
		// mx.mmx_mem = new int [N+1][M+2];
		// mx.imx_mem = new int [N+1][M+2];
		// mx.dmx_mem = new int [N+1][M+2];

		/*
		 * The indirect assignment below looks wasteful; it's actually used for
		 * aligning data on 16-byte boundaries as a cache optimization in the
		 * fast altivec implementation
		 */
		// mx.xmx[0] = mx.xmx_mem[0];
		// mx.mmx[0] = mx.mmx_mem[0];
		// mx.imx[0] = mx.imx_mem[0];
		// mx.dmx[0] = mx.dmx_mem[0];
		mx.maxN = N;
		mx.maxM = M;
		mx.padN = padN;
		mx.padM = padM;

		return mx;
	}

	/*
	 * Function: ResizePlan7Matrix()
	 * 
	 * Purpose: Reallocate a dynamic programming matrix, if necessary, for a
	 * problem of NxM: sequence length N, model size M. (N=1 for small memory
	 * score-only variants; we allocate N+1 rows in the DP matrix.)
	 * 
	 * We know (because of the way hmmsearch and hmmpfam are coded) that only
	 * one of the two dimensions is going to change in size after the first call
	 * to ResizePlan7Matrix(); that is, for hmmsearch, we have one HMM of fixed
	 * size M and our target sequences may grow in N; for hmmpfam, we have one
	 * sequence of fixed size N and our target models may grow in M. What we
	 * have to watch out for is P7SmallViterbi() working on a divide and conquer
	 * problem and passing us N < maxN, M > maxM; we should definitely *not*
	 * reallocate a smaller N. Since we know that only one dimension is going to
	 * grow, we aren't scared of reallocating to maxN,maxM. (If both M and N
	 * could grow, we would be more worried.)
	 * 
	 * Returns individual ptrs to the four matrix components as a convenience.
	 * 
	 * Args: mx - an already allocated model to grow. N - seq length to allocate
	 * for; N+1 rows M - size of model xmx, mmx, imx, dmx - RETURN: ptrs to four
	 * mx components as a convenience
	 * 
	 * Return: (void) mx is (re)allocated here.
	 */
	void ResizePlan7Matrix(int N, int M) {
        //System.out.println("resizing, N:"+N+", M: "+M+", maxN: "+maxN+", maxM: "+maxM);
		if (N <= maxN && M <= maxM)
			return;

        //System.out.println("resizing, N:"+N+", M: "+M);

		if (N > maxN || M > maxM) {
			N += padN;
			maxN = N;

			M += padM;
			maxM = M;

			xmx = new int[maxN + 3][5];
			mmx = new int[maxN + 1][maxM + 2];
			imx = new int[maxN + 1][maxM + 2];
			dmx = new int[maxN + 1][maxM + 2];
		}

		//if (M > maxM) {
			//M += padM;
			//maxM = M;
		//}

		// xmx_mem = new int[maxN+1][5];
		// mmx_mem = new int[maxN+1][maxM+2];
		// imx_mem = new int[maxN+1][maxM+2];
		// dmx_mem = new int[maxN+1][maxM+2];
		//
		// xmx[0] = xmx_mem[0];
		// mmx[0] = mmx_mem[0];
		// imx[0] = imx_mem[0];
		// dmx[0] = dmx_mem[0];
	}
}
