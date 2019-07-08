/**
 * 
 */
package edu.ucf.cs.hmm.visualizer;

import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;

/**
 * 
 * The node and line layout of the hmm
 *
 */
public class HmmLayout {
	int M;

	int size;

	int precision = 2;

	NodeState mNodes[];

	NodeState iNodes[];

	NodeState dNodes[];

	NodeState sNode;

	NodeState nNode;

	NodeState bNode;

	NodeState eNode;

	NodeState cNode;

	NodeState tNode;

	NodeState jNode;

	boolean showEmission = true;

	String[] emissionStrings; // emissionStrings[0]: emission string of N

	// emissionStrings[1..M]: emssion string of each stage
	// emissionStrings[M+1]: emission string of C

	// from m[x] to m[x+1]
	LineTransitionState[] mm;

	// from m[x] to i[x+1]
	LineTransitionState[] mi;

	// from m[x] to d[x+1]
	LineTransitionState[] md;

	// from i[x] to m[x+1]
	LineTransitionState[] im;

	// from i[x] to i[x] (self loop)
	LoopTransitionState[] ii;

	// from d[x] to m[x+1]
	LineTransitionState[] dm;

	// from d[x] to d[x+1]
	LineTransitionState[] dd;

	// from B to M
	BeginTransitionState[] bm;

	// from M to E
	EndTransitionState[] me;

	LineTransitionState sn;

	LoopTransitionState nn;

	LineTransitionState nb;

	LineTransitionState ec;

	LoopTransitionState cc;

	LineTransitionState ct;

	LoopTransitionState jj;

	EJTransitionState ej;

	JBTransitionState jb;

	LineTransitionState bd;

	LineTransitionState de;

	Dimension d = new Dimension();

	public HmmLayout(Hmm hmm, int _size) {
		// initialize memory structure
		size = _size;
		M = hmm.M;

		sNode = new NodeState(NodeState.TYPE_SPECIAL);
		nNode = new NodeState(NodeState.TYPE_LOOP);
		bNode = new NodeState(NodeState.TYPE_SPECIAL);

		mNodes = new NodeState[M + 1];
		iNodes = new NodeState[M];
		dNodes = new NodeState[M + 1];
		for (int i = 1; i <= M; i++) {
			mNodes[i] = new NodeState(NodeState.TYPE_M);
			if (i != M)
				iNodes[i] = new NodeState(NodeState.TYPE_I);
			dNodes[i] = new NodeState(NodeState.TYPE_D);
		}

		eNode = new NodeState(NodeState.TYPE_SPECIAL);
		cNode = new NodeState(NodeState.TYPE_LOOP);
		tNode = new NodeState(NodeState.TYPE_SPECIAL);
		jNode = new NodeState(NodeState.TYPE_LOOP);

		sn = new LineTransitionState();
		nn = new LoopTransitionState();
		nb = new LineTransitionState();
		ec = new LineTransitionState();
		cc = new LoopTransitionState();
		ct = new LineTransitionState();
		jj = new LoopTransitionState();
		ej = new EJTransitionState();
		jb = new JBTransitionState();
		bd = new LineTransitionState();
		de = new LineTransitionState();

		mm = new LineTransitionState[M];
		mi = new LineTransitionState[M];
		md = new LineTransitionState[M];
		im = new LineTransitionState[M];
		ii = new LoopTransitionState[M];
		dm = new LineTransitionState[M];
		dd = new LineTransitionState[M];
		bm = new BeginTransitionState[M + 1];
		me = new EndTransitionState[M + 1];
		for (int i = 1; i <= M; i++) {
			if (i != M) {
				mm[i] = new LineTransitionState();
				mi[i] = new LineTransitionState();
				md[i] = new LineTransitionState();
				im[i] = new LineTransitionState();
				ii[i] = new LoopTransitionState();
				dm[i] = new LineTransitionState();
				dd[i] = new LineTransitionState();
			}
			bm[i] = new BeginTransitionState();
			me[i] = new EndTransitionState();
		}
		clearPath();
		setHmmValue(hmm);
		defaultLayout();
		reNumber();
		generateShape(_size);
	}

	/**
	 * copy the probability of hmm model, will use this to determine thickness of the line and label on the line
	 * @param hmm
	 */
	void setHmmValue(Hmm hmm) {
		for (int i = 1; i <= M; i++) {
			mNodes[i].consensus = hmm.consensus.substring(i - 1, i);
			mNodes[i].p = hmm.consensusP[i];
			if (i != M) {
				mm[i].p = hmm.t[i][0];
				mi[i].p = hmm.t[i][1];
				md[i].p = hmm.t[i][2];
				im[i].p = hmm.t[i][3];
				ii[i].p = hmm.t[i][4];
				dm[i].p = hmm.t[i][5];
				dd[i].p = hmm.t[i][6];
			}
			bm[i].p = hmm.begin[i];
			me[i].p = hmm.end[i];
		}
		sNode.label = "S";
		nNode.label = "N";
		bNode.label = "B";
		eNode.label = "E";
		cNode.label = "C";
		tNode.label = "T";
		jNode.label = "J";
		sn.p = 1.0;
		nn.p = hmm.xt[Hmm.XTN][Hmm.LOOP];  // these indices are backwards
		nb.p = hmm.xt[Hmm.XTN][Hmm.MOVE];
		bd.p = hmm.tbd1;
		de.p = 1.0;
		ec.p = hmm.xt[Hmm.XTE][Hmm.MOVE];
		cc.p = hmm.xt[Hmm.XTC][Hmm.LOOP];
		ct.p = hmm.xt[Hmm.XTC][Hmm.MOVE];
		ej.p = hmm.xt[Hmm.XTE][Hmm.LOOP];
		jj.p = hmm.xt[Hmm.XTJ][Hmm.LOOP];
		jb.p = hmm.xt[Hmm.XTJ][Hmm.MOVE];
	}

	/**
	 * Clear Viterbi path
	 *
	 */
	void clearPath() {
		for (int i = 1; i <= M; i++) {
			mNodes[i].inPath = false;
			if (i != M)
				iNodes[i].inPath = false;
			dNodes[i].inPath = false;
			if (i != M) {
				mm[i].inPath = false;
				mi[i].inPath = false;
				md[i].inPath = false;
				im[i].inPath = false;
				ii[i].inPath = false;
				dm[i].inPath = false;
				dd[i].inPath = false;
			}
			bm[i].inPath = false;
			me[i].inPath = false;
		}
		sNode.inPath = false;
		nNode.inPath = false;
		bNode.inPath = false;
		eNode.inPath = false;
		cNode.inPath = false;
		tNode.inPath = false;
		jNode.inPath = false;
		sn.inPath = false;
		nb.inPath = false;
		nn.inPath = false;
		bd.inPath = false;
		de.inPath = false;
		ec.inPath = false;
		ej.inPath = false;
		cc.inPath = false;
		ct.inPath = false;
		jb.inPath = false;
		jj.inPath = false;
		emissionStrings = null;
	}

	Dimension getDimension() {
		return d;
	}

	/**
	 * Draw the hmm layout
	 * @param g2d graphics context
	 * @param clip clip of the displayable region, clip for time saving
	 */
	void drawWithLabel(Graphics2D g2d, Rectangle2D clip) {
		// draw all connection
		for (int i = 1; i <= M; i++) {
			if (i != M) {
				mm[i].draw(g2d, clip);
				mi[i].draw(g2d, clip);
				md[i].draw(g2d, clip);
				dd[i].draw(g2d, clip);
				im[i].draw(g2d, clip);
				dm[i].draw(g2d, clip);
				ii[i].draw(g2d, clip);
			}

			bm[i].draw(g2d, clip);
			me[i].draw(g2d, clip);
		}
		nn.draw(g2d, clip);
		cc.draw(g2d, clip);
		jj.draw(g2d, clip);
		sn.draw(g2d, clip);
		nb.draw(g2d, clip);
		ec.draw(g2d, clip);
		ct.draw(g2d, clip);
		ej.draw(g2d, clip);
		jb.draw(g2d, clip);
		bd.draw(g2d, clip);
		de.draw(g2d, clip);

		// draw all node
		sNode.draw(g2d, clip);
		nNode.draw(g2d, clip);
		bNode.draw(g2d, clip);

		for (int i = 1; i <= M; i++) {
			mNodes[i].draw(g2d, clip);
			if (i != M)
				iNodes[i].draw(g2d, clip);
			dNodes[i].draw(g2d, clip);
		}
		eNode.draw(g2d, clip);
		cNode.draw(g2d, clip);
		tNode.draw(g2d, clip);
		jNode.draw(g2d, clip);

		// draw all text
		for (int i = 1; i <= M; i++) {
			if (i != M) {
				mm[i].drawText(g2d, clip, precision);
				mi[i].drawText(g2d, clip, precision);
				md[i].drawText(g2d, clip, precision);
				dd[i].drawText(g2d, clip, precision);
				im[i].drawText(g2d, clip, precision);
				dm[i].drawText(g2d, clip, precision);
				ii[i].drawText(g2d, clip, precision);
			}
		}
		nn.drawText(g2d, clip, precision);
		cc.drawText(g2d, clip, precision);
		jj.drawText(g2d, clip, precision);
		sn.drawText(g2d, clip, precision);
		nb.drawText(g2d, clip, precision);
		ec.drawText(g2d, clip, precision);
		ct.drawText(g2d, clip, precision);
		ej.drawText(g2d, clip, precision);
		jb.drawText(g2d, clip, precision);
		bd.drawText(g2d, clip, precision);
		de.drawText(g2d, clip, precision);

		sNode.drawText(g2d, clip);
		nNode.drawText(g2d, clip);
		bNode.drawText(g2d, clip);

		for (int i = 1; i <= M; i++) {
			mNodes[i].drawText(g2d, clip);
			if (i != M)
				iNodes[i].drawText(g2d, clip);
			dNodes[i].drawText(g2d, clip);
		}
		eNode.drawText(g2d, clip);
		cNode.drawText(g2d, clip);
		tNode.drawText(g2d, clip);
		jNode.drawText(g2d, clip);

		if (showEmission && emissionStrings != null) {
			FontMetrics metrics = g2d.getFontMetrics();
			int txtHeight = metrics.getHeight();
			int txtWidth;

			txtWidth = metrics.stringWidth(emissionStrings[0]);
			g2d.drawString(emissionStrings[0], nNode.position.x - txtWidth / 2,
					txtHeight);

			for (int i = 1; i <= M; i++) {
				if (emissionStrings[i].startsWith("-")) {
					String s = emissionStrings[i];
					if (s.length() > 3)
						s = s.substring(0, 1) + "..."
								+ s.substring(s.length() - 1) + "["
								+ s.length() + "]";
					txtWidth = metrics.stringWidth(s);
					g2d.drawString(s, dNodes[i].position.x - txtWidth / 2,
							txtHeight);
				} else {
					String s1 = null;
					String s2 = null;
					if (emissionStrings[i].indexOf("(") != -1) {
						s1 = emissionStrings[i].substring(0, emissionStrings[i]
								.indexOf("("));
						s2 = emissionStrings[i].substring(emissionStrings[i]
								.indexOf("("));
					} else
						s1 = emissionStrings[i];

					if (s1 != null && !s1.equals("")) {
						if (s1.length() > 3)
							s1 = s1.substring(0, 1) + "..."
									+ s1.substring(s1.length() - 1) + "("
									+ s1.length() + ")";
						txtWidth = metrics.stringWidth(s1);
						g2d.drawString(s1, mNodes[i].position.x - txtWidth / 2,
								txtHeight);
					}

					if (s2 != null && !s2.equals("")) {
						if (s2.length() > 3 * 3) // 3 alphabet, every one in
													// the form (x)
							s2 = s2.substring(0, 3) + "..."
									+ s2.substring(s2.length() - 3) + "["
									+ s2.length() / 3 + "]";
						txtWidth = metrics.stringWidth(s2);
						g2d.drawString(s2, iNodes[i].position.x - txtWidth / 2,
								txtHeight);
					}

				}
			}

			txtWidth = metrics.stringWidth(emissionStrings[M + 1]);
			g2d.drawString(emissionStrings[M + 1], cNode.position.x - txtWidth
					/ 2, txtHeight);
		}
	}

	/**
	 * Response to zoom operation
	 * @param factor new size = old size * factor
	 */
	void resize(double factor) {
		size = (int) Math.ceil((double) size * factor);
		calculateDimension();
		for (int i = 1; i <= M; i++) {
			if (i != M) {
				mm[i].start.x *= factor;
				mm[i].start.y *= factor;
				mm[i].end.x *= factor;
				mm[i].end.y *= factor;
				mi[i].start.x *= factor;
				mi[i].start.y *= factor;
				mi[i].end.x *= factor;
				mi[i].end.y *= factor;
				md[i].start.x *= factor;
				md[i].start.y *= factor;
				md[i].end.x *= factor;
				md[i].end.y *= factor;
				dd[i].start.x *= factor;
				dd[i].start.y *= factor;
				dd[i].end.x *= factor;
				dd[i].end.y *= factor;
				im[i].start.x *= factor;
				im[i].start.y *= factor;
				im[i].end.x *= factor;
				im[i].end.y *= factor;
				dm[i].start.x *= factor;
				dm[i].start.y *= factor;
				dm[i].end.x *= factor;
				dm[i].end.y *= factor;
				ii[i].center.x *= factor;
				ii[i].center.y *= factor;
			}

			bm[i].start.x *= factor;
			bm[i].start.y *= factor;
			bm[i].point1.x *= factor;
			bm[i].point1.y *= factor;
			bm[i].end.x *= factor;
			bm[i].end.y *= factor;
			me[i].start.x *= factor;
			me[i].start.y *= factor;
			me[i].point1.x *= factor;
			me[i].point1.y *= factor;
			me[i].end.x *= factor;
			me[i].end.y *= factor;
		}
		nn.center.x *= factor;
		nn.center.y *= factor;
		cc.center.x *= factor;
		cc.center.y *= factor;
		jj.center.x *= factor;
		jj.center.y *= factor;
		sn.start.x *= factor;
		sn.start.y *= factor;
		sn.end.x *= factor;
		sn.end.y *= factor;
		nb.start.x *= factor;
		nb.start.y *= factor;
		nb.end.x *= factor;
		nb.end.y *= factor;
		ec.start.x *= factor;
		ec.start.y *= factor;
		ec.end.x *= factor;
		ec.end.y *= factor;
		ct.start.x *= factor;
		ct.start.y *= factor;
		ct.end.x *= factor;
		ct.end.y *= factor;
		ej.start.x *= factor;
		ej.start.y *= factor;
		ej.point1.x *= factor;
		ej.point1.y *= factor;
		ej.point2.x *= factor;
		ej.point2.y *= factor;
		ej.end.x *= factor;
		ej.end.y *= factor;
		jb.start.x *= factor;
		jb.start.y *= factor;
		jb.point1.x *= factor;
		jb.point1.y *= factor;
		jb.point2.x *= factor;
		jb.point2.y *= factor;
		jb.end.x *= factor;
		jb.end.y *= factor;
		bd.start.x *= factor;
		bd.start.y *= factor;
		bd.end.x *= factor;
		bd.end.y *= factor;
		de.start.x *= factor;
		de.start.y *= factor;
		de.end.x *= factor;
		de.end.y *= factor;

		sNode.position.x *= factor;
		sNode.position.y *= factor;
		nNode.position.x *= factor;
		nNode.position.y *= factor;
		bNode.position.x *= factor;
		bNode.position.y *= factor;

		for (int i = 1; i <= M; i++) {
			mNodes[i].position.x *= factor;
			mNodes[i].position.y *= factor;
			if (i != M) {
				iNodes[i].position.x *= factor;
				iNodes[i].position.y *= factor;
			}
			dNodes[i].position.x *= factor;
			dNodes[i].position.y *= factor;
		}
		eNode.position.x *= factor;
		eNode.position.y *= factor;
		cNode.position.x *= factor;
		cNode.position.y *= factor;
		tNode.position.x *= factor;
		tNode.position.y *= factor;
		jNode.position.x *= factor;
		jNode.position.y *= factor;

		generateShape(size);
	}

	/**
	 * Generate shape of all state and attached lines and arrows
	 * @param _size
	 */
	void generateShape(int _size) {
		size = _size;
		sNode.generateShape(size);
		nNode.generateShape(size);
		bNode.generateShape(size);

		for (int i = 1; i <= M; i++) {
			mNodes[i].generateShape(size);
			if (i != M)
				iNodes[i].generateShape(size);
			dNodes[i].generateShape(size);
		}
		eNode.generateShape(size);
		cNode.generateShape(size);
		tNode.generateShape(size);
		jNode.generateShape(size);

		for (int i = 1; i <= M; i++) {
			if (i != M) {
				mm[i].generateShape(size);
				mi[i].generateShape(size);
				md[i].generateShape(size);
				dd[i].generateShape(size);
				im[i].generateShape(size);
				dm[i].generateShape(size);
				ii[i].generateShape(size);
			}

			bm[i].generateShape(size);
			me[i].generateShape(size);
		}
		nn.generateShape(size);
		cc.generateShape(size);
		jj.generateShape(size);
		sn.generateShape(size);
		nb.generateShape(size);
		ec.generateShape(size);
		ct.generateShape(size);
		ej.generateShape(size);
		jb.generateShape(size);
		bd.generateShape(size);
		de.generateShape(size);
	}

	void calculateDimension() {
		d.width = (M + 7) * size * 4;
		d.height = size * 16;
	}

	/**
	 * After insertion and deletion, the state number will mess up, use renumber to relabel them to 1..N
	 *
	 */
	void reNumber() {
		for (int i = 1; i <= M; i++) {
			mNodes[i].label = "M" + i;
			dNodes[i].label = "D" + i;
			if (i != M)
				iNodes[i].label = "I" + i;
		}
	}

	/**
	 * Rearrange all states on the screen after drag all nodes around
	 *
	 */
	void defaultLayout() {
		calculateDimension();

		// node layout
		sNode.position.x = size * 2;
		sNode.position.y = size * 9;
		nNode.position.x = size * 6;
		nNode.position.y = size * 9;
		bNode.position.x = size * 10;
		bNode.position.y = size * 9;
		for (int i = 1; i <= M; i++) {
			mNodes[i].position.x = size * 10 + size * 4 * i;
			mNodes[i].position.y = size * 7;
			if (i != M) {
				iNodes[i].position.x = size * 10 + size * 4 * i + size * 2;
				iNodes[i].position.y = size * 3;
			}
			dNodes[i].position.x = size * 10 + size * 4 * i;
			dNodes[i].position.y = size * 11;
		}
		eNode.position.x = size * 10 + size * 4 * M + size * 4;
		eNode.position.y = size * 5;
		cNode.position.x = size * 10 + size * 4 * M + size * 8;
		cNode.position.y = size * 5;
		tNode.position.x = size * 10 + size * 4 * M + size * 12;
		tNode.position.y = size * 5;
		jNode.position.x = d.width / 2;
		jNode.position.y = size * 14;

		makeConnect();
	}

	/**
	 * Adjust the start point and end point of each line to connect nodes
	 *
	 */
	void makeConnect() {
		// transition layout
		for (int i = 1; i <= M; i++) {
			if (i != M) {
				mm[i].start.x = mNodes[i].position.x + size;
				mm[i].start.y = mNodes[i].position.y;
				mm[i].end.x = mNodes[i + 1].position.x - size;
				mm[i].end.y = mNodes[i + 1].position.y;

				mi[i].start.x = mNodes[i].position.x + size;
				mi[i].start.y = mNodes[i].position.y;
				mi[i].end.x = iNodes[i].position.x;
				mi[i].end.y = iNodes[i].position.y + size;
				mi[i].labelPosition = TransitionState.POSITION_START;

				md[i].start.x = mNodes[i].position.x + size;
				md[i].start.y = mNodes[i].position.y;
				md[i].end.x = dNodes[i + 1].position.x - size;
				md[i].end.y = dNodes[i + 1].position.y;
				md[i].labelPosition = TransitionState.POSITION_START;

				dd[i].start.x = dNodes[i].position.x + size;
				dd[i].start.y = dNodes[i].position.y;
				dd[i].end.x = dNodes[i + 1].position.x - size;
				dd[i].end.y = dNodes[i + 1].position.y;

				im[i].start.x = iNodes[i].position.x;
				im[i].start.y = iNodes[i].position.y + size;
				im[i].end.x = mNodes[i + 1].position.x - size;
				im[i].end.y = mNodes[i + 1].position.y;
				im[i].labelPosition = TransitionState.POSITION_START;

				dm[i].start.x = dNodes[i].position.x + size;
				dm[i].start.y = dNodes[i].position.y;
				dm[i].end.x = mNodes[i + 1].position.x - size;
				dm[i].end.y = mNodes[i + 1].position.y;
				dm[i].labelPosition = TransitionState.POSITION_START;

				ii[i].center.x = iNodes[i].position.x;
				ii[i].center.y = iNodes[i].position.y;
			}

			bm[i].start.x = bNode.position.x + size;
			bm[i].start.y = bNode.position.y;
			bm[i].end.x = mNodes[i].position.x - size;
			bm[i].end.y = mNodes[i].position.y;

			me[i].start.x = mNodes[i].position.x + size;
			me[i].start.y = mNodes[i].position.y;
			me[i].end.x = eNode.position.x - size;
			me[i].end.y = eNode.position.y;
		}

		nn.center.x = nNode.position.x;
		nn.center.y = nNode.position.y;

		cc.center.x = cNode.position.x;
		cc.center.y = cNode.position.y;

		jj.center.x = jNode.position.x;
		jj.center.y = jNode.position.y;

		sn.start.x = sNode.position.x + size;
		sn.start.y = sNode.position.y;
		sn.end.x = nNode.position.x - size;
		sn.end.y = nNode.position.y;

		nb.start.x = nNode.position.x + size;
		nb.start.y = nNode.position.y;
		nb.end.x = bNode.position.x - size;
		nb.end.y = bNode.position.y;

		ec.start.x = eNode.position.x + size;
		ec.start.y = eNode.position.y;
		ec.end.x = cNode.position.x - size;
		ec.end.y = cNode.position.y;

		ct.start.x = cNode.position.x + size;
		ct.start.y = cNode.position.y;
		ct.end.x = tNode.position.x - size;
		ct.end.y = tNode.position.y;

		ej.start.x = eNode.position.x + size;
		ej.start.y = eNode.position.y;
		ej.end.x = jNode.position.x + size;
		ej.end.y = jNode.position.y;

		jb.start.x = jNode.position.x - size;
		jb.start.y = jNode.position.y;
		jb.end.x = bNode.position.x - size;
		jb.end.y = bNode.position.y;

		bd.start.x = bNode.position.x + size;
		bd.start.y = bNode.position.y;
		bd.end.x = dNodes[1].position.x - size;
		bd.end.y = dNodes[1].position.y;

		de.start.x = dNodes[M].position.x + size;
		de.start.y = dNodes[M].position.y;
		de.end.x = eNode.position.x - size;
		de.end.y = eNode.position.y;
	}

	/**
	 * Mark the path taken by HMM
	 * @param tr trace of the sequence
	 * @param seqs sequence (display on top of layout)
	 */
	void tagPath(P7Trace tr, String seqs) {
		clearPath();
		emissionStrings = new String[M + 2];
		for (int i = 0; i < M + 2; i++)
			emissionStrings[i] = "";
		int lastState = P7Trace.STBOGUS;
		int state = 0;
		int seq = 0; // track the position in seqs
		for (int i = 0; i < tr.statetype.length; i++) {
			switch (tr.statetype[i]) {
			case P7Trace.STS:
				sNode.inPath = true;
				break;
			case P7Trace.STN:
				nNode.inPath = true;
				if (lastState == P7Trace.STS) {
					sn.inPath = true;
				}
				if (lastState == P7Trace.STN) {
					emissionStrings[0] += seqs.charAt(seq);
					seq++;
					nn.inPath = true;
				}
				break;
			case P7Trace.STB:
				bNode.inPath = true;
				if (lastState == P7Trace.STN) {
					nb.inPath = true;
				}
				state = 0;
				break;
			case P7Trace.STM:
				state++;
				mNodes[state].inPath = true;
				if (lastState == P7Trace.STB) {
					bm[state].inPath = true;
				}
				if (lastState == P7Trace.STM) {
					mm[state - 1].inPath = true;
				}
				if (lastState == P7Trace.STD) {
					dm[state - 1].inPath = true;
				}
				if (lastState == P7Trace.STI) {
					im[state - 1].inPath = true;
				}
				emissionStrings[state] += seqs.charAt(seq);
				seq++;
				break;
			case P7Trace.STI:
				if (lastState == P7Trace.STM) {
					iNodes[state].inPath = true;
					mi[state].inPath = true;
				}
				if (lastState == P7Trace.STI) {
					ii[state].inPath = true;
				}
				emissionStrings[state] += "(" + seqs.charAt(seq) + ")";
				seq++;
				break;
			case P7Trace.STD:
				state++;
				dNodes[state].inPath = true;
				if (lastState == P7Trace.STB) {
					bd.inPath = true;
				}
				if (lastState == P7Trace.STM) {
					md[state - 1].inPath = true;
				}
				if (lastState == P7Trace.STD) {
					dd[state - 1].inPath = true;
				}
				emissionStrings[state] += "-";
				break;
			case P7Trace.STE:
				eNode.inPath = true;
				if (lastState == P7Trace.STM) {
					me[state].inPath = true;
				}
				if (lastState == P7Trace.STD) {
					de.inPath = true;
				}
				break;
			case P7Trace.STC:
				cNode.inPath = true;
				if (lastState == P7Trace.STE) {
					ec.inPath = true;
				}
				if (lastState == P7Trace.STC) {
					cc.inPath = true;
					emissionStrings[M + 1] += seqs.charAt(seq);
					seq++;
				}
				break;
			case P7Trace.STT:
				tNode.inPath = true;
				if (lastState == P7Trace.STC) {
					ct.inPath = true;
				}
				break;
			}
			lastState = tr.statetype[i];
		}
	}

	/**
	 * Move a node to a new position
	 * @param s node to move
	 * @param newx
	 * @param newy
	 */
	void moveNode(NodeState s, int newx, int newy) {
		s.position.x = newx;
		s.position.y = newy;
		s.generateShape(size);
		if (s.type == NodeState.TYPE_M) {
			for (int i = 1; i <= M; i++) {
				if (s == mNodes[i]) {
					if (i != 1) {
						mm[i - 1].end.x = mNodes[i].position.x - size;
						mm[i - 1].end.y = mNodes[i].position.y;
						mm[i - 1].generateShape(size);
						im[i - 1].end.x = mNodes[i].position.x - size;
						im[i - 1].end.y = mNodes[i].position.y;
						im[i - 1].generateShape(size);
						dm[i - 1].end.x = mNodes[i].position.x - size;
						dm[i - 1].end.y = mNodes[i].position.y;
						dm[i - 1].generateShape(size);
					}
					if (i != M) {
						mm[i].start.x = mNodes[i].position.x + size;
						mm[i].start.y = mNodes[i].position.y;
						mm[i].generateShape(size);
						mi[i].start.x = mNodes[i].position.x + size;
						mi[i].start.y = mNodes[i].position.y;
						mi[i].generateShape(size);
						md[i].start.x = mNodes[i].position.x + size;
						md[i].start.y = mNodes[i].position.y;
						md[i].generateShape(size);
					}
					bm[i].end.x = mNodes[i].position.x - size;
					bm[i].end.y = mNodes[i].position.y;
					bm[i].generateShape(size);
					me[i].start.x = mNodes[i].position.x + size;
					me[i].start.y = mNodes[i].position.y;
					me[i].generateShape(size);
				}
			}
		} else if (s.type == NodeState.TYPE_I) {
			for (int i = 1; i <= M - 1; i++) {
				if (s == iNodes[i]) {
					mi[i].end.x = iNodes[i].position.x;
					mi[i].end.y = iNodes[i].position.y + size;
					mi[i].generateShape(size);

					im[i].start.x = iNodes[i].position.x;
					im[i].start.y = iNodes[i].position.y + size;
					im[i].generateShape(size);

					ii[i].center.x = iNodes[i].position.x;
					ii[i].center.y = iNodes[i].position.y;
					ii[i].generateShape(size);
				}
			}
		} else if (s.type == NodeState.TYPE_D) {
			for (int i = 1; i <= M; i++) {
				if (s == dNodes[i]) {
					if (i != 1) {
						md[i - 1].end.x = dNodes[i].position.x - size;
						md[i - 1].end.y = dNodes[i].position.y;
						md[i - 1].generateShape(size);

						dd[i - 1].end.x = dNodes[i].position.x - size;
						dd[i - 1].end.y = dNodes[i].position.y;
						dd[i - 1].generateShape(size);
					}
					if (i != M) {
						dm[i].start.x = dNodes[i].position.x + size;
						dm[i].start.y = dNodes[i].position.y;
						dm[i].generateShape(size);

						dd[i].start.x = dNodes[i].position.x + size;
						dd[i].start.y = dNodes[i].position.y;
						dd[i].generateShape(size);
					}
					if (i == 1) {
						bd.end.x = dNodes[1].position.x - size;
						bd.end.y = dNodes[1].position.y;
						bd.generateShape(size);
					}
					if (i == M) {
						de.start.x = dNodes[M].position.x + size;
						de.start.y = dNodes[M].position.y;
						de.generateShape(size);
					}
				}
			}
		} else if (s.type == NodeState.TYPE_SPECIAL
				|| s.type == NodeState.TYPE_LOOP) {
			if (s == sNode) {
				sn.start.x = sNode.position.x + size;
				sn.start.y = sNode.position.y;
				sn.generateShape(size);
			} else if (s == nNode) {
				sn.end.x = nNode.position.x - size;
				sn.end.y = nNode.position.y;
				sn.generateShape(size);
				nn.center.x = nNode.position.x;
				nn.center.y = nNode.position.y;
				nn.generateShape(size);
				nb.start.x = nNode.position.x + size;
				nb.start.y = nNode.position.y;
				nb.generateShape(size);
			} else if (s == bNode) {
				bd.start.x = bNode.position.x + size;
				bd.start.y = bNode.position.y;
				bd.generateShape(size);
				nb.end.x = bNode.position.x - size;
				nb.end.y = bNode.position.y;
				nb.generateShape(size);
				jb.end.x = bNode.position.x - size;
				jb.end.y = bNode.position.y;
				jb.generateShape(size);
				for (int i = 1; i <= M; i++) {
					bm[i].start.x = bNode.position.x + size;
					bm[i].start.y = bNode.position.y;
					bm[i].generateShape(size);
				}
			} else if (s == eNode) {
				de.end.x = eNode.position.x - size;
				de.end.y = eNode.position.y;
				de.generateShape(size);
				ec.start.x = eNode.position.x + size;
				ec.start.y = eNode.position.y;
				ec.generateShape(size);
				ej.start.x = eNode.position.x + size;
				ej.start.y = eNode.position.y;
				ej.generateShape(size);
				for (int i = 1; i <= M; i++) {
					me[i].end.x = eNode.position.x - size;
					me[i].end.y = eNode.position.y;
					me[i].generateShape(size);
				}
			} else if (s == cNode) {
				ec.end.x = cNode.position.x - size;
				ec.end.y = cNode.position.y;
				ec.generateShape(size);
				cc.center.x = cNode.position.x;
				cc.center.y = cNode.position.y;
				cc.generateShape(size);
				ct.start.x = cNode.position.x;
				ct.start.y = cNode.position.y;
				ct.generateShape(size);
			} else if (s == jNode) {
				ej.end.x = jNode.position.x + size;
				ej.end.y = jNode.position.y;
				ej.generateShape(size);
				jj.center.x = jNode.position.x;
				jj.center.y = jNode.position.y;
				jj.generateShape(size);
				jb.start.x = jNode.position.x - size;
				jb.start.y = jNode.position.y;
				jb.generateShape(size);
			} else if (s == tNode) {
				ct.end.x = tNode.position.x - size;
				ct.end.y = tNode.position.y;
				ct.generateShape(size);
			}
		}
	}

	/**
	 * Delete a State in layout, a state include one matching node, one deletion node and one insertion node
	 * @param d which state to delete
	 * @param command the delelted node will be saved in command if not null
	 */
	void deleteState(int d, HmmCommand command) {
		clearPath();
		M--;

		// generate new data structure
		NodeState newmNodes[] = new NodeState[M + 1];
		NodeState newiNodes[] = new NodeState[M];
		NodeState newdNodes[] = new NodeState[M + 1];

		LineTransitionState[] newmm = new LineTransitionState[M];
		LineTransitionState[] newmi = new LineTransitionState[M];
		LineTransitionState[] newmd = new LineTransitionState[M];
		LineTransitionState[] newim = new LineTransitionState[M];
		LoopTransitionState[] newii = new LoopTransitionState[M];
		LineTransitionState[] newdm = new LineTransitionState[M];
		LineTransitionState[] newdd = new LineTransitionState[M];
		BeginTransitionState[] newbm = new BeginTransitionState[M + 1];
		EndTransitionState[] newme = new EndTransitionState[M + 1];

		// copy from the old data structure from node d (node to delete)
		for (int i = 1; i < d; i++) {
			newmNodes[i] = mNodes[i];
			newdNodes[i] = dNodes[i];
			newbm[i] = bm[i];
			newme[i] = me[i];
			if (i != M) // for the last node after deletion, don't have
						// following structure, so do not copy
			{
				newiNodes[i] = iNodes[i];
				newmm[i] = mm[i];
				newmi[i] = mi[i];
				newmd[i] = md[i];
				newim[i] = im[i];
				newii[i] = ii[i];
				newdm[i] = dm[i];
				newdd[i] = dd[i];
			}
		}

		// let's save the deleted node to command history
		if (command != null) {
			command.oldState.mNodes = mNodes[d];
			command.oldState.dNodes = dNodes[d];

			command.oldState.bm = bm[d];
			command.oldState.me = me[d];
			if (d != M + 1) // for last node before deletion, we don't have
							// following structure; on the other hand, the
			// according structure of the previous node will deleted also, so we
			// keep the data structure of the previous node
			{
				command.oldState.iNodes = iNodes[d];
				command.oldState.mm = mm[d];
				command.oldState.mi = mi[d];
				command.oldState.md = md[d];
				command.oldState.im = im[d];
				command.oldState.ii = ii[d];
				command.oldState.dm = dm[d];
				command.oldState.dd = dd[d];
			} else {
				command.oldState.iNodes = iNodes[d - 1];
				command.oldState.mm = mm[d - 1];
				command.oldState.mi = mi[d - 1];
				command.oldState.md = md[d - 1];
				command.oldState.im = im[d - 1];
				command.oldState.ii = ii[d - 1];
				command.oldState.dm = dm[d - 1];
				command.oldState.dd = dd[d - 1];
			}
		}

		// copy the rest nodes
		for (int i = d; i <= M; i++) {
			newmNodes[i] = mNodes[i + 1];
			newdNodes[i] = dNodes[i + 1];
			newbm[i] = bm[i + 1];
			newme[i] = me[i + 1];
			if (i != M) {
				newiNodes[i] = iNodes[i + 1];
				newmm[i] = mm[i + 1];
				newmi[i] = mi[i + 1];
				newmd[i] = md[i + 1];
				newim[i] = im[i + 1];
				newii[i] = ii[i + 1];
				newdm[i] = dm[i + 1];
				newdd[i] = dd[i + 1];
			}
		}
		mNodes = newmNodes;
		iNodes = newiNodes;
		dNodes = newdNodes;

		mm = newmm;
		mi = newmi;
		md = newmd;
		im = newim;
		ii = newii;
		dm = newdm;
		dd = newdd;
		bm = newbm;
		me = newme;

		makeConnect();
	}

	/**
	 * Insert a state to layout, include a matching node, a deletion node and an insertion node
	 * @param d the state to insert into 
	 * @param before if command is null, new state data will be copied from the previous state or from next state, and all nodes after the inserted node will shift to right one position
	 * @param command if not null, the new node take the data inside command
	 */
	void insertState(int d, boolean before, HmmCommand command) {
		clearPath();
		M++;

		int copyFrom;
		if (before)
			copyFrom = d - 1;
		else
			copyFrom = d;

		NodeState newmNodes[] = new NodeState[M + 1];
		NodeState newiNodes[] = new NodeState[M];
		NodeState newdNodes[] = new NodeState[M + 1];

		LineTransitionState[] newmm = new LineTransitionState[M];
		LineTransitionState[] newmi = new LineTransitionState[M];
		LineTransitionState[] newmd = new LineTransitionState[M];
		LineTransitionState[] newim = new LineTransitionState[M];
		LoopTransitionState[] newii = new LoopTransitionState[M];
		LineTransitionState[] newdm = new LineTransitionState[M];
		LineTransitionState[] newdd = new LineTransitionState[M];
		BeginTransitionState[] newbm = new BeginTransitionState[M + 1];
		EndTransitionState[] newme = new EndTransitionState[M + 1];

		for (int i = 1; i < d; i++) {
			newmNodes[i] = mNodes[i];
			newdNodes[i] = dNodes[i];
			newbm[i] = bm[i];
			newme[i] = me[i];
			if (i != M - 1) {
				newiNodes[i] = iNodes[i];
				newmm[i] = mm[i];
				newmi[i] = mi[i];
				newmd[i] = md[i];
				newim[i] = im[i];
				newii[i] = ii[i];
				newdm[i] = dm[i];
				newdd[i] = dd[i];
			}
		}

		// new node, if command is not null (pass as a reference), use the value
		// in command, other wise,
		// copy from the node next to it, if "before" is true, copy from
		// previous node, otherwise, copy from next node
		if (command != null) {
			newmNodes[d] = command.oldState.mNodes;
			newdNodes[d] = command.oldState.dNodes;
			newbm[d] = command.oldState.bm;
			newme[d] = command.oldState.me;
			if (d != M) {
				newiNodes[d] = command.oldState.iNodes;
				newmm[d] = command.oldState.mm;
				newmi[d] = command.oldState.mi;
				newmd[d] = command.oldState.md;
				newim[d] = command.oldState.im;
				newii[d] = command.oldState.ii;
				newdm[d] = command.oldState.dm;
				newdd[d] = command.oldState.dd;
			} else {
				newiNodes[d - 1] = command.oldState.iNodes;
				newmm[d - 1] = command.oldState.mm;
				newmi[d - 1] = command.oldState.mi;
				newmd[d - 1] = command.oldState.md;
				newim[d - 1] = command.oldState.im;
				newii[d - 1] = command.oldState.ii;
				newdm[d - 1] = command.oldState.dm;
				newdd[d - 1] = command.oldState.dd;
			}
		} else {
			newmNodes[d] = (NodeState) mNodes[copyFrom].clone();
			newdNodes[d] = (NodeState) dNodes[copyFrom].clone();
			newmNodes[d].label = "#" + mNodes[copyFrom].label;
			newdNodes[d].label = "#" + dNodes[copyFrom].label;

			newbm[d] = (BeginTransitionState) bm[copyFrom].clone();
			newme[d] = (EndTransitionState) me[copyFrom].clone();

			if (copyFrom == M - 1)
				copyFrom--;
			if (d != M) {
				newiNodes[d] = (NodeState) iNodes[copyFrom].clone();
				newiNodes[d].label = "#" + iNodes[copyFrom].label;

				newmm[d] = (LineTransitionState) mm[copyFrom].clone();
				newmi[d] = (LineTransitionState) mi[copyFrom].clone();
				newmd[d] = (LineTransitionState) md[copyFrom].clone();
				newim[d] = (LineTransitionState) im[copyFrom].clone();
				newii[d] = (LoopTransitionState) ii[copyFrom].clone();
				newdm[d] = (LineTransitionState) dm[copyFrom].clone();
				newdd[d] = (LineTransitionState) dd[copyFrom].clone();
			}

			if (before)
				copyFrom = d - 1;
			else
				copyFrom = d;

			if (d == M && before) {
				newiNodes[M - 1] = (NodeState) iNodes[M - 2].clone();
				newiNodes[M - 1].label = "#" + iNodes[M - 2].label;

				newmm[M - 1] = (LineTransitionState) mm[M - 2].clone();
				newmi[M - 1] = (LineTransitionState) mi[M - 2].clone();
				newmd[M - 1] = (LineTransitionState) md[M - 2].clone();
				newim[M - 1] = (LineTransitionState) im[M - 2].clone();
				newii[M - 1] = (LoopTransitionState) ii[M - 2].clone();
				newdm[M - 1] = (LineTransitionState) dm[M - 2].clone();
				newdd[M - 1] = (LineTransitionState) dd[M - 2].clone();
			}
		}
		for (int i = d + 1; i <= M; i++) {
			newmNodes[i] = mNodes[i - 1];
			newdNodes[i] = dNodes[i - 1];
			newbm[i] = bm[i - 1];
			newme[i] = me[i - 1];
			if (i != M) {
				newiNodes[i] = iNodes[i - 1];
				newmm[i] = mm[i - 1];
				newmi[i] = mi[i - 1];
				newmd[i] = md[i - 1];
				newim[i] = im[i - 1];
				newii[i] = ii[i - 1];
				newdm[i] = dm[i - 1];
				newdd[i] = dd[i - 1];
			}
		}
		mNodes = newmNodes;
		iNodes = newiNodes;
		dNodes = newdNodes;

		mm = newmm;
		mi = newmi;
		md = newmd;
		im = newim;
		ii = newii;
		dm = newdm;
		dd = newdd;
		bm = newbm;
		me = newme;

		this.d.width += 40;

		// if new insertion (not recover from an old delete operation), move all
		// nodes to make place for this new node
		if (command == null) {
			if (before) {
				moveNode(mNodes[d], mNodes[d].position.x + size * 4,
						mNodes[d].position.y);
				moveNode(dNodes[d], dNodes[d].position.x + size * 4,
						dNodes[d].position.y);
				if (copyFrom == M - 1) {
					moveNode(iNodes[d - 1],
							iNodes[d - 1].position.x + size * 4,
							iNodes[d - 1].position.y);
				} else {
					moveNode(iNodes[d], iNodes[d].position.x + size * 4,
							iNodes[d].position.y);
				}
			} else {
				if (copyFrom == M - 1) {
					moveNode(iNodes[d], iNodes[d].position.x + size * 4,
							iNodes[d].position.y);
				}
			}
			for (int i = d + 1; i <= M; i++) {
				moveNode(mNodes[i], mNodes[i].position.x + size * 4,
						mNodes[i].position.y);
				moveNode(dNodes[i], dNodes[i].position.x + size * 4,
						dNodes[i].position.y);
				if (i != M) {
					moveNode(iNodes[i], iNodes[i].position.x + size * 4,
							iNodes[i].position.y);
				}
			}
			moveNode(eNode, eNode.position.x + size * 4, eNode.position.y);
			moveNode(cNode, cNode.position.x + size * 4, cNode.position.y);
			moveNode(tNode, tNode.position.x + size * 4, tNode.position.y);
		}
		makeConnect();
	}
}
