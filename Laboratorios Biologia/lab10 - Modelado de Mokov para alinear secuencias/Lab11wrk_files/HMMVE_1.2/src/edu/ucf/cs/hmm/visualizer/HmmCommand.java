package edu.ucf.cs.hmm.visualizer;

/**
 * Structure to support undo, redo
 */
class HmmNodeStateInCommand {
	double t[];

	double mat[];

	double ins[];

	double tbd1;

	double begin;

	double end;

	boolean modified[];

	String rf;

	String cs;

	double nullModel[];

	double p1;

	double xt[][];

	double beginArray[];

	NodeState mNodes;

	NodeState dNodes;

	NodeState iNodes;

	LineTransitionState mm;

	LineTransitionState mi;

	LineTransitionState md;

	LineTransitionState im;

	LoopTransitionState ii;

	LineTransitionState dm;

	LineTransitionState dd;

	BeginTransitionState bm;

	EndTransitionState me;

	HmmNodeStateInCommand(Hmm _hmm) {
		t = new double[7];
		modified = new boolean[3];
		xt = new double[4][2];
		beginArray = new double[_hmm.M + 1];
	}
}

public class HmmCommand {
	Hmm hmm;

	HmmLayout layout;

	NullModelLayout nullModelLayout;

	enum Target {
		ProfileHmm, NullHmm, NNode, BNode, ENode, CNode, JNode
	}

	enum CommandType {
		Modify, Insert, Delete
	}

	Target target;

	CommandType commandType;

	int position;

	boolean before; // for insert command only

	HmmNodeStateInCommand oldState;

	HmmNodeStateInCommand newState;

	HmmCommand(Hmm _hmm, HmmLayout _layout, NullModelLayout _nullModelLayout) {
		hmm = _hmm;
		layout = _layout;
		nullModelLayout = _nullModelLayout;
		oldState = new HmmNodeStateInCommand(_hmm);
		newState = new HmmNodeStateInCommand(_hmm);
	}

	/**
	 * Undo this command
	 *
	 */
	void undo() {
		if (target == Target.ProfileHmm) {
			if (commandType == CommandType.Modify) {
				hmm.mat[position] = oldState.mat;
				hmm.begin[position] = oldState.begin;
				hmm.end[position] = oldState.end;
				layout.mNodes[position].modified = oldState.modified[0];
				layout.dNodes[position].modified = oldState.modified[2];
				if (position != hmm.M) {
					layout.iNodes[position].modified = oldState.modified[1];
					hmm.t[position] = oldState.t;
					hmm.ins[position] = oldState.ins;
				}
			} else if (commandType == CommandType.Insert) {
				hmm.deleteState(position, null);
				layout.deleteState(position, null);
			} else if (commandType == CommandType.Delete) {
				hmm.insertState(position, before, this);
				layout.insertState(position, before, this);
			}
		} else if (target == Target.NullHmm) {
			if (commandType == CommandType.Modify) {
				hmm.nullModel = oldState.nullModel;
				hmm.p1 = oldState.p1;
				nullModelLayout.gNode.modified = oldState.modified[0];
			}
		} else if (target == Target.NNode) {
			hmm.xt[Hmm.XTN][Hmm.MOVE] = oldState.xt[Hmm.XTN][Hmm.MOVE];
			hmm.xt[Hmm.XTN][Hmm.LOOP] = oldState.xt[Hmm.XTN][Hmm.LOOP];
			layout.nNode.modified = oldState.modified[0];
		} else if (target == Target.BNode) {
			hmm.tbd1 = oldState.tbd1;
			hmm.begin = oldState.beginArray;
			layout.bNode.modified = oldState.modified[0];
		} else if (target == Target.ENode) {
			hmm.xt[Hmm.XTE][Hmm.MOVE] = oldState.xt[Hmm.XTE][Hmm.MOVE];
			hmm.xt[Hmm.XTE][Hmm.LOOP] = oldState.xt[Hmm.XTE][Hmm.LOOP];
			layout.eNode.modified = oldState.modified[0];
		} else if (target == Target.CNode) {
			hmm.xt[Hmm.XTC][Hmm.MOVE] = oldState.xt[Hmm.XTC][Hmm.MOVE];
			hmm.xt[Hmm.XTC][Hmm.LOOP] = oldState.xt[Hmm.XTC][Hmm.LOOP];
			layout.cNode.modified = oldState.modified[0];
		} else if (target == Target.JNode) {
			hmm.xt[Hmm.XTJ][Hmm.MOVE] = oldState.xt[Hmm.XTJ][Hmm.MOVE];
			hmm.xt[Hmm.XTJ][Hmm.LOOP] = oldState.xt[Hmm.XTJ][Hmm.LOOP];
			layout.jNode.modified = oldState.modified[0];
		}
	}

	/**
	 * redo this command
	 *
	 */
	void redo() {
		if (target == Target.ProfileHmm) {
			if (commandType == CommandType.Modify) {
				hmm.mat[position] = newState.mat;
				hmm.begin[position] = newState.begin;
				hmm.end[position] = newState.end;
				layout.mNodes[position].modified = newState.modified[0];
				layout.dNodes[position].modified = newState.modified[2];
				if (position != hmm.M) {
					layout.iNodes[position].modified = newState.modified[1];
					hmm.t[position] = newState.t;
					hmm.ins[position] = newState.ins;
				}
			} else if (commandType == CommandType.Insert) {
				hmm.insertState(position, before, null);
				layout.insertState(position, before, null);
			} else if (commandType == CommandType.Delete) {
				hmm.deleteState(position, null);
				layout.deleteState(position, null);
			}
		} else if (target == Target.NullHmm) {
			if (commandType == CommandType.Modify) {
				hmm.nullModel = newState.nullModel;
				hmm.p1 = newState.p1;
				nullModelLayout.gNode.modified = newState.modified[0];
			}
		} else if (target == Target.NNode) {
			hmm.xt[Hmm.XTN][Hmm.MOVE] = newState.xt[Hmm.XTN][Hmm.MOVE];
			hmm.xt[Hmm.XTN][Hmm.LOOP] = newState.xt[Hmm.XTN][Hmm.LOOP];
			layout.nNode.modified = newState.modified[0];
		} else if (target == Target.BNode) {
			hmm.tbd1 = newState.tbd1;
			hmm.begin = newState.beginArray;
			layout.bNode.modified = newState.modified[0];
		} else if (target == Target.ENode) {
			hmm.xt[Hmm.XTE][Hmm.MOVE] = newState.xt[Hmm.XTE][Hmm.MOVE];
			hmm.xt[Hmm.XTE][Hmm.LOOP] = newState.xt[Hmm.XTE][Hmm.LOOP];
			layout.eNode.modified = newState.modified[0];
		} else if (target == Target.CNode) {
			hmm.xt[Hmm.XTC][Hmm.MOVE] = newState.xt[Hmm.XTC][Hmm.MOVE];
			hmm.xt[Hmm.XTC][Hmm.LOOP] = newState.xt[Hmm.XTC][Hmm.LOOP];
			layout.cNode.modified = newState.modified[0];
		} else if (target == Target.JNode) {
			hmm.xt[Hmm.XTJ][Hmm.MOVE] = newState.xt[Hmm.XTJ][Hmm.MOVE];
			hmm.xt[Hmm.XTJ][Hmm.LOOP] = newState.xt[Hmm.XTJ][Hmm.LOOP];
			layout.jNode.modified = newState.modified[0];
		}
	}
}
