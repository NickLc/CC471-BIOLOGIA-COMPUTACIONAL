package edu.ucf.cs.hmm.visualizer;

import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;

/**
 * Layout for null model, yes, only two node
 */
public class NullModelLayout {
	int size;

	int precision = 2;

	double[] nullModel;

	NodeState gNode;

	NodeState fNode;

	LineTransitionState gf;

	LoopTransitionState gg;

	Dimension d = new Dimension();

	public NullModelLayout(Hmm hmm, int _size) {
		size = _size;

		gNode = new NodeState(NodeState.TYPE_LOOP);
		fNode = new NodeState(NodeState.TYPE_SPECIAL);

		gf = new LineTransitionState();
		gg = new LoopTransitionState();

		gNode.label = "G";
		fNode.label = "F";

		gg.p = hmm.p1;
		gf.p = 1 - hmm.p1;

		defaultLayout();
		generateShape(_size);
	}

	Dimension getDimension() {
		return d;
	}

	void calculateDimension() {
		d.width = size * 8;
		d.height = size * 5;
	}

	/**
	 * rearrange nodes to default position
	 *
	 */
	void defaultLayout() {
		calculateDimension();
		gNode.position.x = size * 2;
		gNode.position.y = size * 3;
		fNode.position.x = size * 6;
		fNode.position.y = size * 3;

		gg.center.x = gNode.position.x;
		gg.center.y = gNode.position.y;

		gf.start.x = gNode.position.x + size;
		gf.start.y = gNode.position.y;
		gf.end.x = fNode.position.x - size;
		gf.end.y = fNode.position.y;
	}

	/**
	 * Respond to zoom operation
	 * @param factor new size = old size * factor
	 */
	void resize(double factor) {
		size = (int) Math.ceil((double) size * factor);
		calculateDimension();
		gNode.position.x *= factor;
		gNode.position.y *= factor;
		fNode.position.x *= factor;
		fNode.position.y *= factor;
		gg.center.x *= factor;
		gg.center.y *= factor;
		gf.start.x *= factor;
		gf.start.y *= factor;
		gf.end.x *= factor;
		gf.end.y *= factor;
		generateShape(size);
	}

	/**
	 * Generate shape of all state and attached lines and arrows
	 * @param _size
	 */
	void generateShape(int _size) {
		size = _size;
		gNode.generateShape(size);
		fNode.generateShape(size);
		gg.generateShape(size);
		gf.generateShape(size);
	}
	/**
	 * copy the probability of hmm model, will use this to determine thickness of the line and label on the line
	 * @param hmm
	 */
	void setHmmValue(Hmm hmm) {
		gf.p = hmm.p1;
		gg.p = 1 - hmm.p1;
	}
	/**
	 * Draw the hmm layout
	 * @param g2d graphics context
	 * @param clip clip of the displayable region, clip for time saving
	 */
	void drawWithLabel(Graphics2D g2d, Rectangle2D clip) {
		gg.draw(g2d, clip);
		gf.draw(g2d, clip);
		gNode.draw(g2d, clip);
		fNode.draw(g2d, clip);
		gg.drawText(g2d, clip, precision);
		gf.drawText(g2d, clip, precision);
		gNode.drawText(g2d, clip);
		fNode.drawText(g2d, clip);
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
		gg.center.x = gNode.position.x;
		gg.center.y = gNode.position.y;
		gf.start.x = gNode.position.x + size;
		gf.start.y = gNode.position.y;
		gf.end.x = fNode.position.x - size;
		gf.end.y = fNode.position.y;

		generateShape(size);
	}
}
