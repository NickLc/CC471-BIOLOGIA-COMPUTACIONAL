package edu.ucf.cs.hmm.visualizer;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;

/**
 * Shape for all node, include matching, deletion, insertion
 * 
 */
public class NodeState {
	public NodeState(int _type) {
		type = _type;
		position = new Point();
		p = 0;
		inPath = false;
		show = true;
		showLabel = true;
		showConsensus = false;
		modified = false;
		label = "";
		consensus = "";
		shape = null;
	}

	int type;

	int size;

	Point position;

	double p;

	boolean modified;

	boolean inPath;

	boolean show;

	boolean showLabel;

	boolean showConsensus;

	String label;

	String consensus;

	Shape shape;

	static final int TYPE_UNKNOWN = -1;

	static final int TYPE_M = 1;

	static final int TYPE_I = 2;

	static final int TYPE_D = 3;

	static final int TYPE_SPECIAL = 4;

	static final int TYPE_LOOP = 5;

	static final Color COLOR_BORDER = Color.BLACK;

	static final Color COLOR_TEXT = Color.BLACK;

	static final Color COLOR_M = new Color(230, 255, 255);

	static final Color COLOR_I = new Color(255, 230, 255);

	static final Color COLOR_D = new Color(255, 255, 230);

	static final Color COLOR_SPECIAL = new Color(230, 230, 255);

	static final Color COLOR_LOOP = new Color(230, 230, 230);

	static final Color COLOR_SELECTED = Color.RED;

	static final Color COLOR_INPATH = new Color(150, 150, 255);

	static final Color[] COLORS = { null, COLOR_M, COLOR_I, COLOR_D,
			COLOR_SPECIAL, COLOR_LOOP };

	/**
	 * Generate shape of the node providing the center
	 * @param x center position of the node
	 * @param y 
	 * @return
	 */
	Shape generateShapeEx(int x, int y) {
		Shape s = null;
		switch (type) {
		case TYPE_M:
			s = new Rectangle(x - size, y - size, size * 2, size * 2);
			break;
		case TYPE_I:
			s = new Polygon(new int[] { x - size, x, x + size, x }, new int[] {
					y, y - size, y, y + size }, 4);
			break;
		case TYPE_D:
			s = new Ellipse2D.Float(x - size, y - size, size * 2, size * 2);
			break;
		case TYPE_SPECIAL:
			s = new Ellipse2D.Float(x - size, y - size, size * 2, size * 2);
			break;
		case TYPE_LOOP:
			s = new Polygon(new int[] { x - size, x, x + size, x }, new int[] {
					y, y - size, y, y + size }, 4);
			break;
		}
		return s;
	}

	/**
	 * Generate node shape, center is set before
	 * @param _size
	 */
	void generateShape(int _size) {
		size = _size;
		shape = generateShapeEx(position.x, position.y);
	}

	public void draw(Graphics2D g2d, Rectangle2D clip) {
		if (show && (clip == null || shape.intersects(clip))) {
			if (p >= 0.5 && type == TYPE_M)
				g2d.setStroke(new BasicStroke(3));
			else
				g2d.setStroke(new BasicStroke());
			if (!inPath)
				g2d.setColor(COLORS[type]);
			else
				g2d.setColor(COLOR_INPATH);
			g2d.fill(shape);
			g2d.setColor(COLOR_BORDER);
			g2d.draw(shape);
		}
	}

	/**
	 * When mouse is inside the shape, change the color of border
	 * @param g2d graphics context
	 * @param selected whether or not the node is currently in selected state (mouse in)
	 */
	public void toggleSelected(Graphics2D g2d, boolean selected) {
		if (p >= 0.5 && type == TYPE_M)
			g2d.setStroke(new BasicStroke(3));
		else
			g2d.setStroke(new BasicStroke());
		if (selected) {
			g2d.setColor(COLOR_SELECTED);
		} else
			g2d.setColor(COLOR_BORDER);
		g2d.draw(shape);
	}

	/**
	 * Draw lable of the node
	 * @param g2d graphics context
	 * @param clip  of the displayable region, clip for efficiency
	 */
	public void drawText(Graphics2D g2d, Rectangle2D clip) {
		if (show && (clip == null || shape.intersects(clip))) {
			String l;
			if (showConsensus)
				l = consensus;
			else
				l = label;
			if (modified)
				l += "*";
			FontMetrics metrics = g2d.getFontMetrics();
			int txtWidth = metrics.stringWidth(l);
			int txtHeight = metrics.getHeight();
			g2d.setColor(COLOR_TEXT);
			g2d.drawString(l, position.x - txtWidth / 2, position.y + txtHeight
					/ 2);
		}
	}
	/**
	 * Determine if a point is inside node
	 * @param e a point
	 */
	public boolean contains(Point e) {
		if (!show || shape == null)
			return false;
		return shape.contains(e);
	}

	protected Object clone() {
		NodeState result = new NodeState(type);
		result.type = type;
		result.size = size;
		result.position = (Point) position.clone();
		result.p = p;
		result.modified = modified;
		result.inPath = inPath;
		result.show = show;
		result.showLabel = showLabel;
		result.showConsensus = showConsensus;
		result.label = label;
		result.consensus = consensus;
		result.shape = null;
		return result;
	}
}
