package edu.ucf.cs.hmm.visualizer;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.Shape;
import java.awt.geom.Rectangle2D;

/**
 * Abstract class for all connection lines
 */
abstract public class TransitionState {
	double p;

	boolean inPath;

	boolean show;

	boolean showP;

	int labelPosition;

	int size;

	static final int POSITION_DEFAULT = 0;

	static final int POSITION_START = 1;

	static final int POSITION_END = 2;

	static final Color COLOR_LINE = Color.BLACK;

	static final Color COLOR_TEXT = Color.BLACK;

	static final Color COLOR_LIGHT_LINE = Color.LIGHT_GRAY;

	static final Color COLOR_BE = new Color(100, 255, 100);

	static final Color COLOR_INPATH = Color.BLUE;

	public TransitionState() {
		p = 0;
		inPath = false;
		show = true;
		showP = true;
		labelPosition = POSITION_DEFAULT;
	}

	abstract void generateShape(int _size);

	protected Shape generateArrow(int x1, int y1, int x2, int y2) {
		double theta = Math.atan(((double) y2 - y1) / (x2 - x1));
		Point c = null;
		if (x2 > x1) {
			c = new Point(x2 - (int) (8 * Math.cos(theta)), y2
					- (int) (8 * Math.sin(theta)));
		} else {
			c = new Point(x2 + (int) (8 * Math.cos(theta)), y2
					- (int) (8 * Math.sin(theta)));
		}
		Point highend = null;
		Point lowend = null;
		highend = new Point(c.x - (int) (4 * Math.sin(theta)), c.y
				+ (int) (4 * Math.cos(theta)));
		lowend = new Point(c.x + (int) (4 * Math.sin(theta)), c.y
				- (int) (4 * Math.cos(theta)));
		Shape arrow = new Polygon(new int[] { x2, highend.x, lowend.x },
				new int[] { y2, highend.y, lowend.y }, 3);
		return arrow;
	}

	protected void setPen(Graphics2D g2d) {
		if (inPath) {
			g2d.setStroke(new BasicStroke(4));
			g2d.setColor(COLOR_INPATH);
			return;
		}
		if (p >= 0.75) {
			g2d.setColor(COLOR_LINE);
			g2d.setStroke(new BasicStroke(3));
		} else if (p >= 0.5) {
			g2d.setColor(COLOR_LINE);
			g2d.setStroke(new BasicStroke(2));
		} else if (p >= 0.25) {
			g2d.setStroke(new BasicStroke());
			g2d.setColor(COLOR_LINE);

		} else {
			g2d.setStroke(new BasicStroke());
			g2d.setColor(COLOR_LIGHT_LINE);
		}

	}

	abstract public void draw(Graphics2D g2d, Rectangle2D clip);

	abstract public void drawText(Graphics2D g2d, Rectangle2D clip,
			int precision);

	abstract public boolean contains(Point e);

	abstract protected Object clone();
}
