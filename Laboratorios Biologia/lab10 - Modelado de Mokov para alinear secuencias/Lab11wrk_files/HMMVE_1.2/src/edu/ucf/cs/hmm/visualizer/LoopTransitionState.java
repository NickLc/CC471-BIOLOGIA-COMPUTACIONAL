package edu.ucf.cs.hmm.visualizer;

import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
/**
 * Transition which appears with one loop and an arrow
 *
 */
public class LoopTransitionState extends TransitionState {
	Point center;

	Shape circleShape;

	Shape arrowShape;

	public LoopTransitionState() {
		super();
		center = new Point();
	}
	
	/**
	 * Generate line, array shape
	 * @param _size size of shape
	 */
	void generateShape(int _size) {
		size = _size;
		circleShape = new Ellipse2D.Float(center.x - size / 2, center.y - size
				* 2, size, size * 2);
		Point highend = new Point(center.x - size / 2 - 4, center.y - size / 2
				- 8);
		Point lowend = new Point(center.x - size / 2 + 4, center.y - size / 2
				- 8);
		arrowShape = new Polygon(new int[] { center.x - size / 2, highend.x,
				lowend.x }, new int[] { center.y - size / 2, highend.y,
				lowend.y }, 3);
	}

	/**
	 * draw lines and arrows
	 * @param g2d graphics context
	 * @param clip clip of the displayable region, clip for efficiency
	 */
	public void draw(Graphics2D g2d, Rectangle2D clip) {
		if (show) {
			setPen(g2d);
			if (clip == null || circleShape.intersects(clip))
				g2d.draw(circleShape);
			if (clip == null || arrowShape.intersects(clip))
				g2d.fill(arrowShape);
		}
	}

	/**
	 * draw text on the line
	 * @param g2d graphics context
	 * @param clip clip of the displayable region, clip for time saving
	 * @param precision precision of the probability to display
	 */
	public void drawText(Graphics2D g2d, Rectangle2D clip, int precision) {
		if (show) {
			g2d.setColor(COLOR_TEXT);
			FontMetrics metrics = g2d.getFontMetrics();
			String label = String.format("%.2f", p);
			int txtWidth = metrics.stringWidth(label);
			int txtHeight = metrics.getHeight();
			g2d.drawString(String.format("%." + precision + "f", p), center.x
					- size / 2, (int) (center.y - size * 1.2));
		}
	}

	/**
	 * Determine if a point is on the line, not used yet
	 * @param e a point
	 */
	public boolean contains(Point e) {
		if (!show || circleShape == null)
			return false;
		return (circleShape.contains(e) || arrowShape.contains(e));
	}

	protected Object clone() {
		LoopTransitionState result = new LoopTransitionState();
		result.p = p;
		result.inPath = inPath;
		result.show = show;
		result.showP = showP;
		result.labelPosition = labelPosition;
		result.size = size;
		result.center = (Point) center.clone();
		result.circleShape = null;
		result.arrowShape = null;
		return result;
	}
}
