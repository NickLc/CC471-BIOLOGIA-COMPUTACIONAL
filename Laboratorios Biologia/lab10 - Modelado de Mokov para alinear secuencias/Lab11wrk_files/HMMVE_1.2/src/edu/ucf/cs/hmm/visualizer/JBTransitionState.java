package edu.ucf.cs.hmm.visualizer;

import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Shape;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
/**
 * Lines from J state to B state
 *
 */
public class JBTransitionState extends TransitionState {
	Point start;

	Point point1;

	Point point2;

	Point end;

	Shape line1Shape;

	Shape line2Shape;

	Shape line3Shape;

	Shape arrowShape;

	public JBTransitionState() {
		super();
		start = new Point();
		point1 = new Point();
		point2 = new Point();
		end = new Point();
	}
	/**
	 * Generate line, array shape
	 * @param _size size of shape
	 */
	void generateShape(int _size) {
		size = _size;
		line1Shape = new Line2D.Float(start.x, start.y, end.x - size, start.y);
		line2Shape = new Line2D.Float(end.x - size, start.y, end.x - size,
				end.y + size);
		line3Shape = new Line2D.Float(end.x - size, end.y + size, end.x, end.y);
		arrowShape = generateArrow(end.x - size, end.y + size, end.x, end.y);
	}
	/**
	 * draw lines and arrows
	 * @param g2d graphics context
	 * @param clip clip of the displayable region, clip for efficiency
	 */
	public void draw(Graphics2D g2d, Rectangle2D clip) {
		if (show) {
			setPen(g2d);
			if (clip == null || line1Shape.intersects(clip))
				g2d.draw(line1Shape);
			if (clip == null || line2Shape.intersects(clip))
				g2d.draw(line2Shape);
			if (clip == null || line3Shape.intersects(clip))
				g2d.draw(line3Shape);
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
			g2d.drawString(String.format("%." + precision + "f", p), end.x
					+ size * 2, start.y - size / 2);
		}
	}

	/**
	 * Determine if a point is on the line, not used yet
	 * @param e a point
	 */
	public boolean contains(Point e) {
		if (!show || line1Shape == null)
			return false;
		return (line1Shape.contains(e) || line2Shape.contains(e)
				|| line3Shape.contains(e) || arrowShape.contains(e));
	}

	protected Object clone() {
		JBTransitionState result = new JBTransitionState();
		result.p = p;
		result.inPath = inPath;
		result.show = show;
		result.showP = showP;
		result.labelPosition = labelPosition;
		result.size = size;
		result.start = (Point) start.clone();
		result.point1 = (Point) point1.clone();
		result.point2 = (Point) point2.clone();
		result.end = (Point) end.clone();
		result.line1Shape = null;
		result.line2Shape = null;
		result.line3Shape = null;
		result.arrowShape = null;
		return result;
	}
}
