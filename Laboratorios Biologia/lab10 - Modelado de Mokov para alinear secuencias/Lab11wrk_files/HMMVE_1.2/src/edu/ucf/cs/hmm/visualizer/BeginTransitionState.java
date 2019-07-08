package edu.ucf.cs.hmm.visualizer;

import java.awt.BasicStroke;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Shape;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;

/**
 * Lines from B state to M state
 *
 */
public class BeginTransitionState extends TransitionState {
	Point start;

	Point point1;

	Point end;

	Shape line1Shape;

	Shape line2Shape;

	Shape arrowShape;

	public BeginTransitionState() {
		super();
		start = new Point();
		point1 = new Point();
		end = new Point();
	}

	/**
	 * Generate line, array shape
	 * @param _size size of shape
	 */
	void generateShape(int _size) {
		size = _size;
		line1Shape = new Line2D.Float(start.x, start.y,
				(int) (end.x - size * 1.5), start.y);
		line2Shape = new Line2D.Float((int) (end.x - size * 1.5), start.y,
				end.x, end.y);
		arrowShape = generateArrow((int) (end.x - size * 1.5), start.y, end.x,
				end.y);
	}

	/**
	 * draw lines and arrows
	 * @param g2d graphics context
	 * @param clip clip of the displayable region, clip for time saving
	 */
	public void draw(Graphics2D g2d, Rectangle2D clip) {
		if (show) {
			if (!inPath) {
				g2d.setColor(COLOR_BE);
				g2d.setStroke(new BasicStroke());
			} else {
				g2d.setColor(COLOR_INPATH);
				g2d.setStroke(new BasicStroke(4));
			}
			if (clip == null || line1Shape.intersects(clip))
				g2d.draw(line1Shape);
			if (clip == null || line2Shape.intersects(clip))
				g2d.draw(line2Shape);
			if (clip == null || arrowShape.intersects(clip))
				g2d.fill(arrowShape);
		}
	}

	/**
	 * draw text on the line
	 * @param g2d graphics context
	 * @param clip clip of the displayable region, clip for efficiency
	 * @param precision precision of the probability to display
	 */
	public void drawText(Graphics2D g2d, Rectangle2D clip, int precision) {
		if (show) {
			g2d.setColor(COLOR_TEXT);
			FontMetrics metrics = g2d.getFontMetrics();
			String label = String.format("%.2f", p);
			int txtWidth = metrics.stringWidth(label);
			int txtHeight = metrics.getHeight();
			g2d.drawString(String.format("%." + precision + "f", p),
					(int) (start.x - size * 1.5), start.y - size);
		}
	}

	/**
	 * Determine if a point is on the line, not used yet
	 * @param e a point
	 */
	public boolean contains(Point e) {
		if (!show || line1Shape == null)
			return false;
		return (line1Shape.contains(e) || line2Shape.contains(e) || arrowShape
				.contains(e));
	}

	protected Object clone() {
		BeginTransitionState result = new BeginTransitionState();
		result.p = p;
		result.inPath = inPath;
		result.show = show;
		result.showP = showP;
		result.labelPosition = labelPosition;
		result.size = size;
		result.start = (Point) start.clone();
		result.point1 = (Point) point1.clone();
		result.end = (Point) end.clone();
		result.line1Shape = null;
		result.line2Shape = null;
		result.arrowShape = null;
		return result;
	}
}