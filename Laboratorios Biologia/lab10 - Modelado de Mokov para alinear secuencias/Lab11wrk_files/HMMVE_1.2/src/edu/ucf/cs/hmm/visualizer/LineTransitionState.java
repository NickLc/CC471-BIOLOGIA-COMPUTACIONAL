package edu.ucf.cs.hmm.visualizer;

import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Shape;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;

/**
 * Transition which appears with one line and an arrow
 *
 */
public class LineTransitionState extends TransitionState {
	Point start;

	Point end;

	Shape lineShape;

	Shape arrowShape;

	public LineTransitionState() {
		super();
		start = new Point();
		end = new Point();
	}

	/**
	 * Generate line, array shape
	 * @param _size size of shape
	 */
	void generateShape(int _size) {
		size = _size;
		lineShape = new Line2D.Float(start.x, start.y, end.x, end.y);
		arrowShape = generateArrow(start.x, start.y, end.x, end.y);
	}
	/**
	 * draw lines and arrows
	 * @param g2d graphics context
	 * @param clip clip of the displayable region, clip for efficiency
	 */
	public void draw(Graphics2D g2d, Rectangle2D clip) {
		if (show) {
			setPen(g2d);
			if (clip == null || lineShape.intersects(clip))
				g2d.draw(lineShape);
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
			String label = String.format("%." + precision + "f", p);
			int txtWidth = metrics.stringWidth(label);
			int txtHeight = metrics.getHeight();
			switch (labelPosition)

			{
			case POSITION_DEFAULT:
				g2d.drawString(label,
						(int) ((start.x + end.x) / 2 - txtWidth / 2),
						(int) ((start.y + end.y) / 2 - size * 0.2));
				break;
			case POSITION_START:
				g2d.drawString(label,
						(int) (start.x + (end.x - start.x) / 4 - txtWidth / 4),
						(int) (start.y + (end.y - start.y) / 4));
				break;
			case POSITION_END:
				g2d.drawString(label, (int) (start.x + (end.x - start.x) * 3
						/ 4 - txtWidth * 3 / 4),
						(int) (start.y + (end.y - start.y) * 3 / 4));
				break;
			}
		}
	}

	/**
	 * Determine if a point is on the line, not used yet
	 * @param e a point
	 */
	public boolean contains(Point e) {
		if (!show || lineShape == null)
			return false;
		return (lineShape.contains(e) || arrowShape.contains(e));
	}

	protected Object clone() {
		LineTransitionState result = new LineTransitionState();
		result.p = p;
		result.inPath = inPath;
		result.show = show;
		result.showP = showP;
		result.labelPosition = labelPosition;
		result.size = size;
		result.start = (Point) start.clone();
		result.end = (Point) end.clone();
		result.lineShape = null;
		result.arrowShape = null;
		return result;
	}
}
