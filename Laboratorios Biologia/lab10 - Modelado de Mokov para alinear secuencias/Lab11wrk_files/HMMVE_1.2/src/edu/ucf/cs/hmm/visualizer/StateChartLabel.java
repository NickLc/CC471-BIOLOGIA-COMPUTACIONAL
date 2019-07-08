package edu.ucf.cs.hmm.visualizer;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;

import javax.swing.JLabel;

/**
 * Pie chart for emission probability
 */
public class StateChartLabel extends JLabel {
	private int alphabetSize;

	private String[] alphabets;

	private double[] state;

	private String name;

	private Color[] colors = { new Color(150, 255, 255),
			new Color(255, 255, 150), new Color(255, 150, 255),
			new Color(150, 150, 150), new Color(150, 150, 255),
			new Color(255, 150, 150), new Color(150, 255, 150),
			new Color(75, 255, 255), new Color(255, 75, 255),
			new Color(255, 255, 75), new Color(75, 75, 255),
			new Color(255, 75, 75), new Color(75, 255, 75),
			new Color(200, 255, 255), new Color(255, 255, 200),
			new Color(255, 200, 255), new Color(200, 200, 255),
			new Color(200, 255, 200), new Color(255, 200, 200),
			new Color(200, 200, 200) };

	public void setAlphabet(String[] _alphabets) {
		alphabets = _alphabets;
		alphabetSize = alphabets.length;
	}

	public void setState(String _name, double[] _state) {
		name = _name;
		state = _state;
	}

	static Point cart(Point c, double radius, double deg) {
		double d = deg * Math.PI / 180.0;
		return new Point((int) (c.x + radius * Math.cos(d)), c.y
				+ (int) (radius * Math.sin(d)));
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		Graphics2D g2d = (Graphics2D) g;

		if (state != null) {
			g2d.setColor(Color.BLUE);
			g2d.drawString(name, 4, 14);

			int width = getSize().width - 20;
			int height = getSize().height - 20;
			Point center = new Point(getSize().width / 2, getSize().height / 2);
			int x = 10;
			int y = 10;
			int radius = width / 2;
			double startAngle = 90;

			for (int i = 0; i < alphabetSize; i++) {
				g2d.setColor(ColorManager.getColor(alphabets[i]));
				g2d.fillArc(x, y, width, height, (int) startAngle,
						(int) (state[i] * 360));
				g2d.setColor(Color.BLACK);
				startAngle += state[i] * 360;
			}
			startAngle = 90;
			for (int i = 0; i < alphabetSize; i++) {
				if (state[i] > 0.08) {
					int angle = (int) startAngle + (int) ((state[i] * 360) / 2);
					g2d.drawString(String.format("%.2f", state[i]), cart(
							center, radius / 2, -angle).x, cart(center,
							radius / 2, -angle).y);
					g2d.drawString(alphabets[i],
							cart(center, radius, -angle).x, cart(center,
									radius, -angle).y);
				}
				startAngle += state[i] * 360;
			}
		}
	}
}
