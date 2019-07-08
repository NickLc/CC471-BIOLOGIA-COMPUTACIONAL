package edu.ucf.cs.hmm.visualizer;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import javax.imageio.ImageIO;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

/**
 * Panel to show Hmm Logo
 */
public class LogoPanel extends JPanel {
	Hmm hmm;

	int size;

	BufferedImage image;

	int totalWidth;

	JScrollPane container;

	boolean bitmapMode = false;

	File tempLogoFile = null;

	public LogoPanel(Hmm _hmm, int _size, JScrollPane _container) {
		hmm = _hmm;
		size = _size;
		container = _container;
	}

	public void resize(double factor) {
		size = (int) Math.round(size * factor);
		generateLogo();
	}

	/**
	 * Generate logo image to a temporary file (I tried memory image, will fail in some cases for unknown reasons)
	 *
	 */
	public void generateLogo() {
		double[] height = getHeights(hmm);
		int[] width = new int[hmm.M + 1];

		// h[][0]--h(M)
		// h[][1]--h(D)
		// h[][2]--h(I)
		double[][] h = new double[hmm.M + 1][3];

		// c[][0]--c(M)
		// c[][1]--c(I)
		double[][] c = new double[hmm.M + 1][2];

		// calculate h (hit) and c (contribution)
		h[1][0] = hmm.begin[1];
		h[1][1] = 1 - h[1][0];
		h[1][2] = h[1][0] * hmm.t[1][Hmm.TMI];
		c[1][0] = h[1][0];
		c[1][1] = h[1][2] / (1 - hmm.t[1][Hmm.TII]);

		for (int i = 2; i <= hmm.M; i++) {
			double mu = hmm.t[i - 1][Hmm.TIM] / (1 - hmm.t[i - 1][Hmm.TII]);
			h[i][0] = h[i - 1][0]
					* (hmm.t[i - 1][Hmm.TMM] + hmm.t[i - 1][Hmm.TMI] * mu)
					+ (1 - h[i - 1][0]) * hmm.t[i - 1][Hmm.TDM];
			if (i != hmm.M) {
				h[i][1] = 1 - h[i][0];
				h[i][2] = h[i][0] * hmm.t[i][Hmm.TMI];
			}
			c[i][0] = h[i][0];
			if (i != hmm.M)
				c[i][1] = h[i][2] / (1 - hmm.t[i][Hmm.TII]);
		}

		width[0] = size;
		for (int i = 1; i <= hmm.M; i++) {
			width[i] = width[i - 1] + (int) Math.round(size * c[i][0]);
			width[i] += (int) Math.round(size * c[i][1]) == 0 ? 1 : (int) Math
					.round(size * c[i][1]);
		}
		totalWidth = width[hmm.M] + size;

		createImage(width, height, c);
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		Graphics2D g2d = (Graphics2D) g;
		Rectangle2D.Float clip = new Rectangle2D.Float();
		clip.x = container.getViewport().getViewPosition().x;
		clip.y = container.getViewport().getViewPosition().y;
		clip.width = container.getViewport().getExtentSize().width;
		clip.height = container.getViewport().getExtentSize().height;
		BufferedImage subImage;
		if (!bitmapMode)
			subImage = image.getSubimage((int) clip.x, (int) clip.y, (int) Math
					.min(clip.width, image.getWidth() - clip.x), (int) Math
					.min(clip.height, image.getHeight() - clip.y));
		else
			subImage = image;
		g2d.drawImage(subImage, null, (int) clip.x, (int) clip.y);
	}

    private void createImage(int[] width, double[] height, double[][] c)
    {

        image = new BufferedImage(totalWidth, size * 16, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, image.getWidth(), image.getHeight());

        g2d.setColor(Color.BLACK);
        FontMetrics metrics = g2d.getFontMetrics();
        int a = (int) (metrics.getAscent());
        int w;
        g2d.drawLine(size, size, size, size * 14);

        Font logoFont = new Font("Lucida", Font.BOLD, 100);
        g2d.setFont(logoFont);
        metrics = g2d.getFontMetrics();
        a = metrics.getAscent();
        AffineTransform t = new AffineTransform();

        int cumulateWidth = size;
        width[0] = size;
        int cumulateHeight;
        int x;
        int y;
        for(int i = 1; i <= hmm.M; i++)
        {
            Map<Double, String> map = new TreeMap<Double, String>();
            for(int j = 0; j < hmm.alphabetSize; j++)
                map.put(hmm.mat[i][j], hmm.alphabets[j]);
            String key;
            cumulateHeight = 0;
            int thisHeight;
            int thisWidth;
            x = cumulateWidth;
            thisWidth = (int) Math.round(size * c[i][0]);
            cumulateWidth += thisWidth;

            AffineTransform oldT = g2d.getTransform();

            for(Double value : map.keySet())
            {
                key = map.get(value);
                w = metrics.stringWidth(key);
                y = size * 14 - cumulateHeight;
                thisHeight = (int) Math.round(height[i] * size * 14 * value);
                cumulateHeight += thisHeight * 0.7;
                g2d.setColor(ColorManager.getColor(key));
                if(thisHeight == 0)
                    continue;
                t.setToTranslation(x, y);
                // from w, a scale to thisWeight, thisHeight
                t.scale((double) thisWidth / w, (double) thisHeight / a);
                g2d.setTransform(t);
                //System.out.println("key="+key);
                g2d.drawString(key, 0, 0);
            }
            g2d.setTransform(oldT);
            g2d.setColor(new Color(0xFF, 0xCC, 0xCC));
            int insertionWidth;
            if(i != hmm.M)
                insertionWidth = (int) Math.round(size * c[i][1]);
            else
                insertionWidth = 1;
            if(insertionWidth == 0)
                insertionWidth = 1;
            g2d.fillRect(cumulateWidth, size, insertionWidth, size * 13);
            cumulateWidth += insertionWidth;
            width[i] = cumulateWidth;
        }
        g2d.setFont(new Font("Courier", Font.PLAIN, 12));
        metrics = g2d.getFontMetrics();
        for(int i = 1; i <= hmm.M; i++)
        {
            g2d.setColor(Color.BLACK);
            g2d.drawLine(width[i - 1], size * 14, width[i], size * 14);
            g2d.drawLine(width[i], (int) (size * 14 + size * 0.2), width[i], (int) (size * 14 - size * 0.2));
            int wid = metrics.stringWidth("" + i);
            g2d.drawString("" + i, width[i - 1] + (width[i] - width[i - 1]) / 2 - wid / 2, 14 * size + size / 2);
        }
        try
        {
            tempLogoFile = File.createTempFile("Logo", "png");
            ImageIO.write(image, "png", tempLogoFile);
            image = ImageIO.read(tempLogoFile);
        }catch(IOException e)
        {
        }
        setPreferredSize(new Dimension(totalWidth + size * 2, size * 16));
        revalidate();
    }

    public BufferedImage getImage()
    {
        return image;
    }
    public static  double[] getHeights(Hmm hmm)
    {
		double[] height = new double[hmm.M + 1];
        double maxScore = 0;
        for(int i = 1; i <= hmm.M; i++)
        {
            int totalScore = 0;
            for(int j = 0; j < hmm.alphabetSize; j++)
                totalScore += hmm.mat[i][j] * Hmm.prob2Score(hmm.mat[i][j], hmm.nullModel[j]);
            height[i] = totalScore;
            if(totalScore > maxScore)
                maxScore = totalScore;
        }
        for(int i = 1; i <= hmm.M; i++)
            height[i] /= maxScore;
        return height;
    }
}
