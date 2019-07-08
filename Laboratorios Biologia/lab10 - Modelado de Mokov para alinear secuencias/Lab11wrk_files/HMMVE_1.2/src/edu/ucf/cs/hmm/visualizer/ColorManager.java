package edu.ucf.cs.hmm.visualizer;

import java.awt.Color;

/**
 * 
 * Manage colors of amino acid
 *
 */
public class ColorManager {
	static Color[] colors = { new Color(0xFF, 0x99, 0x66),
			new Color(0x00, 0x99, 0x99), new Color(0xFF, 0x00, 0x00),
			new Color(0xCC, 0x00, 0x33), new Color(0x00, 0xFF, 0x00),
			new Color(0xFF, 0xFF, 0x00), new Color(0x66, 0x00, 0x33),
			new Color(0xCC, 0x99, 0x33), new Color(0x66, 0x33, 0x00),
			new Color(0xFF, 0x99, 0x33), new Color(0xCC, 0x99, 0xCC),
			new Color(0x33, 0x66, 0x66), new Color(0x00, 0x99, 0xFF),
			new Color(0x66, 0x66, 0xCC), new Color(0x99, 0x00, 0x00),
			new Color(0x00, 0x00, 0xFF), new Color(0x00, 0xFF, 0xFF),
			new Color(0xFF, 0xCC, 0x33), new Color(0x66, 0xCC, 0x66),
			new Color(0x00, 0x66, 0x00), new Color(0xCC, 0xCC, 0xCC) };

	static Color getColor(String key) {
		switch (key.charAt(0)) {
		case 'A':
			return colors[0];
		case 'C':
			return colors[1];
		case 'D':
			return colors[2];
		case 'E':
			return colors[3];
		case 'F':
			return colors[4];
		case 'G':
			return colors[5];
		case 'H':
			return colors[6];
		case 'I':
			return colors[7];
		case 'K':
			return colors[8];
		case 'L':
			return colors[9];
		case 'M':
			return colors[10];
		case 'N':
			return colors[11];
		case 'P':
			return colors[12];
		case 'Q':
			return colors[13];
		case 'R':
			return colors[14];
		case 'S':
			return colors[15];
		case 'T':
			return colors[16];
		case 'V':
			return colors[17];
		case 'W':
			return colors[18];
		case 'Y':
			return colors[19];
		default:
			return colors[20];
		}
	}
}
