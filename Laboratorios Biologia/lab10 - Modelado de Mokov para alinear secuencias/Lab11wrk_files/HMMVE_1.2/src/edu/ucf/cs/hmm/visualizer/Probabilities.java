package edu.ucf.cs.hmm.visualizer;

public class Probabilities {
	String[] names;

	double[] probs;

	boolean[] modified;

	int offset;

	int size;

	public Probabilities(String[] alphabets) {
		size = alphabets.length;
		names = alphabets;
		modified = new boolean[size];
		for (int i = 0; i < size; i++)
			modified[i] = false;
	}

	public Probabilities(int _size) {
		size = _size;
		names = new String[size];
		probs = new double[size];
		offset = 0;
		modified = new boolean[size];
		for (int i = 0; i < size; i++)
			modified[i] = false;
	}

	public void set(double[] p) {
		probs = p.clone();
	}

	public void add(String s, double p) {
		names[offset] = s;
		probs[offset] = p;
		offset++;
	}

	public void normalize() {
		Hmm.FNorm(probs, size);
	}

	public double sum() {
		double sum = 0;
		for (int i = 0; i < size; i++) {
			sum += probs[i];
		}
		return sum;
	}
}
